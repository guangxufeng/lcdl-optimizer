"""两阶段 LCDL 模型 v2：理想储能、固定交换、情景自适应模式。

Public API: Case, SolverOptions, compute_unified_directrix, compute_directrix,
solve_scenario, run_two_stage, export_result, validate_feedback, dispatch_feedback.

MW / MWh / hour / currency; xi, L, d are dimensionless energy shares.
All recourse bounds are explicit rows for the exact continuous-equivalent feasibility oracle.
No sampled result is ever labelled robust. See docs/algorithm.md for dual derivations.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import copy
import hashlib
import itertools
from pathlib import Path
import time

import gurobipy as gp
from gurobipy import GRB
import numpy as np
from scipy import sparse as sp


class ModelDataError(ValueError):
    pass


@dataclass
class SolverOptions:
    time_limit: float = 60.0  # per master / oracle, not a total wall-clock limit
    max_iterations: int = 30
    feasibility_tol: float = 1e-6  # sum of normalized phase-I violations
    cost_abs_tol: float = 0.01
    cost_rel_tol: float = 1e-3
    upgrade_tol: float = 1e-5
    threads: int = 4
    seed: int = 20260916
    output_flag: bool = False
    log_dir: str | None = None
    two_period_vertex_limit: int = 6  # exact interval reduction ONLY when T=2; 0 disables it

    def __post_init__(self):
        if self.time_limit <= 0 or self.max_iterations < 1 or self.threads < 1:
            raise ValueError("求解时间、迭代次数和线程数必须为正。")
        if min(self.feasibility_tol, self.cost_abs_tol, self.cost_rel_tol, self.upgrade_tol) <= 0:
            raise ValueError("容差必须为正。")
        if self.two_period_vertex_limit<0 or self.two_period_vertex_limit>12:
            raise ValueError("two_period_vertex_limit must lie in 0..12")


def _array(value, shape=None, name="array"):
    a = np.asarray(value, dtype=float)
    if shape is not None and a.shape != shape:
        raise ModelDataError(f"{name}: expected {shape}, got {a.shape}")
    if not np.isfinite(a).all():
        raise ModelDataError(f"{name}: NaN/Inf is not allowed")
    return a


class Case:
    """Load any radial network through the documented schema, with arbitrary bus IDs/order."""
    @classmethod
    def from_json(cls, path):
        return cls(json.loads(Path(path).read_text(encoding="utf-8")))

    def __init__(self, data):
        data = copy.deepcopy(data)
        self.data = data
        if data.get("schema_version") not in (1, 2):
            raise ModelDataError("Unsupported schema_version")
        self.name = str(data["name"])
        self.ids = list(data["bus_ids"])
        self.n = len(self.ids)
        if self.n < 2 or len(set(self.ids)) != self.n:
            raise ModelDataError("bus_ids must be unique, with at least two buses")
        self.index = {b: i for i, b in enumerate(self.ids)}
        if data["root_bus"] not in self.index:
            raise ModelDataError("root_bus is not present in bus_ids")
        self.root = self.index[data["root_bus"]]
        self.dt = float(data["dt_hours"])
        self.kv = float(data["base_kv"])
        if not np.isfinite([self.dt, self.kv]).all() or min(self.dt, self.kv) <= 0:
            raise ModelDataError("dt_hours/base_kv must be positive and finite")
        self.pre = _array(data["dr_pre_mw"], name="dr_pre_mw")
        if self.pre.ndim != 2:
            raise ModelDataError("dr_pre_mw must be K x T")
        self.k, self.t = self.pre.shape
        if self.k < 1 or self.t < 2 or (self.pre < 0).any():
            raise ModelDataError("Need at least one DR user and two time steps; load must be nonnegative")
        self.energy = self.pre.sum(axis=1) * self.dt
        if (self.energy <= 0).any():
            raise ModelDataError("Each DR user must have positive cycle energy")
        self.p = _array(data["rigid_p_mw"], (self.n, self.t), "rigid_p_mw")
        self.q = _array(data["fixed_q_mvar"], (self.n, self.t), "fixed_q_mvar")
        self.renew = _array(data["renewable_p_mw"], (self.n, self.t), "renewable_p_mw")
        self.w = _array(data["dr_allocation"], (self.n, self.k), "dr_allocation")
        if (self.p < 0).any() or (self.renew < 0).any() or (self.w < 0).any():
            raise ModelDataError("Load, renewable and allocation must be nonnegative")
        if not np.allclose(self.w.sum(axis=0), 1, atol=1e-10, rtol=0):
            raise ModelDataError("Every dr_allocation column must sum to 1")
        if any(np.any(a[self.root]) for a in (self.p, self.q, self.renew, self.w)):
            raise ModelDataError("Root is an external exchange bus; connect local resources to a downstream bus")
        self.stores = list(data["storage"])
        self.s = len(self.stores)
        self.storage_map = np.zeros((self.n, self.s))
        for s, b in enumerate(self.stores):
            b.setdefault("eta_charge", 1.0)
            b.setdefault("eta_discharge", 1.0)
            if b["bus"] not in self.index:
                raise ModelDataError("Storage references an unknown bus")
            idx = self.index[b["bus"]]
            if idx == self.root:
                raise ModelDataError("Storage must connect to a non-root bus")
            self.storage_map[idx, s] = 1
            vals = [b[key] for key in ("p_charge_mw", "p_discharge_mw", "e_min_mwh", "e_max_mwh",
                                       "e_initial_mwh", "eta_charge", "eta_discharge", "cost_per_mwh")]
            if not np.isfinite(vals).all() or min(vals) < 0:
                raise ModelDataError("Invalid storage parameters")
            if not (b["e_min_mwh"] <= b["e_initial_mwh"] <= b["e_max_mwh"]):
                raise ModelDataError("Initial storage energy outside bounds")
            if not (0 < b["eta_charge"] <= 1 and 0 < b["eta_discharge"] <= 1):
                raise ModelDataError("Efficiencies must lie in (0,1]")
            if b["eta_charge"] != 1 or b["eta_discharge"] != 1:
                raise ModelDataError("Two-stage model requires ideal efficiency eta_charge=eta_discharge=1; legacy lossy data must be revised explicitly")
        self.v0 = float(data["root_voltage_pu"])
        self.vmin = np.broadcast_to(_array(data["voltage_min_pu"]), (self.n,)).copy()
        self.vmax = np.broadcast_to(_array(data["voltage_max_pu"]), (self.n,)).copy()
        if not np.isfinite(self.v0) or (self.vmin <= 0).any() or (self.vmin > self.vmax).any():
            raise ModelDataError("Invalid voltage bounds")
        if not self.vmin[self.root] <= self.v0 <= self.vmax[self.root]:
            raise ModelDataError("Root voltage outside bounds")
        self.gmin = np.broadcast_to(_array(data["grid_min_mw"]), (self.t,)).copy()
        self.gmax = np.broadcast_to(_array(data["grid_max_mw"]), (self.t,)).copy()
        if (self.gmin > self.gmax).any():
            raise ModelDataError("grid_min_mw > grid_max_mw")
        # Legacy price metadata is not used: the new model has NO supplemental grid purchase.
        self.price = np.zeros(self.t)
        self.prices = _array(data["incentive_prices_per_mwh"], (4,), "incentive_prices_per_mwh")
        self.thresholds = _array(data["response_thresholds"], (self.k, 4), "response_thresholds")
        self.beta = float(data["beta"])
        self.rho = float(data["diagnostic_penalty"])
        if (self.price < 0).any() or (self.prices < 0).any() or (np.diff(self.prices) <= 0).any():
            raise ModelDataError("Prices must be nonnegative and the four incentive prices strictly increasing")
        if (self.thresholds <= 0).any() or (self.thresholds > 1).any() or (np.diff(self.thresholds, axis=1) < 0).any():
            raise ModelDataError("Response thresholds must be in (0,1] and nondecreasing")
        if not np.isfinite([self.beta, self.rho]).all() or min(self.beta, self.rho) <= 0:
            raise ModelDataError("beta and diagnostic_penalty must be positive")
        self.xlo = _array(data["xi_lower"], (self.k, self.t), "xi_lower")
        self.xhi = _array(data["xi_upper"], (self.k, self.t), "xi_upper")
        if (self.xlo > self.xhi).any():
            raise ModelDataError("xi_lower > xi_upper")
        self._network()
        for name,count in (("user_ids",self.k),("storage_ids",self.s)):
            if name in data and (len(data[name])!=count or len(set(data[name]))!=count):
                raise ModelDataError(f"{name} must contain {count} unique identifiers")

    def _network(self):
        branches = self.data["branches"]
        self.m = len(branches)
        if self.m != self.n - 1:
            raise ModelDataError("Require a connected radial network (N-1 in-service branches)")
        adjacency=[[] for _ in self.ids]
        for ell,b in enumerate(branches):
            if b["from_bus"] not in self.index or b["to_bus"] not in self.index:
                raise ModelDataError("Branch references an unknown bus")
            i,j=self.index[b["from_bus"]],self.index[b["to_bus"]]
            if i==j:
                raise ModelDataError("Branches cannot be self-loops; must be oriented on a tree")
            adjacency[i].append((j,ell));adjacency[j].append((i,ell))
        parent={};visited={self.root};queue=[self.root]
        self.r=np.zeros(self.m);self.x=np.zeros(self.m)
        self.line_min=np.zeros(self.m);self.line_max=np.zeros(self.m)
        normalized=[None]*self.m
        for i in queue:
            for j,ell in adjacency[i]:
                if j in visited:
                    continue
                visited.add(j);queue.append(j);parent[j]=(i,ell)
                b=copy.deepcopy(branches[ell])
                symmetric=b.get("p_limit_mw")
                lower=b.get("p_min_mw",-symmetric if symmetric is not None else None)
                upper=b.get("p_max_mw",symmetric)
                if lower is None or upper is None:
                    raise ModelDataError("Declare p_limit_mw or both p_min_mw/p_max_mw for each branch")
                lower,upper=float(lower),float(upper)
                if lower>upper:
                    raise ModelDataError("Branch lower bound exceeds upper bound")
                reverse=b["from_bus"]!=self.ids[i]
                if reverse: lower,upper=-upper,-lower
                b.update(from_bus=self.ids[i],to_bus=self.ids[j],p_min_mw=lower,p_max_mw=upper,
                         input_direction_reversed=reverse)
                normalized[ell]=b
                self.r[ell],self.x[ell]=b["r_ohm"],b["x_ohm"]
                self.line_min[ell],self.line_max[ell]=lower,upper
        if len(visited)!=self.n:
            raise ModelDataError("Disconnected network or cycle detected")
        if not np.isfinite([self.r,self.x,self.line_min,self.line_max]).all() or (self.r<0).any() or (self.x<0).any():
            raise ModelDataError("Require nonnegative finite ohmic impedances and finite power bounds")
        self.limits=np.maximum(abs(self.line_min),abs(self.line_max))
        self.oriented_branches=normalized
        self.path=np.zeros((self.n,self.m))
        for j in range(self.n):
            node=j
            while node!=self.root:
                node,ell=parent[node]
                self.path[j,ell]=1
        self.downstream=self.path.T
        self.vp=2*self.path@np.diag(self.r)@self.downstream/self.kv**2
        self.vq=2*self.path@np.diag(self.x)@self.downstream/self.kv**2


def _model(name, opts):
    m = gp.Model(name)
    m.Params.OutputFlag = 0
    m.Params.LogToConsole = int(opts.output_flag)
    m.Params.OutputFlag = int(opts.output_flag or bool(opts.log_dir))
    m.Params.Threads = opts.threads
    m.Params.Seed = opts.seed
    m.Params.TimeLimit = opts.time_limit
    m.Params.FeasibilityTol = min(1e-8, opts.feasibility_tol / 100)
    m.Params.OptimalityTol = 1e-8
    m.Params.IntFeasTol = 1e-8
    m.Params.MIPGap = opts.cost_rel_tol / 10
    m.Params.MIPGapAbs = min(opts.feasibility_tol / 10, opts.cost_abs_tol / 10)
    m.Params.DualReductions = 0
    if opts.log_dir:
        folder = Path(opts.log_dir)
        folder.mkdir(parents=True, exist_ok=True)
        m.Params.LogFile = str(folder / f"{name}_{time.time_ns()}.log")
    return m


def _bound(m, fallback):
    try:
        return float(m.ObjBound)
    except (gp.GurobiError, AttributeError):
        return float(m.ObjVal) if m.Status == GRB.OPTIMAL else fallback


def network_metrics(case, net_p):
    """Deterministic LinDistFlow check. Does not claim AC feasibility."""
    p = _array(net_p, (case.n, case.t), "net_p")
    flow = case.downstream @ p
    v2 = case.v0**2 - case.vp @ p - case.vq @ case.q
    return {
        "flow_mw":flow, "voltage_squared_pu":v2,
        "flow_mvar":case.downstream @ case.q,
        "net_q_mvar":case.q.copy(),
        "line_upper_margin_mw":case.line_max[:,None]-flow,
        "line_lower_margin_mw":flow-case.line_min[:,None],
        "voltage_lower_margin_squared_pu":v2-case.vmin[:,None]**2,
        "voltage_upper_margin_squared_pu":case.vmax[:,None]**2-v2,
        "voltage_pu":np.sqrt(np.maximum(v2, 0)),
        "min_voltage_pu":float(np.sqrt(max(0, v2.min()))),
        "max_voltage_pu":float(np.sqrt(max(0, v2.max()))),
        "max_line_loading":float(np.max(np.abs(flow) / np.maximum(1e-12,np.where(flow>=0,case.line_max[:,None],-case.line_min[:,None])))),
        "voltage_violation_count":int(np.count_nonzero((v2 < case.vmin[:, None]**2-1e-8) | (v2 > case.vmax[:, None]**2+1e-8))),
        "line_violation_count":int(np.count_nonzero((flow>case.line_max[:,None]+1e-8)|(flow<case.line_min[:,None]-1e-8)))
    }


def compute_unified_directrix(case):
    """Equations (1)--(6), independently usable before the nodal QP."""
    t,k,dt=case.t,case.k,case.dt
    eg = dt * (case.p.sum() - case.renew.sum()) + case.energy.sum()
    denominator = eg + case.energy.sum()
    if abs(denominator) < 1e-10:
        raise ModelDataError("Unified directrix denominator Eg + sum(Ek) is zero")
    u = (2 * eg / t + dt * (case.renew.sum(axis=0) - case.p.sum(axis=0))) / denominator
    ug = 2 / t - u
    g0 = eg / dt * ug
    return {"U":u,"UG":ug,"Eg_mwh":float(eg),"grid_plan_mw":g0,
            "user_energy_mwh":case.energy.copy(),"total_dr_energy_mwh":float(case.energy.sum()),
            "unified_denominator_mwh":float(denominator),"rigid_total_mw":case.p.sum(axis=0),
            "renewable_total_mw":case.renew.sum(axis=0)}


def compute_directrix(case, opts=None):
    """Equations (7)--(16): convex nodal QP, with separately exposed unified intermediates."""
    opts=opts or SolverOptions()
    t,k,dt=case.t,case.k,case.dt
    unified=compute_unified_directrix(case)
    u,ug,g0=unified["U"],unified["UG"],unified["grid_plan_mw"]
    eg,denominator=unified["Eg_mwh"],unified["unified_denominator_mwh"]
    if np.any(g0 < case.gmin-1e-8) or np.any(g0 > case.gmax+1e-8):
        raise ModelDataError("Unified exchange plan violates interconnection limits; change declared case assumptions")
    if not np.isclose(u.sum(),1,atol=1e-8) or not np.isclose(ug.sum(),1,atol=1e-8):
        raise ModelDataError("Unified directrix normalization failed")
    m = _model("nodal_directrix", opts)
    l = m.addMVar((k, t), lb=0, ub=1, name="L")
    alpha = l - u
    m.addConstr(l.sum(axis=1) == 1)
    m.addConstr(case.energy @ l == case.energy.sum()*u)
    wpower = case.w * (case.energy / dt)
    fixed = case.p - case.renew
    flow = case.downstream @ fixed + (case.downstream @ wpower) @ l
    v2 = case.v0**2 - case.vp @ fixed - (case.vp @ wpower) @ l - case.vq @ case.q
    m.addConstr(flow <= case.line_max[:, None])
    m.addConstr(flow >= case.line_min[:, None])
    m.addConstr(v2 >= case.vmin[:, None]**2)
    m.addConstr(v2 <= case.vmax[:, None]**2)
    m.setObjective((alpha * alpha).sum())
    m.optimize()
    if m.Status != GRB.OPTIMAL:
        status = m.Status
        m.dispose()
        raise ModelDataError(f"Nodal directrix was not solved to optimality (Gurobi status {status})")
    result = {"L":l.X, "U":u, "UG":ug, "Eg_mwh":float(eg), "grid_plan_mw":g0,
              "alpha":l.X-u, "target_dr_mw":case.energy[:,None]*l.X/dt,
              "user_energy_mwh":case.energy.copy(), "total_dr_energy_mwh":float(case.energy.sum()),
              "unified_denominator_mwh":float(denominator),
              "rigid_total_mw":case.p.sum(axis=0), "renewable_total_mw":case.renew.sum(axis=0),
              "alpha_squared":float(np.square(l.X-u).sum()),
              "solver_status":int(m.Status), "runtime_seconds":float(m.Runtime)}
    net = fixed + wpower @ l.X
    result["network"] = network_metrics(case, net)
    result["network"]["net_p_mw"] = net
    result["energy_residual_mwh"] = dt*result["target_dr_mw"].sum(axis=1)-case.energy
    result["aggregate_alpha_residual_mwh"] = case.energy @ result["alpha"]
    result["power_balance_residual_mw"] = net.sum(axis=0)-g0
    m.dispose()
    return result


class UncertaintySet:
    def __init__(self, case, directrix, levels):
        self.case = case
        self.l = _array(directrix["L"], (case.k, case.t), "L")
        levels = np.asarray(levels)
        if levels.shape != (case.k,) or not np.equal(levels, np.floor(levels)).all() or ((levels < 1) | (levels > 4)).any():
            raise ModelDataError("levels must contain one integer in 1..4 per user")
        self.levels = levels.astype(int)
        self.radius2 = -np.log(case.thresholds[np.arange(case.k), self.levels-1]) / case.beta
        self.lo, self.hi = np.maximum(case.xlo, -self.l), case.xhi.copy()
        if (self.lo > self.hi).any():
            raise ModelDataError("Empty uncertainty set: box contradicts nonnegative load")
        for k in np.flatnonzero(self.radius2==0):
            if (self.lo[k]>0).any() or (self.hi[k]<0).any():
                raise ModelDataError("Empty uncertainty set: minimum norm is positive but the response budget is zero")
            self.lo[k]=0;self.hi[k]=0
        centers = []
        for k in range(case.k):
            lo, hi = self.lo[k], self.hi[k]
            if lo.sum() > 1e-12 or hi.sum() < -1e-12:
                raise ModelDataError("Empty uncertainty set: box contradicts energy conservation")
            # Projection of zero onto {lo <= x <= hi, sum x=0}; unique minimum-norm point.
            left, right = float(lo.min()-1), float(hi.max()+1)
            for _ in range(100):
                mid = (left+right)/2
                if np.clip(mid, lo, hi).sum() > 0:
                    right = mid
                else:
                    left = mid
            center = np.clip((left+right)/2, lo, hi)
            if center @ center > self.radius2[k] + 1e-12:
                raise ModelDataError("Empty uncertainty set: minimum norm exceeds response budget")
            centers.append(center)
        self.center = np.array(centers)

    def contains(self, xi, tol=1e-8):
        x = np.asarray(xi, dtype=float)
        return bool(x.shape == self.lo.shape and np.isfinite(x).all()
                    and np.all(x >= self.lo-tol) and np.all(x <= self.hi+tol)
                    and np.max(np.abs(x.sum(axis=1))) <= tol
                    and np.all((x*x).sum(axis=1) <= self.radius2+tol))

    def parameters(self):
        c=self.case
        return {"levels":self.levels.copy(),"response_thresholds":c.thresholds[np.arange(c.k),self.levels-1],
                "incentive_prices_per_mwh":c.prices[self.levels-1],"beta":c.beta,
                "squared_deviation_budget":self.radius2.copy(),"declared_lower":c.xlo.copy(),
                "declared_upper":c.xhi.copy(),"effective_lower":self.lo.copy(),"effective_upper":self.hi.copy(),
                "minimum_norm_feasible_point":self.center.copy(),"minimum_squared_norm":np.sum(self.center**2,axis=1)}

    def add_to_model(self, model):
        x = model.addMVar(self.lo.shape, lb=self.lo, ub=self.hi, name="xi")
        model.addConstr(x.sum(axis=1) == 0)
        for k in range(self.case.k):
            model.addConstr(x[k] @ x[k] <= self.radius2[k])
        return x.reshape(-1)

    def project(self, value):
        """Return an admissible witness after solver rounding, without changing the set."""
        result = np.asarray(value, float).copy()
        for k in range(self.case.k):
            x = result[k]
            left, right = float((x-self.hi[k]).min()-1), float((x-self.lo[k]).max()+1)
            for _ in range(100):
                mid = (left+right)/2
                if np.clip(x-mid, self.lo[k], self.hi[k]).sum() > 0:
                    left = mid
                else:
                    right = mid
            x = np.clip(x-(left+right)/2, self.lo[k], self.hi[k])
            if x @ x > self.radius2[k]:
                center = self.center[k]
                direction = x-center
                aa, bb = direction @ direction, 2*(center @ direction)
                cc = center @ center-self.radius2[k]
                factor = (-bb+np.sqrt(max(0, bb*bb-4*aa*cc)))/(2*aa) if aa > 0 else 0
                x = center + max(0, min(1, factor)) * (1-1e-10) * direction
            result[k] = x
        return result

    def _ray(self, direction, fraction=1.0):
        result = self.center.copy()
        for k in range(self.case.k):
            d = direction[k] - direction[k].mean()
            norm = np.linalg.norm(d)
            if norm < 1e-14:
                continue
            d /= norm
            c = self.center[k]
            radius = -(c @ d) + np.sqrt(max(0, (c @ d)**2 + self.radius2[k] - c @ c))
            pos, neg = d > 1e-14, d < -1e-14
            if pos.any():
                radius = min(radius, np.min((self.hi[k, pos]-c[pos]) / d[pos]))
            if neg.any():
                radius = min(radius, np.min((self.lo[k, neg]-c[neg]) / d[neg]))
            result[k] += max(0, radius) * fraction * d
        return result

    def samples(self, count=100, seed=0):
        """Reproducible interior/boundary stress samples; NOT a uniform distribution or a proof."""
        rng = np.random.default_rng(seed)
        return [self._ray(rng.normal(size=self.lo.shape), 1.0 if i % 2 else rng.uniform(.2, 1)) for i in range(count)]

    def initial_scenarios(self):
        direction = np.tile(np.cos(2*np.pi*np.arange(self.case.t)/self.case.t), (self.case.k, 1))
        scenarios = [self.center, self._ray(direction), self._ray(-direction)]
        return _deduplicate(scenarios)


def _deduplicate(scenarios):
    out = []
    for x in scenarios:
        if not any(np.max(np.abs(x-a)) < 1e-9 for a in out):
            out.append(np.array(x, copy=True))
    return out


class RecourseSystem:
    """Ideal-storage exact continuous equivalent: Ay <= b+Bxi, Cy=e+Fxi.

    Z has zero columns solely for the internal legacy matrix-call convention.
    There is no common mode and no supplemental grid variable.
    """
    def __init__(self, case, directrix, diagnostic=False):
        self.case, self.directrix, self.diagnostic = case, directrix, diagnostic
        t, k, s = case.t, case.k, case.s
        st, kt = s*t, k*t
        self.nz, self.nx, self.ny = 0, kt, 2*st+kt
        self.ch = slice(0, st)
        self.dis = slice(st, 2*st)
        self.d = slice(2*st, self.ny)
        identity = sp.eye(self.ny, format="csr")
        ch, dis, d = (identity[sl] for sl in (self.ch, self.dis, self.d))
        ax, bx, zx, rx, scales = [], [], [], [], []
        ce, fe, re, eqscales = [], [], [], []

        def inequality(a, rhs, b=None, z=None):
            a = sp.csr_matrix(a)
            n = a.shape[0]
            rhs = np.broadcast_to(np.asarray(rhs, float), (n,)).copy()
            b = sp.csr_matrix((n, kt)) if b is None else sp.csr_matrix(b)
            z = sp.csr_matrix((n, 0)) if z is None else sp.csr_matrix(z)
            scale = np.maximum(1, np.maximum(np.abs(rhs), np.asarray(abs(a).max(axis=1).toarray()).ravel()))
            scale = np.maximum(scale, np.asarray(abs(b).max(axis=1).toarray()).ravel())
            inv = sp.diags(1/scale)
            ax.append(inv @ a); bx.append(inv @ b); zx.append(inv @ z); rx.append(rhs/scale); scales.append(scale)

        def equality(a, rhs, f=None):
            a = sp.csr_matrix(a)
            n = a.shape[0]
            rhs = np.broadcast_to(np.asarray(rhs, float), (n,)).copy()
            f = sp.csr_matrix((n, kt)) if f is None else sp.csr_matrix(f)
            scale = np.maximum(1, np.maximum(np.abs(rhs), np.asarray(abs(a).max(axis=1).toarray()).ravel()))
            scale = np.maximum(scale, np.asarray(abs(f).max(axis=1).toarray()).ravel())
            inv = sp.diags(1/scale)
            ce.append(inv @ a); fe.append(inv @ f); re.append(rhs/scale); eqscales.append(scale)

        self.c = np.zeros(self.ny)
        for j, b in enumerate(case.stores):
            self.c[j*t:(j+1)*t] = self.c[st+j*t:st+(j+1)*t] = case.dt*b["cost_per_mwh"]
        self.quad = np.zeros(self.ny)
        if diagnostic:
            self.quad[self.d] = case.rho
        if s:
            pc = np.repeat([float(b["p_charge_mw"]) for b in case.stores], t)
            pd = np.repeat([float(b["p_discharge_mw"]) for b in case.stores], t)
            inequality(-ch, 0); inequality(ch, pc)
            inequality(-dis, 0); inequality(dis, pd)
            etac = sp.eye(st,format="csr")
            etad = sp.eye(st,format="csr")
            accumulator = sp.kron(sp.eye(s), np.tril(np.ones((t, t))), format="csr")
            es = case.dt * accumulator @ (etac @ ch - etad @ dis)
            lower = np.repeat([b["e_min_mwh"]-b["e_initial_mwh"] for b in case.stores], t)
            upper = np.repeat([b["e_max_mwh"]-b["e_initial_mwh"] for b in case.stores], t)
            inequality(es, upper); inequality(-es, -lower)
            equality(es[np.arange(s)*t+t-1], 0)
        g0 = directrix["grid_plan_mw"]
        inequality(-d, directrix["L"].ravel(), b=sp.eye(kt))
        equality(sp.kron(sp.eye(k), np.ones((1, t))) @ d, 0)
        if not diagnostic:
            equality(d, 0)
        wd = sp.kron(case.w * (case.energy/case.dt), sp.eye(t), format="csr")
        ws = sp.kron(case.storage_map, sp.eye(t), format="csr")
        self.py = ws @ (ch-dis) + wd @ d
        self.px = wd
        self.p0 = case.p - case.renew + (case.w*(case.energy/case.dt)) @ directrix["L"]
        total = sp.kron(np.ones((1, case.n)), sp.eye(t), format="csr")
        # Fixed external exchange, equation (29). No additional power-purchase resource.
        equality(total @ self.py, g0-self.p0.sum(axis=0), f=-total @ wd)
        flow = sp.kron(case.downstream, sp.eye(t), format="csr")
        f0 = flow @ self.p0.ravel()
        upper_limit = np.repeat(case.line_max, t)
        lower_limit = np.repeat(case.line_min, t)
        inequality(flow @ self.py, upper_limit-f0, b=-flow @ wd)
        inequality(-flow @ self.py, -lower_limit+f0, b=flow @ wd)
        # Root voltage is prescribed and was checked when loading Case. Skip its zero rows here.
        indices = [i for i in range(case.n) if i != case.root]
        vp = sp.kron(case.vp[indices], sp.eye(t), format="csr")
        vq = sp.kron(case.vq[indices], sp.eye(t), format="csr")
        vbase = case.v0**2 - vp @ self.p0.ravel() - vq @ case.q.ravel()
        inequality(vp @ self.py, vbase-np.repeat(case.vmin[indices]**2, t), b=-vp @ wd)
        inequality(-vp @ self.py, np.repeat(case.vmax[indices]**2, t)-vbase, b=vp @ wd)
        self.A, self.B, self.Z = (sp.vstack(parts, format="csr") for parts in (ax, bx, zx))
        self.b, self.row_scales = np.concatenate(rx), np.concatenate(scales)
        self.C, self.F = (sp.vstack(parts, format="csr") for parts in (ce, fe))
        self.e, self.eq_scales = np.concatenate(re), np.concatenate(eqscales)

    def add(self, model, xi, z):
        y = model.addMVar(self.ny, lb=-GRB.INFINITY, name="recourse")
        model.addConstr(self.A @ y <= self.b + self.B @ np.asarray(xi).ravel() + self.Z @ z)
        model.addConstr(self.C @ y == self.e + self.F @ np.asarray(xi).ravel())
        return y

    def objective(self, y):
        value = self.c @ y
        if self.diagnostic:
            value = value + self.case.rho * (y[self.d] @ y[self.d])
        return value

    def violation(self, y, xi, z):
        return float(max(np.max(self.A @ y-self.b-self.B @ xi.ravel()-self.Z @ z, initial=0),
                         np.max(np.abs(self.C @ y-self.e-self.F @ xi.ravel()), initial=0)))

    def unpack(self, y, xi):
        c = self.case
        y = self.canonicalize(y)
        ch, dis = y[self.ch].reshape(c.s, c.t), y[self.dis].reshape(c.s, c.t)
        energy = np.zeros_like(ch)
        for s, b in enumerate(c.stores):
            energy[s] = b["e_initial_mwh"] + np.cumsum(c.dt*(b["eta_charge"]*ch[s]-dis[s]/b["eta_discharge"]))
        net_p = (self.p0.ravel()+self.py @ y+self.px @ xi.ravel()).reshape(c.n, c.t)
        capacity=np.array([b["e_max_mwh"] for b in c.stores])[:,None]
        initial=np.array([b["e_initial_mwh"] for b in c.stores])[:,None]
        actual=c.energy[:,None]/c.dt*(self.directrix["L"]+xi)
        auxiliary=y[self.d].reshape(c.k,c.t)
        return {"charge_mw":ch, "discharge_mw":dis, "energy_mwh":energy,
                "mode_z":(ch>0).astype(int), "storage_net_injection_mw":dis-ch,
                "energy_with_initial_mwh":np.concatenate([initial,energy],axis=1),
                "soc_fraction":np.divide(energy,capacity,out=np.zeros_like(energy),where=capacity>0),
                "grid_mw":self.directrix["grid_plan_mw"].copy(),
                "xi":xi.copy(), "response_dr_mw":actual,
                "response_degree":np.exp(-c.beta*np.sum(xi**2,axis=1)),
                "squared_response_deviation":np.sum(xi**2,axis=1),
                "throughput_mwh_per_storage":c.dt*(ch+dis).sum(axis=1),
                "storage_cost_per_device":c.dt*(ch+dis).sum(axis=1)*np.array([b["cost_per_mwh"] for b in c.stores]),
                "adjusted_dr_mw":actual+c.energy[:,None]/c.dt*auxiliary,
                "auxiliary_power_mw":c.energy[:,None]/c.dt*auxiliary,
                "cost_components":{"storage_throughput":float(self.c@y),
                                   "diagnostic_penalty":float((self.quad*y*y).sum())},
                "power_balance_residual_mw":net_p.sum(axis=0)-self.directrix["grid_plan_mw"],
                "terminal_energy_residual_mwh":energy[:,-1]-initial.ravel(),
                "d":y[self.d].reshape(c.k, c.t), "net_p_mw":net_p,
                "cost":float(self.c@y + (self.quad*y*y).sum()), **network_metrics(c, net_p)}

    def canonicalize(self, y):
        """Recover exact scenario-adaptive binary feasibility without changing any net injection."""
        y=np.asarray(y).copy()
        net=y[self.dis]-y[self.ch]
        y[self.ch],y[self.dis]=np.maximum(-net,0),np.maximum(net,0)
        return y


def solve_recourse(system, xi, z=None, opts=None, phase_one=False, explicit_mip=False):
    opts = opts or SolverOptions()
    if z is not None and np.size(z):
        raise ModelDataError("v2 optimizes modes separately for every response; fixed/common z is not accepted")
    z=np.empty(0)
    xi=_array(xi,(system.case.k,system.case.t),"xi")
    m = _model("phase_one" if phase_one else "feedback_dispatch", opts)
    if phase_one:
        y = m.addMVar(system.ny, lb=-GRB.INFINITY)
        slack = m.addMVar(len(system.b), lb=0)
        residual = m.addMVar(len(system.e), lb=0)
        m.addConstr(system.A @ y <= system.b+system.B @ xi.ravel()+system.Z @ z+slack)
        eq = system.C @ y-system.e-system.F @ xi.ravel()
        m.addConstr(eq <= residual); m.addConstr(-eq <= residual)
        m.setObjective(slack.sum()+residual.sum())
    else:
        y = system.add(m, xi, z)
        if explicit_mip and system.case.s:
            mode=m.addMVar((system.case.s,system.case.t),vtype=GRB.BINARY,name="scenario_mode")
            pc=np.array([b["p_charge_mw"] for b in system.case.stores])[:,None]
            pd=np.array([b["p_discharge_mw"] for b in system.case.stores])[:,None]
            m.addConstr(y[system.ch].reshape(mode.shape)<=pc*mode)
            m.addConstr(y[system.dis].reshape(mode.shape)<=pd*(1-mode))
        m.setObjective(system.objective(y))
    m.optimize()
    result = {"solver_status":int(m.Status), "status":"infeasible" if m.Status == GRB.INFEASIBLE else "unverified"}
    if m.Status == GRB.OPTIMAL:
        yy=y.X if phase_one else system.canonicalize(y.X)
        result.update(status="optimal", objective=float(m.ObjVal), y=yy,
                      runtime_seconds=float(m.Runtime),lower_bound=_bound(m,-float("inf")),
                      max_normalized_violation=system.violation(yy, xi, z) if not phase_one else None)
        if not phase_one:
            result["dispatch"] = system.unpack(yy, xi)
            result["objective"]=result["dispatch"]["cost"]
            if result["max_normalized_violation"]>opts.feasibility_tol:
                result["status"]="unverified"
        else:
            result.update(inequality_slack=slack.X,equality_slack=residual.X,
                          scope="Phase-I diagnostic variables; not an executable dispatch")
    m.dispose()
    return result


def solve_master(system, scenarios, opts):
    """With adaptive modes there are no shared decisions: finite max Q is separable."""
    rows=[solve_recourse(system,x,opts=opts) for x in scenarios]
    result={"solver_status":GRB.OPTIMAL,"status":"unverified",
            "lower_bound":max((r.get("lower_bound",-float("inf")) for r in rows),default=0),
            "scenario_results":[{"xi":x.copy(),**r} for x,r in zip(scenarios,rows)]}
    if any(r["status"]=="infeasible" for r in rows):
        result.update(status="proven_infeasible",solver_status=GRB.INFEASIBLE)
    elif all(r["status"]=="optimal" for r in rows):
        result.update(status="candidate",ys=[r["y"] for r in rows],
                      objective=max(r["objective"] for r in rows),
                      scenario_dispatches=[r["dispatch"] for r in rows])
    else:
        result["solver_status"]=next(r["solver_status"] for r in rows if r["status"]!="optimal")
    return result


def global_oracle(system, uncertainty, z, opts, feasibility=True):
    """Exact dual reformulation; NonConvex=2. No artificial dual big-M bounds.

    Phase I: max b(x,z)'lambda + e(x)'nu, A'lambda+C'nu=0,
             -1<=lambda<=0, -1<=nu<=1.
    Cost: maximize primal cost over primal feasibility and exact convex-QP KKT conditions.
          c+A'lambda+C'nu+2*diag(quad)*y=0, lambda>=0, nu free,
          slack=b(x,z)-Ay>=0, lambda'slack=0. No finite multiplier cutoffs.
    A finite valid upper bound is required to certify; an incumbent is only a witness.
    """
    if system.case.t==2 and system.case.k<=opts.two_period_vertex_limit:
        # Exact reduction, NOT endpoint sampling of a general ellipsoid:
        # xi=(a,-a), so 2*a^2<=B makes each user's feasible set an interval.
        radius=np.sqrt(uncertainty.radius2/2)
        lower=np.maximum.reduce([uncertainty.lo[:,0],-uncertainty.hi[:,1],-radius])
        upper=np.minimum.reduce([uncertainty.hi[:,0],-uncertainty.lo[:,1],radius])
        start=time.perf_counter();values=[];witnesses=[];complete=True
        for bits in itertools.product((0,1),repeat=system.case.k):
            aa=np.where(bits,upper,lower);xx=np.column_stack([aa,-aa])
            r=solve_recourse(system,xx,opts=opts,phase_one=feasibility)
            if r["status"]!="optimal":
                complete=False
            else:
                values.append(r["objective"]);witnesses.append(xx)
        best=int(np.argmax(values)) if values else None
        return {"solver_status":GRB.OPTIMAL if complete else GRB.TIME_LIMIT,
                "upper_bound":max(values) if complete and values else float("inf"),
                "witness_value":values[best] if best is not None else None,
                "xi":witnesses[best] if best is not None else None,
                "runtime_seconds":time.perf_counter()-start,
                "certificate":"exact two-period interval vertices; convex equivalent value function",
                "evaluated_vertex_count":2**system.case.k}
    if feasibility and system.diagnostic:
        # Model-specific exact certificate: d(xi)=-xi, ch=dis=0 restores the nodal baseline.
        # Verify the policy against the SAME matrices; do not assume it is valid for arbitrary inputs.
        ac = (-system.A[:,system.d]-system.B).toarray()
        ec = (-system.C[:,system.d]-system.F).toarray()
        for coeff in (ac,ec):
            view=coeff.reshape(len(coeff),system.case.k,system.case.t)
            view-=view.mean(axis=2,keepdims=True)  # sum_t xi_kt=0 exactly
        lo,hi=uncertainty.lo.ravel(),uncertainty.hi.ravel()
        worst_ineq=-system.b-system.Z @ z+np.maximum(ac*lo,ac*hi).sum(axis=1)
        upper_eq=-system.e+np.maximum(ec*lo,ec*hi).sum(axis=1)
        lower_eq=-system.e+np.minimum(ec*lo,ec*hi).sum(axis=1)
        certificate=float(np.maximum(worst_ineq,0).sum()+np.maximum(abs(upper_eq),abs(lower_eq)).sum())
        if certificate <= opts.feasibility_tol:
            return {"solver_status":GRB.OPTIMAL,"upper_bound":certificate,"witness_value":0.0,
                    "xi":uncertainty.center.copy(),"runtime_seconds":0.0,
                    "certificate":"verified affine policy d(xi)=-xi, ch=dis=0 over the entire uncertainty set"}
    m = _model("global_feasibility" if feasibility else "global_cost", opts)
    m.Params.NonConvex = 2
    x = uncertainty.add_to_model(m)
    x.Start = uncertainty.center.ravel()
    if feasibility:
        lam = m.addMVar(len(system.b), lb=-1, ub=0, name="inequality_dual")
        nu = m.addMVar(len(system.e), lb=-1, ub=1, name="equality_dual")
        m.addConstr(system.A.T @ lam + system.C.T @ nu == 0)
        obj = (system.b+system.Z @ z) @ lam + system.e @ nu
        obj += (system.B @ x) @ lam + (system.F @ x) @ nu
    else:
        lam = m.addMVar(len(system.b), lb=0, name="inequality_dual")
        nu = m.addMVar(len(system.e), lb=-GRB.INFINITY, name="equality_dual")
        # Explicit valid physical bounds keep global cost bounds finite without bounding multipliers.
        # d is a difference of two nonnegative vectors each summing to one, hence -1 <= d <= 1.
        lower, upper = np.zeros(system.ny), np.zeros(system.ny)
        c = system.case
        upper[system.ch] = np.repeat([b["p_charge_mw"] for b in c.stores],c.t)
        upper[system.dis] = np.repeat([b["p_discharge_mw"] for b in c.stores],c.t)
        if system.diagnostic:
            lower[system.d], upper[system.d] = -1, 1
        y = m.addMVar(system.ny,lb=lower,ub=upper,name="optimal_recourse")
        slack = m.addMVar(len(system.b),lb=0,name="physical_slack")
        m.addConstr(system.A @ y + slack == system.b+system.Z @ z+system.B @ x)
        m.addConstr(system.C @ y == system.e+system.F @ x)
        stationarity = system.c + system.A.T @ lam + system.C.T @ nu + 2*(sp.diags(system.quad) @ y)
        m.addConstr(stationarity == 0)
        # Nonnegative summands: one equality is equivalent to every complementary pair being zero.
        m.addConstr(lam @ slack == 0)
        obj = system.objective(y)
    m.setObjective(obj, GRB.MAXIMIZE)
    m.optimize()
    # Numeric/suboptimal/infeasible statuses never establish a global certificate.
    usable_status = m.Status in (GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.NODE_LIMIT, GRB.ITERATION_LIMIT, GRB.WORK_LIMIT)
    result = {"solver_status":int(m.Status), "upper_bound":_bound(m, float("inf")) if usable_status else float("inf"),
              "witness_value":None, "xi":None, "runtime_seconds":float(m.Runtime)}
    if m.SolCount and usable_status:
        xx = x.X.reshape(uncertainty.lo.shape)
        xx = uncertainty.project(xx)
        if uncertainty.contains(xx, tol=1e-10):
            # Raw dual objective refers to the unrounded optimizer, not necessarily this repaired witness.
            result.update(witness_value=float(m.ObjVal), xi=xx)
    m.dispose()
    return result


def solve_robust(case, directrix, levels, opts=None, diagnostic=False, retained_scenarios=None):
    opts = opts or SolverOptions()
    uncertainty = UncertaintySet(case, directrix, levels)
    system = RecourseSystem(case, directrix, diagnostic)
    scenarios = uncertainty.initial_scenarios()
    if retained_scenarios:
        scenarios.extend(x for x in retained_scenarios if uncertainty.contains(x, tol=1e-9))
    scenarios = _deduplicate(scenarios)
    history = []
    result = {"status":"unverified", "robust_certified":False, "diagnostic":diagnostic,
              "levels":np.asarray(levels), "history":history, "scenarios":scenarios,
              "uncertainty":uncertainty.parameters(),"scenario_results":[],
              "mode_policy":"scenario_adaptive","feasibility_certified":False}
    for iteration in range(opts.max_iterations):
        master = solve_master(system, scenarios, opts)
        result["scenario_results"]=master["scenario_results"]
        event = {"iteration":iteration+1, "scenario_count":len(scenarios),
                 "master_status":master["solver_status"], "lower_bound":master["lower_bound"],
                 "scenario_results":master["scenario_results"]}
        history.append(event)
        if master["status"] == "proven_infeasible":
            return {**result, "status":"proven_infeasible", "reason":"Finite subset master is proven infeasible"}
        if master["status"] != "candidate":
            return {**result, "reason":"Master has no verified incumbent"}
        z = np.empty(0)
        result.update(lower_bound=master["lower_bound"], master_cost=master["objective"])
        feasibility = global_oracle(system, uncertainty, z, opts, True)
        event["feasibility_upper_bound"] = feasibility["upper_bound"]
        event["feasibility_solver_status"] = feasibility["solver_status"]
        event["feasibility_oracle"]=feasibility
        witness = feasibility["xi"]
        if witness is not None:
            check = solve_recourse(system, witness, z, opts, phase_one=True)
            event["feasibility_witness"] = check.get("objective")
            event["feasibility_witness_check"] = check
            if check["status"] == "optimal" and check["objective"] > opts.feasibility_tol:
                extended = _deduplicate(scenarios+[witness])
                if len(extended) == len(scenarios):
                    return {**result, "reason":"Repeated violating witness; numerical consistency needs review"}
                scenarios[:] = extended
                continue
        if feasibility["upper_bound"] > opts.feasibility_tol or not np.isfinite(feasibility["upper_bound"]):
            return {**result, "reason":"Global feasibility upper bound did not close; do not upgrade incentives"}
        result["feasibility_certified"]=True
        cost = global_oracle(system, uncertainty, z, opts, False)
        ub = cost["upper_bound"]
        event.update(cost_upper_bound=ub, cost_solver_status=cost["solver_status"])
        event["cost_oracle"]=cost
        if cost["xi"] is not None:
            result["worst_cost_witness"]={"xi":cost["xi"],**solve_recourse(system,cost["xi"],opts=opts)}
        gap = ub-master["lower_bound"]
        tol = opts.cost_abs_tol+opts.cost_rel_tol*max(1, abs(master["lower_bound"]))
        result.update(upper_bound=ub, gap=gap)
        if np.isfinite(ub) and -tol <= gap <= tol:
            # Re-optimize every scenario independently for the diagnostic indicator.
            recourses = [solve_recourse(system, x, z, opts) for x in scenarios]
            if not all(r["status"] == "optimal" and r["max_normalized_violation"] <= opts.feasibility_tol for r in recourses):
                return {**result, "reason":"Scenario-adaptive recourse verification failed"}
            indicator = np.max([np.max(np.abs(r["dispatch"]["d"]), axis=1) for r in recourses], axis=0)
            return {**result, "status":"robust_optimal", "robust_certified":True,
                    "diagnostic_indicator":indicator, "reason":"Global feasibility and cost bounds verified within tolerances"}
        if cost["xi"] is None:
            return {**result, "reason":"Global cost search has no admissible witness / no finite closing upper bound"}
        extended = _deduplicate(scenarios+[cost["xi"]])
        if len(extended) == len(scenarios):
            return {**result, "reason":"Cost witness repeated while global bounds remain open"}
        scenarios[:] = extended
    return {**result, "reason":"Scenario-generation iteration limit reached"}


def run_incentive_loop(case, opts=None, directrix=None):
    opts = opts or SolverOptions()
    directrix = directrix or compute_directrix(case, opts)
    levels = np.ones(case.k, dtype=int)
    rounds, retained = [], []
    for _ in range(3*case.k+1):
        try:
            operational = solve_robust(case, directrix, levels, opts, retained_scenarios=retained)
        except ModelDataError as exc:
            return {"directrix":directrix,"levels":levels.copy(),"rounds":rounds,
                    "status":"invalid_uncertainty","robust_certified":False,"reason":str(exc)}
        item = {"levels":levels.copy(), "operational":operational}
        rounds.append(item)
        base = {"directrix":directrix, "levels":levels.copy(), "rounds":rounds,
                "status":operational["status"], "robust_certified":False,
                "feasibility_certified":operational.get("feasibility_certified",False),
                "cost_optimality_certified":operational["status"]=="robust_optimal"}
        if operational["status"] == "robust_optimal":
            return {**base, "robust_certified":True}
        if operational["status"] != "proven_infeasible":
            return {**base, "reason":operational["reason"]}
        diagnostic = solve_robust(case, directrix, levels, opts, diagnostic=True,
                                  retained_scenarios=operational["scenarios"])
        item["diagnostic"] = diagnostic
        if diagnostic["status"] != "robust_optimal":
            return {**base, "status":"proven_infeasible" if diagnostic["status"] == "proven_infeasible" else "unverified",
                    "reason":"Diagnostic model: " + diagnostic["reason"]}
        upgrade = (diagnostic["diagnostic_indicator"] > opts.upgrade_tol) & (levels < 4)
        if not upgrade.any():
            needs_upgrade=diagnostic["diagnostic_indicator"] > opts.upgrade_tol
            return {**base, "status":"no_feasible_scheme" if needs_upgrade.any() else "stalled",
                    "reason":"Required users are at maximum tier" if needs_upgrade.any() else "No diagnostic indicator exceeds upgrade tolerance"}
        levels[upgrade] += 1
        item["upgrade_mask"]=upgrade.copy()
        item["next_levels"]=levels.copy()
        retained = diagnostic["scenarios"]  # next solve filters them against the smaller set
    raise RuntimeError("Internal error: incentive loop exceeded its finite tier bound")


def validate_feedback(case, directrix, levels, feedback_p_mw):
    p = _array(feedback_p_mw, (case.k, case.t), "feedback_p_mw")
    uncertainty = UncertaintySet(case, directrix, levels)
    xi = case.dt*p/case.energy[:, None]-directrix["L"]
    feedback_energy=case.dt*p.sum(axis=1)
    return {"xi":xi,"in_set":uncertainty.contains(xi,tol=1e-9),
            "energy_changed":not np.allclose(feedback_energy,case.energy,rtol=1e-9,atol=1e-9),
            "feedback_energy_mwh":feedback_energy,"energy_residual_mwh":feedback_energy-case.energy,
            "response_degree":np.exp(-case.beta*np.sum(xi**2,axis=1)),
            "budget_margin":uncertainty.radius2-np.sum(xi**2,axis=1),
            "lower_margin":xi-uncertainty.lo,"upper_margin":uncertainty.hi-xi,
            "nonnegative_load_margin_mw":p.copy(),"zero_sum_residual":xi.sum(axis=1)}


def dispatch_feedback(case, directrix, levels, feedback_p_mw, opts=None):
    """Optimize adaptive modes after checking the entire confirmed feedback cycle."""
    opts = opts or SolverOptions()
    check=validate_feedback(case,directrix,levels,feedback_p_mw)
    if check["energy_changed"]:
        return {"status":"energy_changed","requires_stage1_rebuild":True,"feedback_check":check,
                "planned_energy_mwh":case.energy.copy(),"reason":"Cycle energy changed; rebuild both stages per section 3.3"}
    if not check["in_set"]:
        return {"status":"out_of_set","feedback_check":check,"reason":"Feedback outside declared set; update parameters and repeat verification"}
    system = RecourseSystem(case, directrix, diagnostic=False)
    result = solve_recourse(system, check["xi"], opts=opts)
    result["feedback_check"]=check
    return result


def solve_scenario(case, directrix, xi, opts=None, diagnostic=False, explicit_mip=False):
    """Optimize one confirmed scenario with adaptive modes; explicit_mip independently checks equivalence."""
    system=RecourseSystem(case,directrix,diagnostic)
    return solve_recourse(system,xi,opts=opts,explicit_mip=explicit_mip)


def run_two_stage(case, opts=None):
    """Stable complete-result API. Returns intermediate stages even if global verification is unfinished."""
    opts=opts or SolverOptions()
    directrix={}
    try:
        directrix=compute_unified_directrix(case)
        directrix=compute_directrix(case,opts)
    except ModelDataError as exc:
        stage2={"status":"stage1_failed","robust_certified":False,"rounds":[],"reason":str(exc)}
    else:
        stage2=run_incentive_loop(case,opts,directrix)
    return {"result_schema_version":2,"model":"ideal_storage_adaptive_modes_fixed_exchange",
            "status":stage2["status"],"robust_certified":stage2["robust_certified"],
            "reason":stage2.get("reason","Global feasibility and cost bounds verified"),
            "solution_method":"exact_continuous_equivalent_with_scenario_mode_recovery",
            "feasibility_certified":stage2.get("feasibility_certified",False),
            "cost_optimality_certified":stage2.get("cost_optimality_certified",False),
            "input_case":copy.deepcopy(case.data),"solver_options":vars(opts).copy(),
            "axes":{"bus_ids":case.ids.copy(),"user_ids":case.data.get("user_ids",list(range(case.k))),
                    "storage_ids":case.data.get("storage_ids",list(range(case.s))),
                    "time_start_hours":np.arange(case.t)*case.dt,
                    "state_time_hours":np.arange(case.t+1)*case.dt,
                    "branches":case.oriented_branches,"array_order":"entity,time"},
            "parameters":{"user_energy_mwh":case.energy.copy(),"dr_allocation":case.w.copy(),
                          "downstream_matrix":case.downstream.copy(),"path_matrix":case.path.copy(),
                          "voltage_active_sensitivity":case.vp.copy(),"voltage_reactive_sensitivity":case.vq.copy(),
                          "fixed_reactive_flow_mvar":case.downstream@case.q,
                          "assumptions":{"ideal_efficiency":True,"fixed_exchange":True,"scenario_adaptive_modes":True,
                                         "lossless_lindistflow":True,"full_cycle_feedback_before_dispatch":True}},
            "stage1":directrix,"stage2":stage2}


def export_result(result, path, include_matrices=False):
    """Write portable strict JSON; optionally also expose both recourse systems as sparse NPZ files."""
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    artifact=copy.deepcopy(result)
    artifact["source_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if include_matrices and "L" in result["stage1"]:
        case=Case(result["input_case"])
        artifacts={}
        for diagnostic in (False,True):
            system=RecourseSystem(case,result["stage1"],diagnostic)
            label="diagnostic" if diagnostic else "operational"
            folder=path.parent/(path.stem+"_matrices")/label
            folder.mkdir(parents=True,exist_ok=True)
            for name in ("A","B","C","F"):
                sp.save_npz(folder/(name+".npz"),getattr(system,name))
            np.savez(folder/"vectors.npz",b=system.b,e=system.e,c=system.c,quad=system.quad,
                     row_scales=system.row_scales,eq_scales=system.eq_scales)
            artifacts[label]={"directory":str(folder.resolve()),
                              "y_slices":{name:[getattr(system,name).start,getattr(system,name).stop]
                                          for name in ("ch","dis","d")},
                              "xi_shape":[case.k,case.t],"constraints":"A@y<=b+B@xi; C@y=e+F@xi"}
        artifact["matrix_exports"]=artifacts
    path.write_text(json.dumps(json_ready(artifact),ensure_ascii=False,indent=2,allow_nan=False),encoding="utf-8")
    return path.resolve()


def json_ready(value):
    """Portable strict JSON: unavailable/infinite solver bounds are null, never NaN/Infinity."""
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, dict):
        return {key:json_ready(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value
