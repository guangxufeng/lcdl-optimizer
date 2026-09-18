import copy

import numpy as np
import pytest

from lcdl_algorithm import (Case, ModelDataError, SolverOptions, compute_directrix, UncertaintySet,
                            RecourseSystem, solve_recourse, solve_master, global_oracle,
                            solve_robust, run_incentive_loop, dispatch_feedback, network_metrics)


def tiny_data(dt=1.0, deterministic=False):
    return {
        "schema_version":1, "name":"two_bus_analytic", "bus_ids":[10, 99], "root_bus":10,
        "base_kv":12.66, "dt_hours":dt,
        "branches":[{"from_bus":10,"to_bus":99,"r_ohm":.1,"x_ohm":.05,"p_limit_mw":5}],
        "rigid_p_mw":[[0,0],[.6,.6]], "fixed_q_mvar":[[0,0],[.1,.1]],
        "renewable_p_mw":[[0,0],[0,0]], "dr_pre_mw":[[.4,.4]], "dr_allocation":[[0],[1]],
        "storage":[], "root_voltage_pu":1, "voltage_min_pu":.9, "voltage_max_pu":1.1,
        "grid_min_mw":0, "grid_max_mw":5,
        "incentive_prices_per_mwh":[1,2,3,4], "response_thresholds":[[1,1,1,1] if deterministic else [.9,.95,.98,.99]],
        "beta":1, "xi_lower":[[-.05,-.05]], "xi_upper":[[.05,.05]], "diagnostic_penalty":100
    }


def options():
    return SolverOptions(time_limit=15, max_iterations=5, threads=2,two_period_vertex_limit=0)




@pytest.mark.parametrize("dt", [.25,1,6])
def test_units_and_noncontiguous_bus_ids(dt):
    c = Case(tiny_data(dt))
    d = compute_directrix(c, options())
    assert np.allclose(d["L"], .5)
    assert np.allclose(d["grid_plan_mw"],1)
    assert d["Eg_mwh"] == pytest.approx(2*dt)
    m = network_metrics(c,np.array([[0,0],[1,1]]))
    assert m["voltage_squared_pu"][1,0] == pytest.approx(1-2*(.1+.05*.1)/12.66**2)




def test_invalid_tree_rejected():
    data = tiny_data()
    data["branches"][0]["to_bus"] = 10
    with pytest.raises(ModelDataError, match="oriented"):
        Case(data)


def test_empty_uncertainty_is_not_robust_feasibility():
    data = tiny_data()
    data["xi_lower"] = [[.01,.01]]
    c = Case(data)
    d = compute_directrix(c, options())
    with pytest.raises(ModelDataError, match="Empty uncertainty"):
        UncertaintySet(c,d,[1])
    data = tiny_data(deterministic=True)
    data["xi_lower"], data["xi_upper"] = [[.01,-.02]], [[.02,-.01]]
    c = Case(data)
    d = compute_directrix(c, options())
    with pytest.raises(ModelDataError, match="minimum norm"):
        UncertaintySet(c,d,[1])




def test_global_feasibility_oracle_matches_analytic_counterexample():
    c = Case(tiny_data()); d = compute_directrix(c, options())
    u = UncertaintySet(c,d,[1]); s = RecourseSystem(c,d)
    oracle = global_oracle(s,u,np.empty(0),options(),True)
    # Without storage and with a fixed grid, either nonzero deviation requires a violation.
    explicit = solve_recourse(s,np.array([[.05,-.05]]),np.empty(0),options(),True)
    assert explicit["objective"] == pytest.approx(.08,abs=1e-7)
    assert oracle["upper_bound"] == pytest.approx(.08,abs=1e-6)
    assert oracle["xi"] is not None
    witness = solve_recourse(s,oracle["xi"],np.empty(0),options(),True)
    assert witness["objective"] == pytest.approx(.08,abs=1e-6)


def test_diagnostic_cost_oracle_matches_analytic_maximum():
    c = Case(tiny_data()); d = compute_directrix(c, options())
    u = UncertaintySet(c,d,[1]); s = RecourseSystem(c,d,diagnostic=True)
    # No storage and fixed external exchange force d=-xi.
    # Max cost = rho * (0.05^2 + 0.05^2) = 0.5.
    oracle = global_oracle(s,u,np.empty(0),options(),False)
    assert oracle["upper_bound"] == pytest.approx(.5,abs=.005)
    assert oracle["xi"] is not None
    primal = solve_recourse(s,oracle["xi"],np.empty(0),options())
    assert primal["objective"] == pytest.approx(.5,abs=.005)


def test_no_storage_nonzero_uncertainty_is_proven_infeasible():
    c = Case(tiny_data()); d = compute_directrix(c, options())
    r = solve_robust(c,d,[1],options())
    assert r["status"] == "proven_infeasible"
    assert not r["robust_certified"]


def test_singleton_set_certifies_and_dispatches():
    c = Case(tiny_data(deterministic=True)); d = compute_directrix(c,options())
    result = run_incentive_loop(c,options(),d)
    assert result["robust_certified"]
    r = dispatch_feedback(c,d,[1],[[.4,.4]],options())
    assert r["status"] == "optimal"
    assert r["dispatch"]["cost"] == pytest.approx(0,abs=1e-7)
    r = dispatch_feedback(c,d,[1],[[.41,.39]],options())
    assert r["status"] == "out_of_set"


def test_ideal_storage_adaptive_modes_and_fixed_grid():
    data = tiny_data()
    data["storage"] = [{"bus":99,"p_charge_mw":1,"p_discharge_mw":1,"e_min_mwh":0,"e_max_mwh":2,
                        "e_initial_mwh":1,"eta_charge":1,"eta_discharge":1,"cost_per_mwh":1}]
    c=Case(data); d=compute_directrix(c,options()); s=RecourseSystem(c,d)
    xi=np.array([[-.05,.05]])
    r=solve_recourse(s,xi,opts=options())
    assert r["status"] == "optimal"
    schedule=r["dispatch"]
    assert schedule["energy_mwh"][0,-1] == pytest.approx(1,abs=1e-8)
    assert np.allclose(schedule["grid_mw"],d["grid_plan_mw"])
    assert schedule["mode_z"].tolist()==[[1,0]]
    reverse=solve_recourse(s,-xi,opts=options())
    assert reverse["dispatch"]["mode_z"].tolist()==[[0,1]]
    explicit=solve_recourse(s,xi,opts=options(),explicit_mip=True)
    assert explicit["objective"]==pytest.approx(r["objective"],abs=1e-7)
    assert "supplemental_grid_mw" not in schedule
    data["storage"][0]["eta_charge"]=.95
    with pytest.raises(ModelDataError,match="ideal efficiency"):
        Case(data)


def test_unverified_result_never_triggers_upgrade(monkeypatch):
    import lcdl_algorithm as algorithm
    calls=[]
    def timeout(*args,**kwargs):
        calls.append(kwargs)
        return {"status":"unverified","reason":"global search time limit","scenarios":[]}
    monkeypatch.setattr(algorithm,"solve_robust",timeout)
    c=Case(tiny_data()); d=compute_directrix(c,options())
    result=run_incentive_loop(c,options(),d)
    assert result["status"] == "unverified"
    assert len(calls) == 1
    assert result["levels"].tolist() == [1]


def test_proven_infeasibility_diagnosis_upgrade_and_execution():
    # Deliberate mathematical regression ONLY: final threshold=1, unlike the IEEE33 default=0.98.
    data=tiny_data()
    data["response_thresholds"]=[[.9,.95,.99,1]]
    c=Case(data); d=compute_directrix(c,options())
    result=run_incentive_loop(c,options(),d)
    assert result["status"] == "robust_optimal"
    assert result["levels"].tolist() == [4]
    assert len(result["rounds"]) == 4
    for item in result["rounds"][:-1]:
        assert item["operational"]["status"] == "proven_infeasible"
        assert item["diagnostic"]["status"] == "robust_optimal"
        assert item["diagnostic"]["diagnostic_indicator"][0] > 0


def test_quiet_console_still_writes_solver_audit_log(tmp_path):
    c=Case(tiny_data())
    compute_directrix(c,SolverOptions(log_dir=str(tmp_path),output_flag=False))
    logs=list(tmp_path.glob("*.log"))
    assert logs and "Optimal objective" in logs[0].read_text()
