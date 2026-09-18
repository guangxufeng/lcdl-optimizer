import copy
import json

import numpy as np
import pytest
from scipy import sparse

from lcdl_algorithm import (Case, ModelDataError, SolverOptions, compute_directrix, compute_unified_directrix,
                            solve_scenario, dispatch_feedback, run_two_stage, export_result, network_metrics,
                            UncertaintySet, solve_master, RecourseSystem)
from test_algorithm import tiny_data


def storage(bus=99,cost=1):
    return {"bus":bus,"p_charge_mw":.3,"p_discharge_mw":.4,"e_min_mwh":.1,"e_max_mwh":2,
            "e_initial_mwh":1,"cost_per_mwh":cost}


@pytest.mark.parametrize("diagnostic",[False,True])
@pytest.mark.parametrize("cost",[0,1,13])
def test_continuous_equivalent_matches_binary_model(diagnostic,cost):
    data=tiny_data();data["storage"]=[storage(cost=cost)]
    c=Case(data);d=compute_directrix(c)
    for a in [-.05,-.012,.021,.05]:
        xi=np.array([[a,-a]])
        continuous=solve_scenario(c,d,xi,diagnostic=diagnostic)
        binary=solve_scenario(c,d,xi,diagnostic=diagnostic,explicit_mip=True)
        assert continuous["status"]==binary["status"]=="optimal"
        assert continuous["objective"]==pytest.approx(binary["objective"],abs=1e-6)
        r=continuous["dispatch"]
        assert np.max(r["charge_mw"]*r["discharge_mw"])==0
        assert np.allclose(r["grid_mw"],d["grid_plan_mw"])
        assert np.max(np.abs(r["terminal_energy_residual_mwh"]))<1e-7


def test_finite_scenarios_do_not_share_modes():
    data=tiny_data();data["storage"]=[storage()]
    c=Case(data);d=compute_directrix(c)
    scenarios=[np.array([[.05,-.05]]),np.array([[-.05,.05]])]
    result=solve_master(RecourseSystem(c,d),scenarios,SolverOptions())
    assert result["status"]=="candidate"
    modes=[r["dispatch"]["mode_z"].tolist() for r in result["scenario_results"]]
    assert modes==[[[0,1]],[[1,0]]]


def test_reversed_branch_and_asymmetric_limits():
    data=tiny_data()
    data["branches"]=[{"from_bus":99,"to_bus":10,"r_ohm":.1,"x_ohm":.05,"p_min_mw":-2,"p_max_mw":.3}]
    c=Case(data);d=compute_directrix(c)
    assert c.line_min.tolist()==[-.3] and c.line_max.tolist()==[2]
    assert c.oriented_branches[0]["from_bus"]==10
    assert c.oriented_branches[0]["input_direction_reversed"]
    assert d["network"]["flow_mw"][0,0]==pytest.approx(1)
    assert d["network"]["line_upper_margin_mw"][0,0]==pytest.approx(1)


@pytest.mark.parametrize("topology",["chain","star"])
def test_different_five_bus_topologies_and_heterogeneous_loads(topology):
    data=tiny_data(dt=.25)
    data["bus_ids"]=[90,5,77,31,12];data["root_bus"]=77
    data["rigid_p_mw"]=[[.1,.14],[.2,.16],[0,0],[.12,.08],[.05,.06]]
    data["fixed_q_mvar"]=[[.02,.03],[.04,.05],[0,0],[.01,.01],[.01,.02]]
    data["renewable_p_mw"]=[[0,0],[.1,.1],[0,0],[0,0],[0,0]]
    data["dr_pre_mw"]=[[.1,.15],[.2,.12]]
    data["dr_allocation"]=[[1,0],[0,0],[0,0],[0,1],[0,0]]
    data["xi_lower"]=[[-.005,-.005]]*2;data["xi_upper"]=[[.005,.005]]*2
    data["response_thresholds"]=[[.9,.95,.98,.99]]*2
    data["storage"]=[storage(90),storage(31)]
    edges=[(77,90),(90,5),(5,31),(31,12)] if topology=="chain" else [(77,90),(77,5),(77,31),(77,12)]
    data["branches"]=[{"from_bus":j,"to_bus":i,"r_ohm":.03,"x_ohm":.01,"p_limit_mw":3} for i,j in edges[::-1]]
    c=Case(data);d=compute_directrix(c)
    assert d["L"].shape==(2,2)
    assert np.allclose(d["power_balance_residual_mw"],0,atol=1e-8)
    r=solve_scenario(c,d,np.array([[.003,-.003],[-.002,.002]]))
    assert r["status"]=="optimal"
    assert r["dispatch"]["flow_mw"].shape==(4,2)
    assert r["dispatch"]["energy_with_initial_mwh"].shape==(2,3)


def test_feedback_energy_change_requires_both_stages():
    c=Case(tiny_data());d=compute_directrix(c)
    result=dispatch_feedback(c,d,[1],[[.4,.5]])
    assert result["status"]=="energy_changed"
    assert result["requires_stage1_rebuild"]


@pytest.mark.parametrize("with_storage",[False,True])
def test_complete_export_and_sparse_matrix_roundtrip(tmp_path,with_storage):
    data=tiny_data(deterministic=True)
    if with_storage: data["storage"]=[storage()]
    c=Case(data)
    r=run_two_stage(c,SolverOptions(time_limit=5))
    assert r["robust_certified"]
    out=export_result(r,tmp_path/"result.json",include_matrices=True)
    saved=json.loads(out.read_text(encoding="utf-8"))
    stage=saved["stage1"]
    for key in ["U","UG","alpha","L","target_dr_mw","user_energy_mwh","network"]:
        assert key in stage
    assert saved["stage2"]["rounds"][0]["operational"]["uncertainty"]["squared_deviation_budget"]==[0]
    row=r["stage2"]["rounds"][0]["operational"]["scenario_results"][0]
    folder=tmp_path/"result_matrices"/"operational"
    a,b=sparse.load_npz(folder/"A.npz"),sparse.load_npz(folder/"B.npz")
    v=np.load(folder/"vectors.npz")
    assert np.max(a@row["y"]-v["b"]-b@row["xi"].ravel())<1e-7
    assert row["dispatch"]["energy_with_initial_mwh"].shape==(int(with_storage),3)


def test_stage1_failure_keeps_unified_parameters(tmp_path):
    data=tiny_data();data["grid_max_mw"]=.5
    result=run_two_stage(Case(data))
    assert result["status"]=="stage1_failed"
    assert np.allclose(result["stage1"]["U"],[.5,.5])
    export_result(result,tmp_path/"failed.json",include_matrices=True)


def test_empty_uncertainty_keeps_stage1_outputs():
    data=tiny_data();data["xi_lower"]=[[.01,.01]]
    result=run_two_stage(Case(data))
    assert result["status"]=="invalid_uncertainty"
    assert "target_dr_mw" in result["stage1"]


def test_disconnected_cycle_is_rejected():
    data=tiny_data()
    data["bus_ids"]=[10,20,30,99]
    for key in ["rigid_p_mw","fixed_q_mvar","renewable_p_mw"]:
        data[key]=[[0,0],[0,0],[0,0],data[key][1]]
    data["dr_allocation"]=[[0],[0],[0],[1]]
    data["branches"]=[{"from_bus":i,"to_bus":j,"r_ohm":.1,"x_ohm":.1,"p_limit_mw":2}
                      for i,j in [(20,30),(30,99),(99,20)]]
    with pytest.raises(ModelDataError,match="Disconnected"):
        Case(data)


def test_perfect_response_rejects_small_nonzero_deviation():
    c=Case(tiny_data(deterministic=True));d=compute_directrix(c)
    u=UncertaintySet(c,d,[1])
    assert not u.contains(np.array([[1e-6,-1e-6]]))
    result=dispatch_feedback(c,d,[1],[[.400001,.399999]])
    assert result["status"]=="out_of_set"


def test_nontrivial_adaptive_robust_certificate_matches_analytic_cost():
    data=tiny_data();data["storage"]=[storage()]
    r=run_two_stage(Case(data),SolverOptions(time_limit=5))
    assert r["robust_certified"]
    op=r["stage2"]["rounds"][0]["operational"]
    assert op["upper_bound"]==pytest.approx(.08,abs=1e-7)
    assert op["history"][0]["cost_oracle"]["evaluated_vertex_count"]==2
