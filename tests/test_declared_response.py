import copy
import json
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from lcdl_algorithm import (Case, SolverOptions, ModelDataError, compute_directrix,
    make_declaration, solve_declared_dispatch, run_two_stage, dispatch_feedback,
    export_result, with_declaration_model, UncertaintySet, DeclaredUncertaintySet,
    declared_global_oracle, DeclarationSystem)


def data():
    return json.loads((Path(__file__).resolve().parents[1]/'data/minimal_case.json').read_text(encoding='utf-8'))


def test_three_round_analytic_robust_bounds_and_feedback():
    r=run_two_stage(Case(data()),feedback_p_mw=[[.47,.33]])
    assert r['status']=='ready_to_execute' and r['executable'] and r['robust_certified']
    rounds=r['stage2']['rounds']
    assert [x['levels'].tolist() for x in rounds]==[[1],[2],[3]]
    assert [x['operational']['feasibility']['upper_bound'] for x in rounds]==pytest.approx([.42,.15,0])
    assert rounds[0]['operational']['joint_diagnostic']['H_mw']==pytest.approx([.42])
    assert rounds[-1]['operational']['cost']['upper_bound']==pytest.approx(.188)
    assert r['execution']['dispatch']['cost']==pytest.approx(.14)
    assert np.allclose(r['execution']['dispatch']['energy_with_initial_mwh'],[[1,.93,1]])


@pytest.mark.parametrize('dt',[.25,1,6])
@pytest.mark.parametrize('binary',[False,True])
def test_units_and_binary_equivalence(dt,binary):
    d=data();d['dt_hours']=dt
    c=Case(d);stage=compute_directrix(c);dec=make_declaration(c,stage,[3])
    r=solve_declared_dispatch(c,stage,dec['declared_p_mw'],allow_virtual=False,explicit_mip=binary)
    assert r['executable']
    assert r['objective']==pytest.approx(.12,abs=1e-7)
    assert r['dispatch']['throughput_mwh_per_storage'][0]==pytest.approx(.12*dt)
    assert r['dispatch']['energy_mwh'][0,0]==pytest.approx(1-.06*dt)


def test_no_storage_arbitrary_tiers_and_zero_radius():
    d=data();d['storage']=[];d['storage_ids']=[]
    d['incentive_prices_per_mwh']=[1,3,7]
    d['response_degrees']=[[0,.5,1]];d['deviation_factors']=[[.1,.2,0]]
    r=run_two_stage(Case(d))
    assert r['status']=='robust_optimal' and not r['executable']
    assert len(r['stage2']['rounds'])==3
    assert r['stage2']['rounds'][-1]['operational']['cost']['upper_bound']==0


def test_eta_objective_not_weighted_by_storage_cost():
    d=data();d['storage'][0]['throughput_cost_coefficient']=1e9
    c=Case(d);stage=compute_directrix(c)
    dec=make_declaration(c,stage,[3]);r=solve_declared_dispatch(c,stage,dec['declared_p_mw'])
    assert r['objective']==pytest.approx(0)
    assert r['dispatch']['cost_components']['storage_throughput']>1


def test_max_tier_failure_and_open_bounds(monkeypatch):
    d=data();d['storage']=[];d['storage_ids']=[]
    r=run_two_stage(Case(d))
    assert r['status']=='no_feasible_scheme' and len(r['stage2']['rounds'])==4
    import lcdl_algorithm as algo
    monkeypatch.setattr(algo,'declared_global_oracle',lambda *a,**k:{'zero_certified':False,'positive_certified':False})
    r=run_two_stage(Case(data()))
    assert r['status']=='unverified' and len(r['stage2']['rounds'])==1
    assert not r['robust_certified']


def test_feedback_requires_band_energy_and_physical_check():
    d=data();d['storage']=[];d['storage_ids']=[]
    c=Case(d);stage=compute_directrix(c);dec=make_declaration(c,stage,[4])
    r=dispatch_feedback(c,stage,dec,[[.41,.39]])
    assert r['feedback_check']['in_set'] and r['status']=='infeasible' and not r['executable']
    assert dispatch_feedback(c,stage,dec,[[.5,.3]])['status']=='out_of_set'
    assert dispatch_feedback(c,stage,dec,[[.4,.5]])['status']=='energy_changed'
    with pytest.raises(ModelDataError,match='saved declaration'):
        dispatch_feedback(c,stage,[4],[[.4,.4]])


def test_provider_and_pending_round():
    calls=[]
    def provider(context):
        calls.append(copy.deepcopy(context))
        if context['round_index']==1: return None
        return {'response_degree':[.2],'deviation_factor':[.05]}
    r=run_two_stage(Case(data()),declaration_provider=provider)
    assert r['status']=='awaiting_declaration' and r['stage2']['requested_levels'].tolist()==[2]
    assert calls[1]['previous_round']['upgrade_mask'].tolist()==[True]


@pytest.mark.parametrize('topology',['chain','star'])
def test_joint_H_uses_one_dispatch_and_symmetric_gap_tie_break(topology):
    d=data();d['bus_ids']=[99,70,10];d['root_bus']=10
    d['rigid_p_mw']=[[.3]*2,[.3]*2,[0]*2]
    d['fixed_q_mvar']=[[.01]*2,[.02]*2,[0]*2];d['renewable_p_mw']=[[0]*2]*3
    d['dr_pre_mw']=[[.7,.1],[.2,.2]];d['dr_allocation']=[[1,0],[0,1],[0,0]]
    d['user_ids']=['a','b'];d['response_degrees']=[[0,.5,.8,1]]*2;d['deviation_factors']=[[.1]*4]*2
    d['storage']=[];d['storage_ids']=[]
    edges=[(99,10),(70,99)] if topology=='chain' else [(70,10),(99,10)]
    d['branches']=[{'from_bus':a,'to_bus':b,'r_ohm':.01,'x_ohm':.01,'p_min_mw':-5,'p_max_mw':5} for a,b in edges]
    c=Case(d);r=run_two_stage(c)
    assert r['status']!='stalled' and not r['robust_certified']
    op=r['stage2']['rounds'][0]['operational']
    assert op['feasibility']['lower_bound']>0
    diag=op['joint_diagnostic']
    assert diag['upgrade_certified']
    assert np.all(r['stage2']['rounds'][0]['upgrade_mask'])
    assert sum(diag['H_mw'])==pytest.approx(op['feasibility']['upper_bound'],abs=1e-6)
    assert diag['H_mw'][0]==pytest.approx(diag['H_mw'][1],abs=1e-6)
    assert np.allclose(diag['H_mw'],diag['dispatch']['epsilon_mw'][c.user_nodes].sum(axis=1))
    assert 'user_diagnostics' not in op and 'H_lower_bound' not in op


def test_general_three_period_vertex_polytope():
    d=data()
    for key in ('rigid_p_mw','fixed_q_mvar','renewable_p_mw'):
        d[key]=[row+[row[-1]] for row in d[key]]
    d['dr_pre_mw']=[[.4,.4,.4]]
    c=Case(d);stage=compute_directrix(c);dec=make_declaration(c,stage,[4])
    u=DeclaredUncertaintySet(c,dec);vertices=u.vertices(256)
    assert len(vertices)==6
    assert all(u.contains(x) for x in vertices)
    r=declared_global_oracle(c,stage,dec)
    assert r['zero_certified'] and r['solver']['evaluated_vertex_count']==6


def test_affine_rhs_and_full_matrix_exports(tmp_path):
    c=Case(data());r=run_two_stage(c,feedback_p_mw=[[.47,.33]])
    p=export_result(r,tmp_path/'result.json',True);saved=json.loads(p.read_text(encoding='utf-8'))
    assert len(saved['matrix_exports'])==5
    for row in saved['stage2']['rounds']:
        diag=row['operational'].get('joint_diagnostic')
        if diag is not None:
            meta=diag['model_export']
            folder=tmp_path/saved['matrix_exports'][meta['base_matrix_key']]['directory']
            a,b,cc,f=[sparse.load_npz(folder/(key+'.npz')) for key in ('A','B','C','F')]
            v=np.load(folder/'vectors.npz');y=np.array(diag['y']);delta=np.array(diag['delta_mw']).ravel()
            assert np.max(a@y-v['b']-b@delta)<1e-6
            assert np.max(np.abs(cc@y-v['e']-f@delta))<1e-6
            assert v['c']@y==pytest.approx(meta['primary_objective'],abs=1e-6)
            assert np.square(diag['dispatch']['epsilon_mw']).sum()==pytest.approx(diag['secondary_objective'])
    for i,row in enumerate(r['stage2']['rounds']):
        for suffix,oracle in [('eta',row['operational']['feasibility'])]+([('cost',row['operational']['cost'])] if 'cost' in row['operational'] else []):
            folder=tmp_path/saved['matrix_exports'][f'round_{i:03d}_{suffix}']['directory']
            a,b,cc,f=[sparse.load_npz(folder/(key+'.npz')) for key in ('A','B','C','F')]
            v=np.load(folder/'vectors.npz')
            for result in oracle['scenario_results']:
                y=result['y'];delta=result['delta_mw'].ravel()
                assert np.max(a@y-v['b']-b@delta)<1e-7
                assert np.max(np.abs(cc@y-v['e']-f@delta))<1e-7
                assert v['c']@y==pytest.approx(result['objective'])


def test_invalid_parameters_and_migration():
    d=data();d['response_degrees']=[[0,.8,.4,1]]
    with pytest.raises(ModelDataError,match='nondecreasing'): Case(d)
    d=data();d['slack_penalty']=100
    with pytest.raises(ModelDataError,match='separates'): Case(d)
    d=data();d['grid_max_mw']=.2
    assert run_two_stage(Case(d))['status']=='stage1_failed'
    from test_algorithm import tiny_data
    source=tiny_data();migrated=with_declaration_model(source,[[0,.5,.8,1]],[[0]*4],[])
    assert source['schema_version']==1 and migrated['schema_version']==3
    assert run_two_stage(Case(migrated))['robust_certified']


def test_optional_cost_is_separate_from_robust_feasibility():
    r=run_two_stage(Case(data()),SolverOptions(optimize_robust_cost=False))
    assert r['status']=='robust_feasible' and r['robust_certified']
    assert not r['cost_optimality_certified'] and not r['executable']


def test_global_dual_matches_exact_positive_eta():
    c=Case(data());d=compute_directrix(c);dec=make_declaration(c,d,[1])
    oracle=declared_global_oracle(c,d,dec,SolverOptions(time_limit=5,declaration_vertex_limit=0))
    assert oracle['positive_certified'] and oracle['optimality_certified']
    assert oracle['lower_bound']==pytest.approx(.42,abs=1e-6)
    assert oracle['upper_bound']==pytest.approx(.42,abs=1e-6)


def test_positive_witness_without_closed_worst_bound_does_not_upgrade(monkeypatch):
    import lcdl_algorithm as algo
    original=algo.declared_global_oracle
    def open_bound(*args,**kwargs):
        value=original(*args,**kwargs)
        value.update(upper_bound=value['lower_bound']+1,optimality_certified=False)
        return value
    monkeypatch.setattr(algo,'declared_global_oracle',open_bound)
    r=run_two_stage(Case(data()))
    assert r['status']=='unverified' and len(r['stage2']['rounds'])==1
    op=r['stage2']['rounds'][0]['operational']
    assert op['feasibility']['positive_certified']
    assert op['joint_diagnostic']['selection_certified']
    assert not op['joint_diagnostic']['upgrade_certified']
    assert 'upgrade_mask' not in r['stage2']['rounds'][0]


def test_joint_tie_break_failure_preserves_diagnostic_without_upgrading(monkeypatch):
    import lcdl_algorithm as algo
    monkeypatch.setattr(algo,'joint_gap_diagnostic',lambda *a,**k:{'selection_certified':False,'status':'unverified'})
    r=run_two_stage(Case(data()))
    assert r['status']=='unverified' and not r['executable']
    assert r['stage2']['rounds'][0]['operational']['status']=='proven_not_robust'
    assert 'upgrade_mask' not in r['stage2']['rounds'][0]
