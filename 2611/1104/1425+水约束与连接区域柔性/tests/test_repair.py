"""Numerical regression tests: gradients, controlled DOF, original replay."""
import sys
from pathlib import Path
import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'source'))
from optimizer import repair, port_neighborhood, water_energy_gradient
from geometry_constraints_20261103.repair import repair as original_repair
from study_common import CONFIG, load


def fixture(tmp_path):
    m=Chem.AddHs(Chem.MolFromSmiles('CCOCC'))
    assert AllChem.EmbedMolecule(m,randomSeed=42)==0
    AllChem.MMFFOptimizeMolecule(m)
    m=Chem.RemoveHs(m)
    p=tmp_path/'environment.pdb'
    p.write_text('ATOM      1  CA  ALA A   1      25.000  25.000  25.000  1.00 20.00           C\nEND\n')
    return m,[0,0,0,1,1],p


def test_water_gradient_smooth_inside_and_outside():
    xyz=np.array([[0.,0.,0.],[4.,2.,1.]])
    wc=np.array([[1.,.5,.4],[-1.,2.,1.]])
    limits=np.full((2,2),2.5)
    e,g=water_energy_gradient(xyz,wc,limits,200.,1e-6)
    assert e>0
    finite=np.zeros_like(g)
    for i in range(xyz.size):
        d=np.zeros(xyz.size);d[i]=1e-6;d=d.reshape(xyz.shape)
        finite.ravel()[i]=(water_energy_gradient(xyz+d,wc,limits,200.,1e-6)[0]-water_energy_gradient(xyz-d,wc,limits,200.,1e-6)[0])/2e-6
    np.testing.assert_allclose(g,finite,rtol=1e-6,atol=1e-5)
    _,g0=water_energy_gradient(np.zeros((1,3)),np.zeros((1,3)),np.ones((1,1)),200.,1e-6)
    assert np.isfinite(g0).all()


@pytest.mark.parametrize('mode',['water','local'])
def test_full_objective_gradient(tmp_path,mode):
    m,f,p=fixture(tmp_path)
    cfg=load(CONFIG)
    w=dict(cfg['water'],atoms=[dict(xyz=(m.GetConformer().GetPositions()[0]+[1.,.2,.1]).tolist(),radius_A=1.55,water_id='synthetic')])
    probes=[]
    def probe(obj,v,bounds):
        rng=np.random.default_rng(42)
        v=v+rng.normal(0,.005,v.shape)
        e,g=obj(v)
        finite=np.zeros_like(g)
        for i in range(len(v)):
            dv=np.zeros_like(v);dv[i]=1e-6
            finite[i]=(obj(v+dv)[0]-obj(v-dv)[0])/2e-6
        np.testing.assert_allclose(g,finite,rtol=2e-4,atol=2e-3)
        probes.append(True)
    repair(m,f,p,params=dict(max_iterations=2),water=w if mode=='water' else None,local=cfg['local'] if mode=='local' else None,objective_probe=probe)
    assert probes==[True]


def test_disabled_extension_matches_original(tmp_path):
    m,f,p=fixture(tmp_path)
    x=m.GetConformer().GetPositions().copy()
    q,a=original_repair(m,f,p,params=dict(max_iterations=20))
    z,b=repair(m,f,p,params=dict(max_iterations=20))
    np.testing.assert_array_equal(q.GetConformer().GetPositions(),z.GetConformer().GetPositions())
    assert a['final_objective']==b['final_objective']
    np.testing.assert_array_equal(m.GetConformer().GetPositions(),x)


def test_port_selection_and_displacement_bound(tmp_path):
    m,f,p=fixture(tmp_path)
    assert port_neighborhood(m,f)==[1,2,3,4]
    q,a=repair(m,f,p,params=dict(max_iterations=20),local=load(CONFIG)['local'])
    assert a['max_local_displacement_A']<=.10+1e-12
    assert a['fixed_region_pair_distance_error_A']<1e-8
    assert a['local_enabled'] and not a['water_enabled']
    assert [(b.GetBeginAtomIdx(),b.GetEndAtomIdx(),str(b.GetBondType())) for b in q.GetBonds()]==[(b.GetBeginAtomIdx(),b.GetEndAtomIdx(),str(b.GetBondType())) for b in m.GetBonds()]
