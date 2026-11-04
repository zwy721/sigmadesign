# Adapted verbatim from the frozen original repair.py except the explicitly
# optional water term and port-neighborhood coordinates. Original file is read-only.

def port_neighborhood(mol,fragment_ids):
    endpoints=set()
    for b in mol.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        if fragment_ids[i]!=fragment_ids[j]:endpoints.update([i,j])
    selected=set(endpoints)
    for i in endpoints:selected.update(a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors())
    return sorted(selected)

def water_energy_gradient(xyz,water_xyz,limits,weight,epsilon):
    import numpy as np
    dp=xyz[:,None,:]-water_xyz[None,:,:]
    dist=np.sqrt(np.sum(dp*dp,axis=-1)+epsilon*epsilon)
    over=np.maximum(limits-dist,0.)
    energy=float(weight*np.square(over).sum())
    gradient=(-2*weight*over[:,:,None]*dp/dist[:,:,None]).sum(1)
    return energy,gradient

"""A-priori rigid-fragment MMFF94s optimization with pocket/pose restraints.

Only fragment SE(3) and added hydrogens are variables. No reference ligand is
accepted by this API. Angles are XYZ Euler parameters with analytic Jacobian.
"""
import json,time
import numpy as np
from pathlib import Path
from scipy.optimize import minimize
from scipy.spatial.distance import cdist
from rdkit import Chem
from rdkit.Chem import AllChem,rdMolTransforms
from denovo_quality.protein import read_protein

PARAMETERS=dict(protocol='a_priori_v1',forcefield='MMFF94s',pose_k=5.,orientation_k=2.,protein_k=200.,protein_radius_ratio=.80,amide_force_constant=2.,amide_window_deg=5.,max_iterations=800,hydrogen_iterations=200,translation_bound_A=2.,rotation_bound_rad=1.5707963267948966,soft_atom_displacement_A=2.,displacement_k=200.,waters_in_optimizer=False,training=False,hyperparameter_selection='Single a-priori setting; no CASF hyperparameter search')

def rotation_derivatives(v):
    a,b,c=v;sa,sb,sc=np.sin(v);ca,cb,cc=np.cos(v)
    X=np.array([[1,0,0],[0,ca,-sa],[0,sa,ca]]);Y=np.array([[cb,0,sb],[0,1,0],[-sb,0,cb]]);Z=np.array([[cc,-sc,0],[sc,cc,0],[0,0,1]])
    dX=np.array([[0,0,0],[0,-sa,-ca],[0,ca,-sa]]);dY=np.array([[-sb,0,cb],[0,0,0],[-cb,0,-sb]]);dZ=np.array([[-sc,-cc,0],[cc,-sc,0],[0,0,0]])
    return Z@Y@X,[Z@Y@dX,Z@dY@X,dZ@Y@X]

def repair(mol,fragment_ids,protein_path,params=None,water=None,local=None,objective_probe=None):
    p=dict(PARAMETERS);p.update(params or {});started=time.perf_counter();cpu_started=time.process_time()
    if fragment_ids is None:raise ValueError('Rigid repair requires original fragment partition')
    n=mol.GetNumAtoms();x0=np.array(mol.GetConformer().GetPositions());groups=[np.where(np.asarray(fragment_ids)==f)[0] for f in sorted(set(fragment_ids))];G=len(groups)
    q=Chem.AddHs(Chem.Mol(mol),addCoords=True)
    if not AllChem.MMFFHasAllMoleculeParams(q):raise ValueError('No complete MMFF94s parameterization')
    props=AllChem.MMFFGetMoleculeProperties(q,mmffVariant='MMFF94s')
    def makeff():return AllChem.MMFFGetMoleculeForceField(q,props,nonBondedThresh=100.,ignoreInterfragInteractions=False)
    ff=makeff()
    for i in range(n):ff.AddFixedPoint(i)
    hstatus=ff.Minimize(maxIts=p['hydrogen_iterations']);initial_energy=float(ff.CalcEnergy());coords0=np.array(q.GetConformer().GetPositions());ff=makeff()
    amide=[]
    for ids in mol.GetSubstructMatches(Chem.MolFromSmarts('[O,S]=[C]-[N]-[!#1]')):
        phi=rdMolTransforms.GetDihedralDeg(mol.GetConformer(),*ids);target=0. if abs(phi)<90 else (180. if phi>=0 else -180.)
        ff.MMFFAddTorsionConstraint(*ids,False,target-p['amide_window_deg'],target+p['amide_window_deg'],p['amide_force_constant']);amide.append(dict(atoms=list(ids),initial_deg=phi,target_deg=target))
    centers=[x0[g].mean(0) for g in groups];locals_=[x0[g]-centers[j] for j,g in enumerate(groups)]
    atoms=[a for a in read_protein(str(protein_path))['atoms'] if a['kind']!='water'];pc=np.array([a['xyz'] for a in atoms]);pt=Chem.GetPeriodicTable();rv=np.array([pt.GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()]);rc=np.array([pt.GetRcovalent(a.GetAtomicNum()) for a in mol.GetAtoms()]);pv=np.array([a['vdw_radius_A'] for a in atoms]);prc=np.array([a['covalent_radius_A'] for a in atoms]);inorganic=np.array([a['kind']=='inorganic_cofactor' for a in atoms])
    limits=p['protein_radius_ratio']*(rv[:,None]+pv[None,:]);limits[:,inorganic]=p['protein_radius_ratio']*(rc[:,None]+prc[None,inorganic])
    # Restrict receptor to atoms that can be reached within conservative bounds.
    keep=(cdist(x0,pc)<10.).any(0);pc=pc[keep];limits=limits[:,keep]
    nh=len(coords0)-n;h_end=G*6+nh*3
    flexible=port_neighborhood(mol,fragment_ids) if local is not None else []
    nf=len(flexible)
    if water is not None:
        wc=np.asarray([w['xyz'] for w in water['atoms']],dtype=float).reshape(-1,3)
        wr=np.asarray([w['radius_A'] for w in water['atoms']],dtype=float)
        water_limits=water['radius_ratio']*(rv[:,None]+wr[None,:])
    v0=np.concatenate([np.zeros(G*6),coords0[n:].ravel(),np.zeros(3*nf)]);evals=0
    def obj(v):
        nonlocal evals
        evals+=1;u=v[:G*6].reshape(G,6);xx=coords0.copy();derivatives=[]
        local_delta=np.zeros((n,3))
        if nf:local_delta[flexible]=v[h_end:].reshape(nf,3)
        rotations=[]
        for j,g in enumerate(groups):
            R,dR=rotation_derivatives(u[j,3:]);base=locals_[j]+local_delta[g] if nf else locals_[j]
            xx[g]=base@R.T+centers[j]+u[j,:3];derivatives.append(dR);rotations.append(R)
        xx[n:]=v[G*6:h_end].reshape(nh,3);flat=tuple(xx.ravel());e=float(ff.CalcEnergy(flat));grad=np.array(ff.CalcGrad(flat)).reshape(-1,3)
        delta=xx[:n]-x0;e+=p['pose_k']*np.square(delta).sum();grad[:n]+=2*p['pose_k']*delta
        disp=np.linalg.norm(delta,axis=1);over=np.maximum(disp-p['soft_atom_displacement_A'],0);e+=p['displacement_k']*(over**2).sum();grad[:n]+=2*p['displacement_k']*over[:,None]*delta/np.maximum(disp[:,None],1e-12)
        dp=xx[:n,None,:]-pc[None,:,:];dist=np.linalg.norm(dp,axis=-1);over=np.maximum(limits-dist,0);e+=p['protein_k']*np.square(over).sum();grad[:n]+=( -2*p['protein_k']*over[:,:,None]*dp/np.maximum(dist[:,:,None],1e-9)).sum(1)
        if water is not None:
            we,wg=water_energy_gradient(xx[:n],wc,water_limits,water['weight_kcal_mol_A2'],water['epsilon_A'])
            e+=we;grad[:n]+=wg
        gg=np.zeros_like(u);local_grad=np.zeros((n,3))
        for j,g in enumerate(groups):
            gg[j,:3]=grad[g].sum(0)
            base=locals_[j]+local_delta[g] if nf else locals_[j]
            for k in range(3):gg[j,3+k]=np.sum(grad[g]*(base@derivatives[j][k].T))
            if nf:local_grad[g]=grad[g]@rotations[j]
        e+=p['orientation_k']*np.square(u[:,3:]).sum();gg[:,3:]+=2*p['orientation_k']*u[:,3:]
        if nf:
            e+=local['position_k_kcal_mol_A2']*np.square(local_delta).sum()
            local_grad+=2*local['position_k_kcal_mol_A2']*local_delta
        return e,np.concatenate([gg.ravel(),grad[n:].ravel(),local_grad[flexible].ravel()])
    bounds=[]
    for j in range(G):bounds.extend([(-p['translation_bound_A'],p['translation_bound_A'])]*3+[(-p['rotation_bound_rad'],p['rotation_bound_rad'])]*3)
    bounds.extend([(None,None)]*(nh*3))
    if nf:bounds.extend([(-local['component_bound_A'],local['component_bound_A'])]*(3*nf))
    initial_objective=obj(v0)[0]
    if objective_probe is not None:objective_probe(obj,v0.copy(),bounds)
    result=minimize(obj,v0,method='L-BFGS-B',jac=True,bounds=bounds,options=dict(maxiter=p['max_iterations'],ftol=1e-10,gtol=1e-4,maxls=40))
    u=result.x[:G*6].reshape(G,6);xx=x0.copy()
    local_delta=np.zeros((n,3))
    if nf:local_delta[flexible]=result.x[h_end:].reshape(nf,3)
    for j,g in enumerate(groups):
        R,_=rotation_derivatives(u[j,3:]);base=locals_[j]+local_delta[g] if nf else locals_[j]
        xx[g]=base@R.T+centers[j]+u[j,:3]
    out=Chem.Mol(mol)
    for i,x in enumerate(xx):out.GetConformer().SetAtomPosition(i,x.tolist())
    maxdist=0.
    for g in groups:maxdist=max(maxdist,float(np.abs(cdist(xx[g],xx[g])-cdist(x0[g],x0[g])).max()))
    if not nf:assert maxdist<1e-8
    meta=dict(parameters=p,success=bool(result.success),termination=str(result.message),iterations=int(result.nit),function_evaluations=evals,initial_H_relaxed_MMFF_energy=initial_energy,initial_objective=initial_objective,final_objective=float(result.fun),heavy_RMSD_no_alignment_A=float(np.sqrt(np.mean(np.sum((xx-x0)**2,axis=1)))),max_atom_displacement_A=float(np.linalg.norm(xx-x0,axis=1).max()),centroid_displacement_A=float(np.linalg.norm(xx.mean(0)-x0.mean(0))),max_fragment_pair_distance_change_A=maxdist,fragment_translation_A=u[:,:3].tolist(),fragment_Euler_rotation_rad=u[:,3:].tolist(),amide_constraints=amide,hydrogen_optimization_status=int(hstatus),seconds=time.perf_counter()-started)
    fixed_pair_error=0.
    for g in groups:
        fixed=np.asarray([i for i in g if i not in flexible],dtype=int)
        if len(fixed)>1:fixed_pair_error=max(fixed_pair_error,float(np.abs(cdist(xx[fixed],xx[fixed])-cdist(x0[fixed],x0[fixed])).max()))
    assert fixed_pair_error<1e-8
    meta.update(cpu_seconds=time.process_time()-cpu_started,water_enabled=water is not None,local_enabled=local is not None,flexible_atom_indices=flexible,local_displacements_A=local_delta.tolist(),max_local_displacement_A=float(np.linalg.norm(local_delta,axis=1).max()),fixed_region_pair_distance_error_A=fixed_pair_error)
    if local is not None:meta['local_parameters']=local
    if water is not None:
        meta['water_parameters']={k:v for k,v in water.items() if k!='atoms'}
        meta['water_ids']=[w['water_id'] for w in water['atoms']]
        meta['initial_water_energy']=water_energy_gradient(x0,wc,water_limits,water['weight_kcal_mol_A2'],water['epsilon_A'])[0]
        meta['final_water_energy']=water_energy_gradient(xx,wc,water_limits,water['weight_kcal_mol_A2'],water['epsilon_A'])[0]
    return out,meta
