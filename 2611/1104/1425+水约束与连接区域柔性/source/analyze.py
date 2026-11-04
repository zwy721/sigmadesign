"""Paired descriptive statistics on all saved candidates; no outcome selection."""
from study_common import *
from collections import Counter
from rdkit import Chem, RDLogger
from denovo_quality.io import load_molecule
from geometry_constraints_20261103.stereo import derive
RDLogger.DisableLog('rdApp.*')

VARIANTS=['A0_B0','A1','B1']

def compare_stereo(before,after):
    aa={(d['type'],d['center']):d for d in before['decisions']}
    bb={(d['type'],d['center']):d for d in after['decisions']}
    flips=lost=regained=0
    for k in aa.keys()|bb.keys():
        a,b=aa.get(k,{}),bb.get(k,{})
        if a.get('reliable') and b.get('reliable') and a['geometry_CIP']!=b['geometry_CIP']:flips+=1
        lost+=bool(a.get('reliable') and not b.get('reliable'))
        regained+=bool(b.get('reliable') and not a.get('reliable'))
    return flips,lost,regained


def main():
    rows=load(OUT/'input_manifest.json')
    raw={key(r):r for r in rows}
    raw_mols={}
    raw_stereo={}
    for r in rows:
        mol,_,_=load_molecule(r)
        raw_mols[key(r)]=mol
        _,_,st=derive(mol,load(r['chemistry_spec_path']))
        raw_stereo[key(r)]=st
    metrics=jsonl(OUT/'quality_A/molecule_metrics.jsonl')+jsonl(OUT/'quality_B/molecule_metrics.jsonl')
    processing=load(OUT/'processing_A.json')+load(OUT/'processing_B.json')
    pi={(r['variant'],key(r)):r for r in processing}
    control={key(r):r for r in metrics if r['metadata']['group']=='A0_B0'}
    baseline_stereo={key(r):load(Path(r['metadata']['sdf_path']).with_name('stereo_decisions.json')) for r in control.values()}
    candidate=[];bonds=[];angles=[];torsions=[];stereos=[];water_contacts=[];water_pairs=[];displacements=[];clashes=[];failures=[]
    write_json(OUT/'raw_geometry_stereo.json',[dict(receptor_id=k[0],seed=k[1],rank=k[2],**v) for k,v in raw_stereo.items()])
    for r in metrics:
        k=key(r);variant=r['metadata']['group'];p=pi[variant,k]
        folder=Path(r['metadata']['sdf_path']).parent
        opt=load(folder/'repair.json');st=load(folder/'stereo_decisions.json')
        mol,info,fids=load_molecule(r['metadata'])
        xyz=mol.GetConformer().GetPositions();x0=raw_mols[k].GetConformer().GetPositions()
        ctrl,_,_=load_molecule(control[k]['metadata']);xc=ctrl.GetConformer().GetPositions()
        delta=xyz-x0;dd=xyz-xc
        ident=dict(variant=variant,receptor_id=k[0],seed=k[1],rank=k[2])
        g=r['modules']['geometry'];pb=r['modules']['posebusters'];protein=r['modules']['protein'];checks=r['checks']
        ws=waters_for(raw[k]);wc=np.asarray([w['xyz'] for w in ws]);wr=np.asarray([w['radius_A'] for w in ws]);pt=Chem.GetPeriodicTable();rv=np.array([pt.GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])
        distances=np.linalg.norm(xyz[:,None]-wc[None,:,:],axis=-1);relative=distances/(rv[:,None]+wr[None,:]);bad=relative<.75
        for j,w in enumerate(ws):
            i=int(np.argmin(relative[:,j]));near=int(np.argmin(distances[:,j]))
            water_contacts.append(dict(**ident,water_id=w['water_id'],prepared_atom_index=w['prepared_atom_index'],water_xyz=w['xyz'],water_radius_A=w['radius_A'],worst_ligand_radius_A=float(rv[i]),worst_radii_sum_A=float(rv[i]+wr[j]),QC_min_distance_A=float(.75*(rv[i]+wr[j])),soft_contact_limit_A=float(.80*(rv[i]+wr[j])),closest_atom=near,min_distance_A=float(distances[near,j]),worst_relative_atom=i,min_relative_distance=float(relative[i,j]),distance_fail=bool(bad[:,j].any()),clashing_pairs=int(bad[:,j].sum())))
            # Reporting subset only; optimizer always uses every catalog water.
            # Include all pairs <=6 A, plus nearest and worst-relative for every water.
            selected=set(np.where(distances[:,j]<=6.)[0].tolist())|{i,near}
            for ii in sorted(selected):
                regularized=float(np.sqrt(distances[ii,j]**2+1e-12))
                lim=float(.80*(rv[ii]+wr[j]));over=max(lim-regularized,0.)
                grad=(-400*over*(xyz[ii]-wc[j])/regularized).tolist()
                water_pairs.append(dict(**ident,water_id=w['water_id'],prepared_atom_index=w['prepared_atom_index'],ligand_atom=ii,ligand_element=mol.GetAtomWithIdx(ii).GetSymbol(),distance_A=float(distances[ii,j]),relative_distance=float(relative[ii,j]),ligand_radius_A=float(rv[ii]),water_radius_A=float(wr[j]),QC_min_distance_A=float(.75*(rv[ii]+wr[j])),soft_contact_limit_A=lim,water_energy_active=variant=='A1',diagnostic_pair_soft_energy_kcal_mol=200*over**2,diagnostic_pair_soft_gradient_kcal_mol_A=grad))
        assert bool(bad.any())==(pb['tests']['minimum_distance_to_waters'] is False)
        fr,lr,rr=compare_stereo(raw_stereo[k],st)
        fc,lc,rc=compare_stereo(baseline_stereo[k],st)
        row=dict(**ident,strict_HQ=checks.get('high_quality') is True,auxiliary_HQ_no_water=aux_hq(r),PB_dock=pb['all_pass'],PB_no_water=all(v is True for n,v in pb['tests'].items() if n not in WATER_TESTS),water_distance_fail=bool(bad.any()),water_volume_fail=pb['tests']['volume_overlap_with_waters'] is False,min_water_distance_A=float(distances.min()),min_relative_water_distance=float(relative.min()),clashing_water_pairs=int(bad.sum()),clashing_water_atoms=int(bad.any(0).sum()),water_count=len(ws),severe_internal_clashes=g['severe_internal_clash_pairs'],severe_protein_clashes=protein['severe_clash_pairs'],raw_nonwater_environment_vdw_clashes=protein['raw_all_nonwater_vdw_severe_pairs'],inorganic_covalent_clashes=protein['inorganic_cofactor_covalent_radius_severe_pairs'],abnormal_bonds=g['abnormal_bond_count'],abnormal_angles=g['abnormal_angle_count'],abnormal_amides=g['abnormal_amide_torsion_count'],abnormal_rings=g['abnormal_ring_count'],RMSD_from_raw_no_alignment_A=float(np.sqrt(np.square(delta).sum(1).mean())),max_atom_displacement_from_raw_A=float(np.linalg.norm(delta,axis=1).max()),centroid_displacement_from_raw_A=float(np.linalg.norm(xyz.mean(0)-x0.mean(0))),RMSD_from_control_no_alignment_A=float(np.sqrt(np.square(dd).sum(1).mean())),max_atom_displacement_from_control_A=float(np.linalg.norm(dd,axis=1).max()),centroid_displacement_from_control_A=float(np.linalg.norm(xyz.mean(0)-xc.mean(0))),flexible_atom_count=len(opt.get('flexible_atom_indices',[])),heavy_atom_count=len(xyz),max_local_displacement_A=opt.get('max_local_displacement_A',0.),max_intrafragment_pair_change_unrounded_A=opt['max_fragment_pair_distance_change_A'],fixed_region_pair_error_unrounded_A=opt.get('fixed_region_pair_distance_error_A',opt['max_fragment_pair_distance_change_A']),optimizer_converged=opt['success'],optimizer_termination=opt['termination'],iterations=opt['iterations'],function_evaluations=opt['function_evaluations'],hydrogen_status=opt['hydrogen_optimization_status'],repair_cpu_seconds=p['repair_cpu_seconds'],repair_wall_seconds=p['repair_wall_seconds'],QC_cpu_seconds=r['quality_cpu_seconds'],QC_wall_seconds=r.get('seconds'),stereo_reliable=st['reliable'],stereo_unavailable=st['unavailable'],stereo_reliable_flips_from_raw=fr,stereo_reliability_lost_from_raw=lr,stereo_reliability_regained_from_raw=rr,stereo_reliable_flips_from_control=fc,stereo_reliability_lost_from_control=lc,stereo_reliability_regained_from_control=rc,failed_checks=r.get('failed_checks'),unavailable_checks=r.get('unavailable_checks'),PB_failed_tests=[n for n,v in pb['tests'].items() if v is False],errors=r.get('errors'),sdf_path=r['metadata']['sdf_path'])
        for n,v in checks.items():row['QC_'+n]=v
        for n,v in pb['tests'].items():row['PB_'+n]=v
        for relation in ['cross_fragment','within_fragment']:
            bb=[b for b in g['bonds'] if b['fragment_relation']==relation]
            row[relation+'_abnormal_bonds']=sum(b['abnormal'] for b in bb)
            row[relation+'_bond_MAE_A']=float(np.mean([abs(b['deviation_A']) for b in bb])) if bb else None
            row[relation+'_bond_max_deviation_A']=max((abs(b['deviation_A']) for b in bb),default=None)
        for b in g['bonds']:bonds.append(dict(**ident,**b))
        for a in g['angles']:
            relation='within_fragment' if len({fids[i] for i in a['atoms']})==1 else 'cross_fragment'
            angles.append(dict(**ident,fragment_relation=relation,**a))
        for relation in ['cross_fragment','within_fragment']:
            aa=[a for a in angles if all(a[n]==ident[n] for n in ident) and a['fragment_relation']==relation]
            row[relation+'_abnormal_angles']=sum(a['abnormal'] for a in aa)
        for scope,records in [('amide_planarity',g['amide_torsions']),('experimental_preference_diagnostic',g['experimental_torsions'])]:
            for t in records:
                relation='within_fragment' if len({fids[i] for i in t['atoms']})==1 else 'cross_fragment'
                torsions.append(dict(**ident,scope=scope,fragment_relation=relation,**t))
        for d in st['decisions']:stereos.append(dict(**ident,**d))
        for i in range(len(xyz)):
            displacements.append(dict(**ident,atom_index=i,fragment_id=fids[i],element=mol.GetAtomWithIdx(i).GetSymbol(),raw_xyz=x0[i].tolist(),repaired_xyz=xyz[i].tolist(),displacement_from_raw_A=float(np.linalg.norm(delta[i])),displacement_from_control_A=float(np.linalg.norm(dd[i])),local_delta_A=opt.get('local_displacements_A',[[0.,0.,0.]]*len(xyz))[i]))
        for c in g['clash_pairs']:clashes.append(dict(**ident,scope='internal',**c))
        for c in protein['clash_pair_details']:clashes.append(dict(**ident,scope='environment',**c))
        candidate.append(row)
        if not row['strict_HQ'] or not row['optimizer_converged'] or row['unavailable_checks'] or row['errors'] or row['stereo_reliable_flips_from_control'] or row['stereo_reliability_lost_from_control']:
            failures.append(row)
    # Unexpected optimizer failures remain in the denominator and explicit failure list.
    for p in processing:
        if p['status']!='ok':
            row=dict(variant=p['variant'],receptor_id=key(p)[0],seed=key(p)[1],rank=key(p)[2],strict_HQ=False,auxiliary_HQ_no_water=False,processing_failed=True,error=p.get('error'))
            candidate.append(row);failures.append(row)
    candidate.sort(key=lambda r:(r['variant'],key(r)))
    savecsv(OUT/'per_candidate_results.csv',candidate)
    savecsv(OUT/'failure_cases.csv',failures)
    for filename,rr in [('bond_details.csv',bonds),('angle_details.csv',angles),('torsion_details.csv',torsions),('stereo_centers.csv',stereos),('water_contacts.csv',water_contacts),('water_pair_details.csv',water_pairs),('atom_displacements.csv',displacements),('clash_details.csv',clashes)]:savecsv(OUT/filename,rr)
    write_json(OUT/'per_candidate_results.json',candidate)
    summary=[];by_receptor=[];by_seed=[]
    def summarize(rr,ident):
        out=dict(ident,n=len(rr))
        count_fields=['strict_HQ','auxiliary_HQ_no_water','PB_dock','PB_no_water','water_distance_fail','water_volume_fail','optimizer_converged','stereo_reliable_flips_from_control','stereo_reliability_lost_from_control','stereo_reliability_regained_from_control','stereo_reliable_flips_from_raw','stereo_reliability_lost_from_raw','stereo_reliability_regained_from_raw','abnormal_bonds','abnormal_angles','abnormal_amides','abnormal_rings','cross_fragment_abnormal_bonds','within_fragment_abnormal_bonds','cross_fragment_abnormal_angles','within_fragment_abnormal_angles','severe_internal_clashes','severe_protein_clashes','repair_cpu_seconds','QC_cpu_seconds']
        for field in count_fields:out[field]=sum(r.get(field,0) or 0 for r in rr)
        for field in ['RMSD_from_raw_no_alignment_A','max_atom_displacement_from_raw_A','centroid_displacement_from_raw_A','RMSD_from_control_no_alignment_A','max_atom_displacement_from_control_A','centroid_displacement_from_control_A','max_local_displacement_A','min_relative_water_distance','cross_fragment_bond_MAE_A','within_fragment_bond_MAE_A','repair_cpu_seconds','repair_wall_seconds','QC_cpu_seconds']:
            v=[r[field] for r in rr if r.get(field) is not None]
            for label,func in [('mean',np.mean),('median',np.median),('max',np.max)]:out[field+'_'+label]=float(func(v)) if v else None
        out['candidates_with_internal_clash']=sum(r.get('severe_internal_clashes',0)>0 for r in rr)
        out['candidates_with_protein_clash']=sum(r.get('severe_protein_clashes',0)>0 for r in rr)
        for field in sorted({f for r in rr for f in r if f.startswith(('QC_','PB_'))}):
            out[field+'_pass']=sum(r.get(field) is True for r in rr)
        return out
    for v in VARIANTS:
        rr=[r for r in candidate if r['variant']==v]
        summary.append(summarize(rr,dict(variant=v)))
        for pid in sorted({r['receptor_id'] for r in rr}):
            by_receptor.append(summarize([r for r in rr if r['receptor_id']==pid],dict(variant=v,receptor_id=pid)))
        for pid,seed in sorted({(r['receptor_id'],r['seed']) for r in rr}):
            by_seed.append(summarize([r for r in rr if r['receptor_id']==pid and r['seed']==seed],dict(variant=v,receptor_id=pid,seed=seed)))
    pairs=[];regressions=[]
    cc={key(r):r for r in candidate if r['variant']=='A0_B0'}
    for v in ['A1','B1']:
        rr=[r for r in candidate if r['variant']==v]
        record=dict(variant=v,n=len(rr))
        for field in ['strict_HQ','auxiliary_HQ_no_water','QC_bond_geometry_valid','QC_angle_geometry_valid','QC_no_severe_internal_clash','QC_no_severe_protein_clash','QC_stereo_definition_complete','QC_stereo_geometry_consistent']:
            record[field+'_gained']=sum(r.get(field) is True and cc[key(r)].get(field) is not True for r in rr)
            record[field+'_lost']=sum(r.get(field) is not True and cc[key(r)].get(field) is True for r in rr)
        for r in rr:
            c=cc[key(r)]
            bad=[f for f in ['strict_HQ','auxiliary_HQ_no_water','QC_bond_geometry_valid','QC_angle_geometry_valid','QC_no_severe_internal_clash','QC_no_severe_protein_clash','QC_stereo_definition_complete','QC_stereo_geometry_consistent'] if c.get(f) is True and r.get(f) is not True]
            if r.get('stereo_reliable_flips_from_control',0):bad.append('reliable_stereo_code_changed')
            if r.get('stereo_reliability_lost_from_control',0):bad.append('stereo_reliability_lost')
            if bad:regressions.append(dict(variant=v,receptor_id=key(r)[0],seed=key(r)[1],rank=key(r)[2],regressed_checks=bad,RMSD_from_control_A=r.get('RMSD_from_control_no_alignment_A'),max_from_control_A=r.get('max_atom_displacement_from_control_A')))
        # Six receptor clusters, all 24 saved candidates travel together.
        group_delta=np.array([sum(int(r['strict_HQ'])-int(cc[key(r)]['strict_HQ']) for r in rr if r['receptor_id']==pid)/24 for pid in sorted({r['receptor_id'] for r in rr})])
        rng=np.random.default_rng(20261104)
        bs=group_delta[rng.integers(0,len(group_delta),size=(5000,len(group_delta)))].mean(1)
        record['descriptive_cluster_bootstrap_HQ_rate_delta_95_interval']=np.quantile(bs,[.025,.975]).tolist()
        record['cluster_note']='Six previously inspected receptors only; illustrative paired cluster resampling, not independent-candidate inference or generalization validation.'
        pairs.append(record)
    for filename,rr in [('summary.csv',summary),('by_receptor.csv',by_receptor),('by_seed.csv',by_seed),('paired_changes.csv',pairs),('regressions.csv',regressions)]:savecsv(OUT/filename,rr)
    write_json(OUT/'summary.json',dict(overall=summary,by_receptor=by_receptor,paired_changes=pairs,additional_neural_NFE=0))
    print('Analysis complete:',[(r['variant'],r['strict_HQ'],r['auxiliary_HQ_no_water']) for r in summary],flush=True)

if __name__=='__main__':main()
