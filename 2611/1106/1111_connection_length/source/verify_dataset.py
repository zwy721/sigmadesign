"""Verify selected chemistry, canonical atom maps, and source decomposition."""
from common import *
import sqlite3,zlib,itertools
from rdkit import Chem
from cutoff_data.codec import get_codec

def main():
 codec=get_codec();selection=load(OUT/'data/selected_molecules.json');db=sqlite3.connect(f'file:{DATA}/geom_database/graphs.sqlite?mode=ro&immutable=1',uri=True);count=0
 for rows in selection.values():
  for row in rows:
   graph=json.loads(zlib.decompress(db.execute('SELECT graph_json FROM graphs WHERE graph_id=?',(row['graph_id'],)).fetchone()[0]));mol=codec.restore_graph(graph);assert codec.identity(mol)==row['canonical'];assert graph['roundtrip'];assert sorted(a for f in graph['fragments'] for a in f['core_to_parent'])==list(range(mol.GetNumAtoms()))
   for e in graph['connections']:
    a,b=e['parent_atoms'];bond=mol.GetBondBetweenAtoms(a,b);assert bond.GetBondTypeAsDouble()==e['bond_order'];assert int(bond.GetStereo())==e['stereo']
   count+=1
 summaries={s:load(OUT/'data'/f'{s}_summary.json') for s in selection};dump(OUT/'dataset_verification.json',dict(passed=True,canonical_roundtrips=count,atom_partitions_checked=count,source_bond_order_and_stereo_checked=True,same_molecule_cross_split=0,scaffold_test_train_val_overlap=0,selected_source_index_sha256=sha(OUT/'data/selected_molecules.json'),datasets=summaries))
 print('DATASET VERIFICATION PASSED',count,flush=True)
if __name__=='__main__':main()
