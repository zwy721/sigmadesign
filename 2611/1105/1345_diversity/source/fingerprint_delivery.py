"""Hash delivered reports, frozen code and all generation artifacts."""
from common import *
from datetime import datetime,timezone
def main():
 paths=list((OUT/'source').glob('*.py'))+list((OUT/'main/structures').rglob('*'))+list((OUT/'pilot/structures').rglob('*'))
 paths += [p for p in OUT.iterdir() if p.suffix in ['.md','.json','.csv','.npz','.png','.svg'] and p.name!='delivery_fingerprints.json']
 data={str(p.relative_to(OUT)):dict(bytes=p.stat().st_size,sha256=sha256(p)) for p in sorted(set(paths)) if p.is_file()}
 write_json(OUT/'delivery_fingerprints.json',dict(timestamp=datetime.now(timezone.utc).isoformat(),scope='Reports, CSV/JSON, source and all pilot/main generation files; active QC logs/caches excluded.',files=data))
 print('Delivery fingerprints',len(data),flush=True)
if __name__=='__main__':main()

