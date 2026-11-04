from study_common import *

if __name__=='__main__':
    target=OUT/'execution_lock.json'
    assert not target.exists(), 'Do not overwrite an execution lock'
    write_json(target,dict(fixed_parameters_sha256=sha256(CONFIG),numerical_sources={str(TASK/f):sha256(TASK/f) for f in ['source/optimizer.py','source/run_experiments.py','source/audit.py','tests/test_repair.py']}))
