#!/usr/bin/env bash
set -euo pipefail
study_task_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
study_root="$(cd -- "$study_task_dir/../../.." && pwd)"
study_output="${1:-$study_task_dir/replications/$(date -u +%Y%m%dT%H%M%SZ)}"
if [[ -e "$study_output/fixed_parameters.json" ]]; then
    echo "Output already contains a frozen experiment: $study_output; choose a new directory." >&2
    exit 1
fi
mkdir -p -- "$study_output"
export WATER_PORT_STUDY_OUTPUT="$(cd -- "$study_output" && pwd)"
study_output="$WATER_PORT_STUDY_OUTPUT"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$study_root/.quality_deps:$study_root/../generation_experiment/vina_dependencies:$study_root"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
study_python="$study_root/../.venv-cutoff-train/bin/python"
cd -- "$study_root"
"$study_python" -B "$study_task_dir/source/audit.py" > "$study_output/input_audit.log" 2>&1
"$study_python" -B "$study_task_dir/source/freeze_execution.py"
"$study_python" -B -m pytest -q -p no:cacheprovider "$study_task_dir/tests/test_repair.py" > "$study_output/regression_tests.log" 2>&1
"$study_python" -B "$study_task_dir/source/run_experiments.py" --experiment A --workers 8 > "$study_output/experiment_A.log" 2>&1
"$study_python" -B "$study_task_dir/source/run_experiments.py" --experiment B --workers 8 > "$study_output/experiment_B.log" 2>&1
"$study_python" -B "$study_task_dir/source/analyze.py" > "$study_output/analysis.log" 2>&1
"$study_python" -B "$study_task_dir/source/verify.py" > "$study_output/verification.log" 2>&1
"$study_python" -B "$study_task_dir/source/build_reports.py" > "$study_output/report_build.log" 2>&1
cp -- "$study_task_dir/README.md" "$study_task_dir/RUN_COMMANDS.md" "$study_output/"
echo "Completed offline reproduction: $study_output"
