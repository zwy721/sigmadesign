# SWITCH t=0.6 de novo quality evaluation

`test_prompt_used.md` is the task snapshot. `manifest.csv` contains 36 predicted molecules, 36 official true-graph full-sampling controls and 12 distinct original crystal structures. No model generation or training is performed.

`config.json` was frozen before the first batch. Numerical thresholds stay unchanged throughout interface corrections. The final typed receptor check distinguishes real protein atoms from metal cofactors, using PoseBusters' covalent-radius convention for inorganic cofactors; raw metal vdW contacts remain visible. Earlier engineering artifacts in `development/` are excluded from final results.

Reproduce (or resume from cache) from the repository root:

```bash
bash de_novo_quality_20261103/run.sh
PYTHONPATH="$PWD/.quality_deps:$PWD/../generation_experiment/vina_dependencies:$PWD" ../.venv-cutoff-train/bin/python -B de_novo_quality_20261103/verify_results.py
```

The complete report is [`results/DE_NOVO_QUALITY_RESULTS.md`](results/DE_NOVO_QUALITY_RESULTS.md), also copied to the repository root. All requested CSV/JSON/PNG outputs and 3D examples are in `results/`. Per-molecule full values, null states and errors are in `molecule_metrics.jsonl`; the CSV is a flattened view.

`training_index.json` indexes all 17,420 actual PDBBind-stage members of the legacy checkpoint with input/source hashes. Its novelty scope excludes GEOM pretraining, connection-head training and official SigmaDock training. No claim of novelty against that full union is made.

Reusable implementation, dependencies, generic manifest examples, quality rules and tests: [`../denovo_quality/README.md`](../denovo_quality/README.md). Entry point: `python evaluate_denovo.py --manifest manifest.csv --output results/`.
