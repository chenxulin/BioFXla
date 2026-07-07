# BioFXIa

FXIa small-molecule scaffold discovery workflow.

## Environment

Run commands from the preconfigured `aidd` environment. Python-side package metadata is kept in `pyproject.toml`, while external executables such as GNINA and REINVENT4 are provided by the local environment/tool checkout.

## Pipeline

1. Collect and standardize public FXIa data.
2. Generate molecules with Mol2Mol, LibInvent, LinkInvent, and REINVENT de novo.
3. Filter generated molecules with RDKit structure standardization, Lipinski Rule-of-Five, clinical-candidate similarity, and patent-proxy similarity.
4. Merge filtered molecules into a de-duplicated candidate pool.
5. Dock candidate molecules with GNINA and merge all pose SDF scores into ranked hit tables.
6. Keep 500 docking top hits, run ADMET-AI on all 500, and rank by GNINA, ADMET, drug-likeness, similarity, and scaffold diversity.
7. Advance the final 50 MD candidates to pose review, MD setup, and later FEP.

## Current Stage

Implemented:

- FXIa public seed table construction.
- REINVENT4 input preparation.
- REINVENT4 scale-up sampling helpers.
- RDKit property and Lipinski Rule-of-Five filtering.
- Multi-model candidate-pool merging.
- GNINA chunk preparation, GPU chunk runner, score merging, and top-hit extraction.
- Docking top-hit pose QC, Murcko/ECFP clustering, and ADMET-AI input generation.
- ADMET-AI prediction integration and final 50-candidate MD ranking.

Not yet implemented:

- MD setup/analysis.
- FEP setup.

## Main Commands

Build processed seed tables:

```bash
mamba run -n aidd python scripts/data/build_fxia_seed_tables.py
```

Prepare REINVENT4 inputs:

```bash
mamba run -n aidd python scripts/generation/prepare_reinvent4_inputs.py --output-dir results/reinvent4_fxia_pilot
```

Run or resume REINVENT4 scale-up:

```bash
mamba run -n aidd bash scripts/generation/run_reinvent4_scaleup.sh
```

Filter one generated CSV:

```bash
mamba run -n aidd python scripts/filtering/filter_generated_molecules.py \
  --input results/reinvent4_fxia_pilot/raw_scale/mol2mol_scale.csv \
  --manifest results/reinvent4_fxia_pilot/generation_input_manifest.json \
  --output-dir results/reinvent4_fxia_pilot/filtered_scale/mol2mol_scale \
  --output-prefix mol2mol_scale
```

Merge filtered candidate pools:

```bash
mamba run -n aidd python scripts/filtering/merge_generation_candidates.py \
  --input mol2mol=results/reinvent4_fxia_pilot/filtered_scale/mol2mol_scale/mol2mol_scale_filtered.csv \
  --output results/candidates/fxia_generation_pool_filtered.csv
```

Prepare all filtered candidates for GNINA chunked docking:

```bash
mamba run -n aidd python scripts/docking/prepare_gnina_ligand_chunks.py \
  --candidate-csv results/candidates/fxia_generation_pool_filtered.csv \
  --output-dir results/docking/gnina_7mbo/inputs/all_filtered_chunks_full
```

Run GNINA chunks on available GPUs:

```bash
mamba run -n aidd python scripts/docking/run_gnina_gpu_chunks.py \
  --chunks-dir results/docking/gnina_7mbo/inputs/all_filtered_chunks_full \
  --output-dir results/docking/gnina_7mbo/all_filtered_gpu_cnnfast
```

Merge GNINA pose SDF scores and select 500 docking top hits:

```bash
mamba run -n aidd python scripts/docking/merge_gnina_scores.py \
  --poses-dir results/docking/gnina_7mbo/all_filtered_gpu_cnnfast/poses \
  --candidate-csv results/candidates/fxia_generation_pool_filtered.csv \
  --ligand-manifest-csv results/docking/gnina_7mbo/inputs/all_filtered_chunks_full/all_filtered_ligand_manifest.csv \
  --output-csv results/docking/gnina_7mbo/all_filtered_gpu_cnnfast/fxia_all_filtered_gnina_ranked.csv \
  --top-hits-csv results/docking/gnina_7mbo/all_filtered_gpu_cnnfast/fxia_all_filtered_gnina_top_hits.csv \
  --summary-json results/docking/gnina_7mbo/all_filtered_gpu_cnnfast/fxia_all_filtered_gnina_summary.json \
  --top-n 500
```

Cluster the 500 docking top hits and prepare ADMET-AI input:

```bash
mamba run -n aidd python scripts/docking/select_docking_hits.py \
  --input-csv results/docking/gnina_7mbo/all_filtered_gpu_cnnfast/fxia_all_filtered_gnina_top_hits.csv \
  --output-dir results/docking/gnina_7mbo/all_filtered_gpu_cnnfast/top_hit_selection
```

Run ADMET-AI on all 500 top hits:

```bash
mamba run -n aidd admet_predict \
  --data_path results/docking/gnina_7mbo/all_filtered_gpu_cnnfast/top_hit_selection/top_hits_admet_ai_input.csv \
  --save_path results/admet/admet_ai/top_hits_admet_ai_predictions.csv \
  --smiles_column smiles \
  --num_workers 8
```

Rank by GNINA, ADMET-AI, drug-likeness, similarity, and scaffold diversity; export 50 MD candidates:

```bash
mamba run -n aidd python scripts/admet/rank_admet_docking_hits.py \
  --docking-csv results/docking/gnina_7mbo/all_filtered_gpu_cnnfast/top_hit_selection/top_hits_clustered.csv \
  --admet-csv results/admet/admet_ai/top_hits_admet_ai_predictions.csv \
  --output-dir results/admet/admet_ai/top_hits_ranked \
  --top-n 50 \
  --write-sdf
```

## Notes

- Clinical candidates are calibration/exclusion anchors only, not direct design templates.
- Patent-proxy filtering is technical triage, not a legal freedom-to-operate conclusion.
- Docking, MD, ADMET, and FEP outputs must remain clearly separated from source-backed activity data.
