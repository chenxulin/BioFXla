# FXIa Molecule Generation And Selection Workflow

This note describes the current executable workflow from generated molecules to MD candidates.

## Canonical Flow

1. Generate molecules with Mol2Mol, LibInvent, LinkInvent, and REINVENT de novo.
2. Filter each generation output with RDKit standardization, Lipinski Rule-of-Five, property limits, clinical-candidate similarity, and patent-proxy similarity.
3. Merge passing molecules into `results/candidates/fxia_generation_pool_filtered.csv`.
4. Prepare all filtered candidates as GNINA ligand chunks.
5. Run GNINA on all chunks and merge pose SDF scores into full ranked docking CSVs.
6. Keep 500 docking top hits after similarity filtering.
7. Cluster/QC those 500 hits and export a compact ADMET-AI input table.
8. Run ADMET-AI on all 500 hits.
9. Rank by GNINA score, ADMET-AI, physicochemical properties, similarity, and scaffold diversity.
10. Export 50 final MD candidates.

## Main Output Files

- Docking full ranking: `results/docking/gnina_7mbo/all_filtered_gpu_cnnfast/fxia_all_filtered_gnina_ranked.csv`
- Docking top hits: `results/docking/gnina_7mbo/all_filtered_gpu_cnnfast/fxia_all_filtered_gnina_top_hits.csv`
- Clustered 500-hit table: `results/docking/gnina_7mbo/all_filtered_gpu_cnnfast/top_hit_selection/top_hits_clustered.csv`
- ADMET-AI input: `results/docking/gnina_7mbo/all_filtered_gpu_cnnfast/top_hit_selection/top_hits_admet_ai_input.csv`
- ADMET-AI predictions: `results/admet/admet_ai/top_hits_admet_ai_predictions.csv`
- ADMET-ranked table: `results/admet/admet_ai/top_hits_ranked/top_hits_admet_ranked.csv`
- Final MD candidate CSV: `results/admet/admet_ai/top_hits_ranked/top_hits_md_candidates.csv`
- Final MD candidate SDF: `results/admet/admet_ai/top_hits_ranked/top_hits_md_candidates.sdf`

## Interpretation Notes

- Clinical candidates remain calibration and exclusion anchors, not direct templates.
- Patent-proxy similarity is technical triage, not a freedom-to-operate conclusion.
- ADMET-AI predictions are in silico prioritization signals; they do not establish safety or developability.
- The 50 MD candidates are computational hypotheses and still require pose review before MD setup.
