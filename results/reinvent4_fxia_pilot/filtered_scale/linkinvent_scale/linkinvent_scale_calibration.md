# FXIa Mol2Mol Pilot Calibration

- Total REINVENT rows: 50878
- Valid generated rows: 50878
- Unique parent InChIKeys: 50807
- Rows passing property filters: 25161
- Rows passing clinical similarity hard filter: 50878
- Rows passing current patent-proxy soft filter: 50878
- Rows passing property + clinical hard filters: 20640
- Rows passing all current filters including patent-proxy soft filter: 20640

## Recommended Gate For First Scale-Up

- MW: 350.0-650.0
- TPSA: 60.0-140.0
- cLogP: 1.0-5.0
- Formal charge: -1 to 1
- SA score: <= 4.0
- ECFP4 Tanimoto to clinical anchors: < 0.45
- ECFP4 Tanimoto to curated patent examples when available: < 0.5
- Current seed-derived patent-proxy soft triage: < 0.85

## Observed Pilot Distributions

| metric | min | p10 | median | p90 | max |
|---|---:|---:|---:|---:|---:|
| mw | 350.169 | 374.441 | 438.553 | 525.535 | 648.562 |
| tpsa | 60.04 | 76.64 | 103.85 | 129.11 | 139.96 |
| clogp | 1.0 | 1.652 | 2.95 | 4.307 | 5.0 |
| formal_charge | -1.0 | 0.0 | 0.0 | 0.0 | 1.0 |
| sa_score | 2.096 | 2.525 | 3.071 | 3.726 | 4.0 |
| lipinski_rule_of_five_violations | 0.0 | 0.0 | 0.0 | 1.0 | 2.0 |
| max_clinical_similarity | 0.066 | 0.111 | 0.143 | 0.18 | 0.284 |
| max_patent_proxy_similarity | 0.154 | 0.214 | 0.271 | 0.343 | 0.525 |

## Caveats

- Patent-proxy similarity currently uses a broad public seed proxy and is a soft triage field, not a legal/IP rejection rule.
- These filters are a first scale-up gate, not evidence of FXIa activity, selectivity, synthesizability, or freedom to operate.
- Clinical candidates remain calibration/exclusion anchors only.
