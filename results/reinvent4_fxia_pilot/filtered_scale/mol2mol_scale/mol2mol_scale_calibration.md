# FXIa Mol2Mol Pilot Calibration

- Total REINVENT rows: 73853
- Valid generated rows: 73853
- Unique parent InChIKeys: 73077
- Rows passing property filters: 4834
- Rows passing clinical similarity hard filter: 67428
- Rows passing current patent-proxy soft filter: 67317
- Rows passing property + clinical hard filters: 793
- Rows passing all current filters including patent-proxy soft filter: 772

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
| mw | 351.193 | 488.954 | 532.988 | 585.068 | 649.917 |
| tpsa | 61.68 | 103.85 | 125.94 | 136.94 | 139.95 |
| clogp | 1.268 | 3.428 | 4.253 | 4.837 | 5.0 |
| formal_charge | -1.0 | 0.0 | 0.0 | 0.0 | 1.0 |
| sa_score | 2.194 | 3.356 | 3.597 | 3.904 | 4.0 |
| lipinski_rule_of_five_violations | 0.0 | 0.0 | 1.0 | 1.0 | 2.0 |
| max_clinical_similarity | 0.118 | 0.16 | 0.256 | 0.368 | 0.525 |
| max_patent_proxy_similarity | 0.333 | 0.604 | 0.732 | 0.833 | 1.0 |

## Caveats

- Patent-proxy similarity currently uses a broad public seed proxy and is a soft triage field, not a legal/IP rejection rule.
- These filters are a first scale-up gate, not evidence of FXIa activity, selectivity, synthesizability, or freedom to operate.
- Clinical candidates remain calibration/exclusion anchors only.
