# FXIa Mol2Mol Pilot Calibration

- Total REINVENT rows: 97283
- Valid generated rows: 97283
- Unique parent InChIKeys: 97283
- Rows passing property filters: 17847
- Rows passing clinical similarity hard filter: 97283
- Rows passing current patent-proxy soft filter: 97283
- Rows passing property + clinical hard filters: 16278
- Rows passing all current filters including patent-proxy soft filter: 16278

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
| mw | 350.127 | 361.509 | 410.43 | 493.708 | 649.498 |
| tpsa | 60.0 | 66.216 | 85.01 | 114.71 | 140.0 |
| clogp | 1.0 | 1.947 | 3.4 | 4.592 | 5.0 |
| formal_charge | -1.0 | 0.0 | 0.0 | 0.0 | 1.0 |
| sa_score | 1.649 | 2.16 | 2.738 | 3.546 | 4.0 |
| lipinski_rule_of_five_violations | 0.0 | 0.0 | 0.0 | 0.0 | 2.0 |
| max_clinical_similarity | 0.045 | 0.101 | 0.133 | 0.176 | 0.309 |
| max_patent_proxy_similarity | 0.075 | 0.141 | 0.173 | 0.213 | 0.327 |

## Caveats

- Patent-proxy similarity currently uses a broad public seed proxy and is a soft triage field, not a legal/IP rejection rule.
- These filters are a first scale-up gate, not evidence of FXIa activity, selectivity, synthesizability, or freedom to operate.
- Clinical candidates remain calibration/exclusion anchors only.
