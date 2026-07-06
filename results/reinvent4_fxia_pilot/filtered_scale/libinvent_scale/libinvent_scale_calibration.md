# FXIa Mol2Mol Pilot Calibration

- Total REINVENT rows: 27556
- Valid generated rows: 27556
- Unique parent InChIKeys: 27547
- Rows passing property filters: 63
- Rows passing clinical similarity hard filter: 27556
- Rows passing current patent-proxy soft filter: 27556
- Rows passing property + clinical hard filters: 0
- Rows passing all current filters including patent-proxy soft filter: 0

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
| mw | 561.69 | 589.7 | 614.71 | 633.521 | 647.736 |
| tpsa | 87.46 | 100.498 | 113.76 | 129.912 | 133.63 |
| clogp | 3.619 | 4.041 | 4.516 | 4.896 | 5.0 |
| formal_charge | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| sa_score | 3.451 | 3.577 | 3.689 | 3.84 | 3.973 |
| lipinski_rule_of_five_violations | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| max_clinical_similarity | 0.101 | 0.103 | 0.11 | 0.12 | 0.139 |
| max_patent_proxy_similarity | 0.319 | 0.333 | 0.345 | 0.357 | 0.37 |

## Caveats

- Patent-proxy similarity currently uses a broad public seed proxy and is a soft triage field, not a legal/IP rejection rule.
- These filters are a first scale-up gate, not evidence of FXIa activity, selectivity, synthesizability, or freedom to operate.
- Clinical candidates remain calibration/exclusion anchors only.
