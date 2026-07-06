import csv
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.filtering.merge_generation_candidates import merge_candidate_files


def write_filtered_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "generated_smiles",
        "canonical_smiles",
        "parent_smiles",
        "parent_inchi_key",
        "mw",
        "tpsa",
        "clogp",
        "hbd",
        "hba",
        "rotatable_bonds",
        "formal_charge",
        "sa_score",
        "lipinski_rule_of_five_violations",
        "passes_lipinski_rule_of_five",
        "max_clinical_similarity",
        "max_patent_proxy_similarity",
        "passes_all_filters",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_merge_candidate_files_deduplicates_by_parent_key_and_tracks_source(tmp_path):
    first = tmp_path / "mol2mol.csv"
    second = tmp_path / "libinvent.csv"
    write_filtered_csv(
        first,
        [
            {
                "generated_smiles": "CCO",
                "canonical_smiles": "CCO",
                "parent_smiles": "CCO",
                "parent_inchi_key": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                "mw": "46.069",
                "tpsa": "20.23",
                "clogp": "-0.001",
                "hbd": "1",
                "hba": "1",
                "rotatable_bonds": "0",
                "formal_charge": "0",
                "sa_score": "1.0",
                "lipinski_rule_of_five_violations": "0",
                "passes_lipinski_rule_of_five": "true",
                "max_clinical_similarity": "0.1",
                "max_patent_proxy_similarity": "0.2",
                "passes_all_filters": "true",
            }
        ],
    )
    write_filtered_csv(
        second,
        [
            {
                "generated_smiles": "CCO",
                "canonical_smiles": "CCO",
                "parent_smiles": "CCO",
                "parent_inchi_key": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                "mw": "46.069",
                "tpsa": "20.23",
                "clogp": "-0.001",
                "hbd": "1",
                "hba": "1",
                "rotatable_bonds": "0",
                "formal_charge": "0",
                "sa_score": "1.0",
                "lipinski_rule_of_five_violations": "0",
                "passes_lipinski_rule_of_five": "true",
                "max_clinical_similarity": "0.1",
                "max_patent_proxy_similarity": "0.2",
                "passes_all_filters": "true",
            },
            {
                "generated_smiles": "CCN",
                "canonical_smiles": "CCN",
                "parent_smiles": "CCN",
                "parent_inchi_key": "QUSNBJAOOMFDIB-UHFFFAOYSA-N",
                "mw": "45.085",
                "tpsa": "26.02",
                "clogp": "-0.035",
                "hbd": "1",
                "hba": "1",
                "rotatable_bonds": "0",
                "formal_charge": "0",
                "sa_score": "1.0",
                "lipinski_rule_of_five_violations": "0",
                "passes_lipinski_rule_of_five": "true",
                "max_clinical_similarity": "0.1",
                "max_patent_proxy_similarity": "0.2",
                "passes_all_filters": "false",
            },
        ],
    )

    output = tmp_path / "pool.csv"
    summary = merge_candidate_files(
        inputs=[("mol2mol", first), ("libinvent", second)],
        output_csv=output,
    )

    assert summary["input_rows"] == 3
    assert summary["passing_rows"] == 2
    assert summary["merged_rows"] == 1
    with output.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 1
    assert rows[0]["source_models"] == "libinvent;mol2mol"
