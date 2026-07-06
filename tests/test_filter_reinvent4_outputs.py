import csv
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.filtering.filter_generated_molecules import filter_reinvent4_outputs


def test_filter_reinvent4_outputs_writes_properties_similarity_and_summary(tmp_path):
    input_csv = tmp_path / "mol2mol.csv"
    with input_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["SMILES", "Input_SMILES", "Tanimoto", "NLL"])
        writer.writeheader()
        writer.writerow(
            {
                "SMILES": "CC(=O)Nc1ccc(C(=O)NCCN2CCOCC2)cc1",
                "Input_SMILES": "CC(=O)Nc1ccc(C(=O)NCCN2CCOCC2)cc1",
                "Tanimoto": "1.0",
                "NLL": "1.2",
            }
        )

    clinical = tmp_path / "clinical.smi"
    patent = tmp_path / "patent.smi"
    clinical.write_text("CC(=O)Nc1ccc(C(=O)NCCN2CCOCC2)cc1\n")
    patent.write_text("c1ccccc1\n")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "outputs": {
                    "clinical_exclusion": str(clinical),
                    "patent_proxy_exclusion": str(patent),
                }
            }
        )
    )

    summary = filter_reinvent4_outputs(input_csv=input_csv, manifest_path=manifest, output_dir=tmp_path / "filtered")

    assert summary["total_rows"] == 1
    assert summary["valid_rows"] == 1
    assert (tmp_path / "filtered/mol2mol_filtered.csv").exists()
    assert (tmp_path / "filtered/mol2mol_calibration.json").exists()
    assert (tmp_path / "filtered/mol2mol_calibration.md").exists()

    with (tmp_path / "filtered/mol2mol_filtered.csv").open(newline="") as fh:
        row = next(csv.DictReader(fh))
    for field in [
        "mw",
        "tpsa",
        "clogp",
        "formal_charge",
        "sa_score",
        "lipinski_rule_of_five_violations",
        "passes_lipinski_rule_of_five",
        "max_clinical_similarity",
        "max_patent_proxy_similarity",
    ]:
        assert row[field] != ""


def test_filter_reinvent4_outputs_uses_strict_lipinski_rule_of_five(tmp_path):
    input_csv = tmp_path / "mol2mol.csv"
    with input_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["SMILES", "Input_SMILES", "Tanimoto", "NLL"])
        writer.writeheader()
        writer.writerow(
            {
                "SMILES": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
                "Input_SMILES": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
                "Tanimoto": "1.0",
                "NLL": "1.2",
            }
        )

    clinical = tmp_path / "clinical.smi"
    patent = tmp_path / "patent.smi"
    clinical.write_text("c1ccncc1\n")
    patent.write_text("c1ccccc1\n")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "outputs": {
                    "clinical_exclusion": str(clinical),
                    "patent_proxy_exclusion": str(patent),
                }
            }
        )
    )

    filter_reinvent4_outputs(input_csv=input_csv, manifest_path=manifest, output_dir=tmp_path / "filtered")

    with (tmp_path / "filtered/mol2mol_filtered.csv").open(newline="") as fh:
        row = next(csv.DictReader(fh))
    assert int(row["lipinski_rule_of_five_violations"]) >= 1
    assert row["passes_lipinski_rule_of_five"] == "false"
    assert row["passes_all_filters"] == "false"


def test_filter_reinvent4_outputs_accepts_output_prefix(tmp_path):
    input_csv = tmp_path / "libinvent_scale.csv"
    with input_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["SMILES", "NLL"])
        writer.writeheader()
        writer.writerow({"SMILES": "CC(=O)Nc1ccc(C(=O)NCCN2CCOCC2)cc1", "NLL": "1.2"})

    clinical = tmp_path / "clinical.smi"
    patent = tmp_path / "patent.smi"
    clinical.write_text("c1ccncc1\n")
    patent.write_text("c1ccccc1\n")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "outputs": {
                    "clinical_exclusion": str(clinical),
                    "patent_proxy_exclusion": str(patent),
                }
            }
        )
    )

    summary = filter_reinvent4_outputs(
        input_csv=input_csv,
        manifest_path=manifest,
        output_dir=tmp_path / "filtered",
        output_prefix="libinvent_scale",
    )

    assert summary["output_prefix"] == "libinvent_scale"
    assert (tmp_path / "filtered/libinvent_scale_filtered.csv").exists()
    assert (tmp_path / "filtered/libinvent_scale_calibration.json").exists()
    assert (tmp_path / "filtered/libinvent_scale_calibration.md").exists()
