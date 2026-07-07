import csv
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.admet.rank_admet_docking_hits import rank_admet_docking_hits


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def test_rank_admet_docking_hits_uses_admet_risk_and_scaffold_diversity(tmp_path):
    docking_fields = [
        "candidate_id",
        "docking_rank",
        "cnn_score",
        "minimized_affinity_kcal_mol",
        "max_clinical_similarity",
        "max_patent_proxy_similarity",
        "lipinski_rule_of_five_violations",
        "murcko_scaffold",
        "ecfp_cluster_id",
        "analysis_smiles",
        "parent_smiles",
        "pose_qc_pass",
        "pose_sdf",
    ]
    docking_rows = [
        {
            "candidate_id": "BAD-DOCKING",
            "docking_rank": "1",
            "cnn_score": "0.99",
            "minimized_affinity_kcal_mol": "-10",
            "max_clinical_similarity": "0.1",
            "max_patent_proxy_similarity": "0.1",
            "lipinski_rule_of_five_violations": "0",
            "murcko_scaffold": "scaffold-a",
            "ecfp_cluster_id": "1",
            "analysis_smiles": "CCO",
            "parent_smiles": "CCO",
            "pose_qc_pass": "true",
            "pose_sdf": "",
        },
        {
            "candidate_id": "GOOD-A1",
            "docking_rank": "2",
            "cnn_score": "0.98",
            "minimized_affinity_kcal_mol": "-9",
            "max_clinical_similarity": "0.1",
            "max_patent_proxy_similarity": "0.1",
            "lipinski_rule_of_five_violations": "0",
            "murcko_scaffold": "scaffold-a",
            "ecfp_cluster_id": "2",
            "analysis_smiles": "CCN",
            "parent_smiles": "CCN",
            "pose_qc_pass": "true",
            "pose_sdf": "",
        },
        {
            "candidate_id": "GOOD-B",
            "docking_rank": "3",
            "cnn_score": "0.97",
            "minimized_affinity_kcal_mol": "-8",
            "max_clinical_similarity": "0.1",
            "max_patent_proxy_similarity": "0.1",
            "lipinski_rule_of_five_violations": "0",
            "murcko_scaffold": "scaffold-b",
            "ecfp_cluster_id": "3",
            "analysis_smiles": "CCC",
            "parent_smiles": "CCC",
            "pose_qc_pass": "true",
            "pose_sdf": "",
        },
    ]
    admet_fields = [
        "candidate_id",
        "smiles",
        "QED",
        "HIA_Hou",
        "Bioavailability_Ma",
        "Solubility_AqSolDB",
        "Caco2_Wang",
        "hERG",
        "AMES",
        "DILI",
        "ClinTox",
        "CYP3A4_Veith",
        "CYP2C9_Veith",
        "PAINS_alert",
        "BRENK_alert",
    ]
    good_admet = {
        "QED": "0.7",
        "HIA_Hou": "0.9",
        "Bioavailability_Ma": "0.7",
        "Solubility_AqSolDB": "-3",
        "Caco2_Wang": "-4.5",
        "hERG": "0.05",
        "AMES": "0.05",
        "DILI": "0.05",
        "ClinTox": "0.05",
        "CYP3A4_Veith": "0.1",
        "CYP2C9_Veith": "0.1",
        "PAINS_alert": "0",
        "BRENK_alert": "0",
    }
    admet_rows = [
        {
            "candidate_id": "BAD-DOCKING",
            "smiles": "CCO",
            **good_admet,
            "hERG": "0.95",
            "AMES": "0.90",
            "DILI": "0.90",
        },
        {"candidate_id": "GOOD-A1", "smiles": "CCN", **good_admet},
        {"candidate_id": "GOOD-B", "smiles": "CCC", **good_admet},
    ]

    docking_csv = tmp_path / "clustered.csv"
    admet_csv = tmp_path / "admet.csv"
    output_dir = tmp_path / "ranked"
    write_csv(docking_csv, docking_fields, docking_rows)
    write_csv(admet_csv, admet_fields, admet_rows)

    summary = rank_admet_docking_hits(
        docking_csv=docking_csv,
        admet_csv=admet_csv,
        output_dir=output_dir,
        top_n=2,
        max_per_scaffold=1,
        max_per_cluster=1,
    )

    selected = read_rows(output_dir / "top_hits_md_candidates.csv")
    assert summary["selected_count"] == 2
    assert [row["candidate_id"] for row in selected] == ["GOOD-A1", "GOOD-B"]
    assert all(row["admet_risk_count"] == "0" for row in selected)
    assert len({row["murcko_scaffold"] for row in selected}) == 2
