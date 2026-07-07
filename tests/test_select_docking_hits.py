import csv
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.docking.select_docking_hits import select_docking_hits


def write_top_hits(path: Path) -> None:
    fieldnames = [
        "docking_rank",
        "candidate_id",
        "source_models",
        "minimized_affinity_kcal_mol",
        "cnn_score",
        "cnn_affinity",
        "mw",
        "tpsa",
        "clogp",
        "max_clinical_similarity",
        "max_patent_proxy_similarity",
        "parent_smiles",
        "canonical_smiles",
        "pose_sdf",
    ]
    rows = [
        ("1", "HIT-001", "reinvent_denovo", "-9.0", "0.99", "7.1", "356", "80", "2.2", "0.10", "0.20", "COc1ccc(C(=O)NCCO)cc1", "COc1ccc(C(=O)NCCO)cc1"),
        ("2", "HIT-002", "reinvent_denovo", "-8.8", "0.98", "7.0", "370", "83", "2.3", "0.11", "0.21", "COc1ccc(C(=O)NCCN)cc1", "COc1ccc(C(=O)NCCN)cc1"),
        ("3", "HIT-003", "mol2mol", "-8.4", "0.97", "6.8", "382", "90", "2.6", "0.12", "0.22", "O=C(NCCO)c1ccncc1", "O=C(NCCO)c1ccncc1"),
        ("4", "HIT-004", "libinvent", "-8.2", "0.96", "6.7", "410", "95", "2.9", "0.13", "0.23", "Cc1nnc(NC(=O)CO)s1", "Cc1nnc(NC(=O)CO)s1"),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(zip(fieldnames[:-1], row), pose_sdf=""))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def test_select_docking_hits_balances_md_subset_by_scaffold(tmp_path):
    input_csv = tmp_path / "top_hits.csv"
    output_dir = tmp_path / "selected"
    write_top_hits(input_csv)

    summary = select_docking_hits(
        input_csv=input_csv,
        output_dir=output_dir,
        admet_n=4,
        pre_admet_md_n=3,
        max_per_scaffold_admet=2,
        max_per_scaffold_md=1,
        max_per_cluster_admet=2,
        max_per_cluster_md=1,
        fingerprint_sim_threshold=0.70,
    )

    assert summary["input_rows"] == 4
    assert summary["admet_subset_count"] == 4
    assert summary["pre_admet_md_subset_count"] == 3

    annotated = read_rows(output_dir / "top_hits_clustered.csv")
    md_rows = read_rows(output_dir / "top_hits_pre_admet_md_candidates.csv")

    for field in [
        "murcko_scaffold",
        "generic_murcko_scaffold",
        "ecfp_cluster_id",
        "pose_qc_pass",
        "pose_qc_flags",
        "selected_for_admet",
        "selected_for_pre_admet_md",
    ]:
        assert field in annotated[0]

    assert md_rows[0]["candidate_id"] == "HIT-001"
    assert len({row["murcko_scaffold"] for row in md_rows}) == 3


def test_select_docking_hits_writes_full_admet_ai_input_when_admet_n_is_zero(tmp_path):
    input_csv = tmp_path / "top_hits.csv"
    output_dir = tmp_path / "selected"
    write_top_hits(input_csv)

    summary = select_docking_hits(
        input_csv=input_csv,
        output_dir=output_dir,
        admet_n=0,
        pre_admet_md_n=0,
        max_per_scaffold_md=1,
        max_per_cluster_md=1,
    )

    assert summary["admet_input_count"] == 4
    assert summary["pre_admet_md_subset_count"] == 0

    admet_rows = read_rows(output_dir / "top_hits_admet_ai_input.csv")
    assert len(admet_rows) == 4
    assert set(admet_rows[0]) >= {"candidate_id", "smiles", "docking_rank", "cnn_score"}
    assert not (output_dir / "top_hits_md_subset.csv").exists()
