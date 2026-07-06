import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.generation.prepare_reinvent4_inputs import prepare_generation_inputs


def test_prepare_generation_inputs_creates_required_seed_package(tmp_path):
    summary = prepare_generation_inputs(output_dir=tmp_path)

    assert 20 <= summary["mol2mol_seed_count"] <= 50
    assert summary["clinical_exclusion_count"] >= 2
    assert summary["libinvent_template_count"] >= 6

    for relative in [
        "mol2mol/fxia_mol2mol_pilot.smi",
        "libinvent/fxia_libinvent_scaffolds.smi",
        "exclusions/fxia_clinical_similarity_exclusion.smi",
        "exclusions/fxia_patent_similarity_exclusion.smi",
        "configs/mol2mol_pilot.toml",
        "configs/mol2mol_scale.toml",
        "configs/libinvent_scale.toml",
        "configs/linkinvent_scale.toml",
        "configs/reinvent_denovo_scale.toml",
        "generation_input_manifest.json",
    ]:
        assert (tmp_path / relative).exists(), relative

    mol2mol_text = (tmp_path / "mol2mol/fxia_mol2mol_pilot.smi").read_text()
    assert "CHEMBL5427131" not in mol2mol_text
    assert "CHEMBL4112929" not in mol2mol_text

    config_text = (tmp_path / "configs/mol2mol_pilot.toml").read_text()
    assert 'device = "cpu"' in config_text
    assert "pubchem_ecfp4_with_count_with_rank_reinvent4_dict_voc.prior" in config_text

    denovo_config_text = (tmp_path / "configs/reinvent_denovo_scale.toml").read_text()
    assert "reinvent_pubchem.prior" in denovo_config_text
