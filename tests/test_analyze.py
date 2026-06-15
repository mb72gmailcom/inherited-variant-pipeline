import json
from pathlib import Path

from inherited.analyze import analyze_vcf, save_results

FIXTURES = Path(__file__).parent / "fixtures"


def test_analyze_vcf_writes_inherited_and_bad(tmp_path):
    dinh, dm_bad, stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        multiallelic=True,
    )

    assert "var_inh" in dinh
    assert dinh["var_inh"]["child1"] == ("0/1", "0/0", "0/1", "30")
    assert "var_bad" in dm_bad
    assert dm_bad["var_bad"]["child1"] == ("1/1", "1/1", "0/1", "30")
    assert "var_common" not in dinh
    assert stats.inherited_entries >= 1
    assert stats.mendelian_bad_entries >= 1

    save_results(tmp_path / "out", dinh, dm_bad, stats)
    inherited = json.loads((tmp_path / "out" / "inherited.json").read_text())
    assert "var_inh" in inherited


def test_analyze_vcf_biallelic_mode():
    dinh, dm_bad, stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        multiallelic=False,
    )

    assert "var_inh" in dinh
    assert "var_bad" in dm_bad
    assert stats.variants_seen == 4
