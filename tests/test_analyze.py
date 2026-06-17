import json
from pathlib import Path

from inherited.analyze import analyze_vcf, save_results, save_run_params, summarize_inherited, summarize_mendelian_bad

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

    per_variant = json.loads((tmp_path / "out" / "inherited_per_variant.json").read_text())
    per_person = json.loads((tmp_path / "out" / "inherited_per_person.json").read_text())
    bad_per_gt = json.loads((tmp_path / "out" / "mendelian_bad_per_gt.json").read_text())
    summary = json.loads((tmp_path / "out" / "stats.json").read_text())

    assert per_variant["var_inh"] == 1
    assert per_person["child1"] == 1
    assert bad_per_gt["1/1:1/1:0/1"] == 1
    assert summary["inherited_variants"] == 1
    assert summary["mendelian_bad_variants"] == 1

    params_path = save_run_params(
        tmp_path / "out" / "chr22",
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        multiallelic=True,
        af_threshold=0.01,
    )
    assert params_path == tmp_path / "out" / "params.json"
    params = json.loads(params_path.read_text())
    assert params["multiallelic"] is True
    assert params["vcf"].endswith("tiny.vcf")
    assert "chromosome_runs" not in params


def test_summarize_inherited_and_mendelian_bad():
    dinh = {
        "v1": {"p1": ("0/1", "0/0", "0/1", "30"), "p2": ("0/1", "0/0", "0/1", "25")},
        "v2": {"p1": ("0/1", "0/1", "0/1", "30")},
    }
    dm_bad = {
        "b1": {
            "p1": ("1/1", "1/1", "0/1", "30"),
            "p2": ("1/1", "1/1", "0/1", "28"),
        }
    }

    per_variant, per_person = summarize_inherited(dinh)
    assert per_variant == {"v1": 2, "v2": 1}
    assert per_person == {"p1": 2, "p2": 1}

    per_gt = summarize_mendelian_bad(dm_bad)
    assert per_gt == {"1/1:1/1:0/1": 2}


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
