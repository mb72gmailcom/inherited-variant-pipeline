import json
from pathlib import Path

from inherited.analyze import analyze_vcf, save_run_params
from inherited.output import parse_trio_calls, read_result_tsv, serialize_trio_calls

FIXTURES = Path(__file__).parent / "fixtures"


def test_analyze_vcf_writes_streamed_results(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=tmp_path / "out",
        multiallelic=True,
        block_size=1,
    )

    inherited_records = read_result_tsv(tmp_path / "out" / "inherited.tsv")
    bad_records = read_result_tsv(tmp_path / "out" / "mendelian_bad.tsv")

    assert len(inherited_records) == 1
    chrom, pos, ref, alt, hits = inherited_records[0]
    assert chrom == "22"
    assert pos == "3000"
    assert ref == "A"
    assert alt == "G"
    assert hits["child1"] == ("0/1", "0/0", "0/1", "30")

    assert len(bad_records) == 1
    _, _, _, _, bad_hits = bad_records[0]
    assert bad_hits["child1"] == ("1/1", "1/1", "0/1", "30")

    assert stats.inherited_entries >= 1
    assert stats.mendelian_bad_entries >= 1

    per_variant = json.loads((tmp_path / "out" / "inherited_per_variant.json").read_text())
    per_person = json.loads((tmp_path / "out" / "inherited_per_person.json").read_text())
    bad_per_gt = json.loads((tmp_path / "out" / "mendelian_bad_per_gt.json").read_text())
    summary = json.loads((tmp_path / "out" / "stats.json").read_text())

    assert per_variant["var_inh"] == 1
    assert per_person["child1"] == 1
    assert bad_per_gt["1/1:1/1:0/1"] == 1
    assert summary["inherited_variants"] == 1
    assert summary["mendelian_bad_variants"] == 1


def test_serialize_and_parse_trio_calls():
    hits = {"child1": ("0/1", "0/0", "0/1", "30")}
    payload = serialize_trio_calls(hits)
    assert payload == "child1=0/1|0/0|0/1|30"
    assert parse_trio_calls(payload) == hits


def test_analyze_vcf_biallelic_mode(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=tmp_path / "out",
        multiallelic=False,
    )

    inherited_records = read_result_tsv(tmp_path / "out" / "inherited.tsv")
    bad_records = read_result_tsv(tmp_path / "out" / "mendelian_bad.tsv")
    assert len(inherited_records) == 1
    assert len(bad_records) == 1
    assert stats.variants_seen == 4
