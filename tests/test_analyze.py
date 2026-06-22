import json
from pathlib import Path

from inherited.analyze import analyze_vcf
from inherited.output import (
    parse_trio_calls,
    read_result_tsv,
    serialize_patient_ids,
    serialize_payload,
    serialize_trio_calls,
)

FIXTURES = Path(__file__).parent / "fixtures"


def test_analyze_vcf_writes_short_format_by_default(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=tmp_path / "out",
        multiallelic=True,
        block_size=1,
    )

    inherited_records = read_result_tsv(tmp_path / "out" / "inherited.tsv", short_format=True)
    bad_records = read_result_tsv(tmp_path / "out" / "mendelian_bad.tsv", short_format=True)

    assert len(inherited_records) == 1
    chrom, pos, ref, alt, patient_ids = inherited_records[0]
    assert chrom == "22"
    assert pos == "3000"
    assert ref == "A"
    assert alt == "G"
    assert patient_ids == ["child1"]

    assert len(bad_records) == 1
    _, _, _, _, bad_patient_ids = bad_records[0]
    assert bad_patient_ids == ["child1"]

    assert stats.inherited_entries >= 1
    assert stats.mendelian_bad_entries >= 1

    header = (tmp_path / "out" / "inherited.tsv").read_text(encoding="utf-8").splitlines()[0]
    assert "PATIENTS" in header


def test_analyze_vcf_writes_full_format_when_disabled(tmp_path):
    analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=tmp_path / "out",
        short_format=False,
    )

    inherited_records = read_result_tsv(tmp_path / "out" / "inherited.tsv", short_format=False)
    chrom, pos, ref, alt, hits = inherited_records[0]
    assert hits["child1"] == ("0/0", "0/1", "0/1", "30")

    header = (tmp_path / "out" / "inherited.tsv").read_text(encoding="utf-8").splitlines()[0]
    assert "TRIO_CALLS" in header


def test_serialize_payload_short_and_full():
    hits = {"child1": ("0/0", "0/1", "0/1", "30")}
    assert serialize_payload(hits, short_format=True) == "child1"
    assert serialize_payload(hits, short_format=False) == "child1=0/0|0/1|0/1|30"
    assert serialize_patient_ids(hits) == "child1"
    assert serialize_trio_calls(hits) == "child1=0/0|0/1|0/1|30"
    assert parse_trio_calls("child1=0/0|0/1|0/1|30") == hits


def test_analyze_vcf_biallelic_mode(tmp_path):
    stats = analyze_vcf(
        vcf_path=FIXTURES / "tiny.vcf",
        af_json_path=FIXTURES / "tiny_af.json",
        family_file=FIXTURES / "families.tsv",
        output_dir=tmp_path / "out",
        multiallelic=False,
    )

    inherited_records = read_result_tsv(tmp_path / "out" / "inherited.tsv", short_format=True)
    bad_records = read_result_tsv(tmp_path / "out" / "mendelian_bad.tsv", short_format=True)
    assert len(inherited_records) == 1
    assert len(bad_records) == 1
    assert stats.variants_seen == 4

    summary = json.loads((tmp_path / "out" / "stats.json").read_text())
    assert summary["inherited_variants"] == 1
