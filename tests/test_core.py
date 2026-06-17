from pathlib import Path

import pytest

from inherited.af import is_rare, load_af_json
from inherited.classify import classify_trio
from inherited.families import load_family_relations
from inherited.genotype import get_good_site, is_good


FIXTURES = Path(__file__).parent / "fixtures"


def test_load_af_json_scalar_and_object(tmp_path):
    table = load_af_json(FIXTURES / "tiny_af.json")
    assert table["var_rare"] == 0.001
    assert is_rare(table, "var_common") is False
    assert is_rare(table, "var_rare") is True


def test_load_af_json_object_value(tmp_path):
    path = tmp_path / "af.json"
    path.write_text('{"k1": {"AF": 0.2, "AF_EUR": 0.01}, "k2": {"AF": 0.005}}')
    table = load_af_json(path)
    assert table["k1"] == 0.01
    assert table["k2"] == 0.005


def test_load_family_relations():
    rel = load_family_relations(FIXTURES / "families.tsv")
    assert rel.trio_cl["child1"] == ("ma1", "fa1")
    assert rel.trios_ids == [["child1", "fa1", "ma1"]]
    assert rel.family_size["fam1"] == 3


def test_is_good_rejects_low_dp():
    assert not is_good("0/1", "5", "2,3", "0,0,0,0", "30", 1)


def test_is_good_rejects_missing_ad():
    assert not is_good("0/1", "30", ".", "0,0,0,0", "30", 1)
    assert not is_good("0/1", "30", "15,.", "0,0,0,0", "30", 1)


def test_get_good_site_handles_missing_ad():
    sample = "0/1:30:.:0,0,0,0:30:0,30,30:."
    ac, gt, gq = get_good_site(sample, 1)
    assert ac == -1


def test_get_good_site_counts_alt():
    sample = "0/1:30:15,15:0,0,0,0:30:0,30,30:."
    ac, gt, gq = get_good_site(sample, 1)
    assert ac == 1
    assert gt == "0/1"
    assert gq == "30"


def test_get_good_site_returns_negative_one_when_not_good():
    sample = "0/1:5:2,3:0,0,0,0:30:0,30,30:."
    ac, gt, gq = get_good_site(sample, 1)
    assert ac == -1
    assert gt == "."
    assert gq == "0"


@pytest.mark.parametrize(
    ("ac", "mac", "fac", "expected"),
    [
        (1, 2, 2, False),
        (2, 2, 2, True),
        (2, 1, 0, False),
        (1, 1, 0, True),
        (2, 0, 1, False),
        (1, 0, 1, True),
    ],
)
def test_classify_trio(ac, mac, fac, expected):
    assert classify_trio(ac, mac, fac) is expected
