from __future__ import annotations

from inherited.constants import (
    DEFAULT_AB,
    DEFAULT_DP,
    DEFAULT_GQ,
)


def _parse_int_field(value: str) -> int | None:
    if value in ("", "."):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _parse_ad_field(ad: str) -> list[int] | None:
    if not ad or ad == ".":
        return None
    values: list[int] = []
    for part in ad.split(","):
        parsed = _parse_int_field(part)
        if parsed is None:
            return None
        values.append(parsed)
    return values or None


def is_good(
    gt: str,
    dp: str,
    ad: str,
    sb: str,
    gq: str,
    alt_index: int = 1,
) -> bool:
    if "." in gt:
        return False

    dp_value = _parse_int_field(dp)
    if dp_value is None or dp_value < DEFAULT_DP:
        return False

    gq_value = _parse_int_field(gq)
    if gq_value is None or gq_value < DEFAULT_GQ:
        return False

    ads = _parse_ad_field(ad)
    if ads is None or len(ads) <= alt_index:
        return False
    if sum(ads) == 0:
        return False
    if ads[alt_index] / sum(ads) < DEFAULT_AB:
        return False
    return True


def get_good_site(sample_field: str, alt_index: int = 1) -> tuple[int, str, str]:
    """Return allele count, GT, and GQ for one alternate allele.

    Returns ``-1`` when the site fails quality filters or FORMAT is incomplete.
    Returns ``0`` for a good-quality homozygous-reference call.
    """
    parts = sample_field.split(":")
    if len(parts) < 6:
        return -1, ".", "0"

    gt, dp, ad, sb, gq = parts[0], parts[1], parts[2], parts[3], parts[4]
    if is_good(gt, dp, ad, sb, gq, alt_index):
        ac = sum(int(g) == alt_index for g in gt.replace("|", "/").split("/"))
        return ac, gt, gq
    return -1, ".", "0"
