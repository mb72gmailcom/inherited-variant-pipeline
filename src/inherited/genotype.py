from __future__ import annotations

from inherited.constants import (
    DEFAULT_AB,
    DEFAULT_DP,
    DEFAULT_GQ,
)


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
    if dp == "." or int(dp) < DEFAULT_DP:
        return False
    if gq == "." or int(gq) < DEFAULT_GQ:
        return False

    ads = [int(x) for x in ad.split(",")]
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
