from __future__ import annotations

from inherited.constants import (
    DEFAULT_AB,
    DEFAULT_DP,
    DEFAULT_GQ,
)


def _parse_int_field(value: str) -> int | None:
    value = value.strip()
    if value in ("", "."):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _parse_ad_field(ad: str, clean_missing_as_zero: bool = False) -> list[int] | None:
    if not ad or ad.strip() == ".":
        return None

    values: list[int] = []
    for part in ad.split(","):
        part = part.strip()
        if part in ("", "."):
            if clean_missing_as_zero:
                values.append(0)
                continue
            return None
        try:
            values.append(int(part))
        except ValueError:
            return None
    return values or None


def is_good(
    gt: str,
    dp: str,
    ad: str,
    sb: str,
    gq: str,
    alt_index: int = 1,
    *,
    clean_missing_ad_as_zero: bool = False,
) -> bool:
    if "." in gt:
        return False

    dp_value = _parse_int_field(dp)
    if dp_value is None or dp_value < DEFAULT_DP:
        return False

    gq_value = _parse_int_field(gq)
    if gq_value is None or gq_value < DEFAULT_GQ:
        return False

    ads = _parse_ad_field(ad, clean_missing_as_zero=clean_missing_ad_as_zero)
    if ads is None or len(ads) <= alt_index:
        return False
    if sum(ads) == 0:
        return False
    if ads[alt_index] / sum(ads) < DEFAULT_AB:
        return False
    return True


def get_good_site(
    sample_field: str,
    alt_index: int = 1,
    *,
    clean_ad: bool = False,
) -> tuple[int, str, str]:
    """Return allele count, GT, and GQ for one alternate allele.

    Returns ``-1`` when the site fails quality filters or FORMAT is incomplete.
    Returns ``0`` for a good-quality homozygous-reference call.

    When ``clean_ad`` is True, missing AD values (``.``) are treated as zero depth
    for that allele. Used on the multiallelic path so one missing ALT depth does
    not reject all alleles at the site.
    """
    parts = sample_field.split(":")
    if len(parts) < 6:
        return -1, ".", "0"

    gt, dp, ad, sb, gq = parts[0], parts[1], parts[2], parts[3], parts[4]
    if is_good(
        gt,
        dp,
        ad,
        sb,
        gq,
        alt_index,
        clean_missing_ad_as_zero=clean_ad,
    ):
        ac = sum(int(g) == alt_index for g in gt.replace("|", "/").split("/"))
        return ac, gt, gq
    return -1, ".", "0"
