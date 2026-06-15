from __future__ import annotations


def classify_trio(ac: int, mac: int, fac: int) -> bool | None:
    """Classify a trio for one alternate allele.

    Returns:
        True  -> inherited bucket (dinh)
        False -> mendelian bad bucket (dm_bad)
        None  -> child has no alternate allele call (skip)
    """
    if ac <= 0:
        return None

    if mac > 0 and fac > 0:
        if mac == 2 and fac == 2 and ac == 1:
            return False
        return True

    if mac > 0:
        return False if ac == 2 else True

    if fac > 0:
        return False if ac == 2 else True

    return None


def classify_trio_genotypes(
    ac: int,
    mac: int,
    fac: int,
    m_gt: str,
    f_gt: str,
    c_gt: str,
    c_gq: str,
) -> tuple[str, str, str, str] | None:
    """Return the genotype tuple to store, or None if the site should be skipped."""
    if ac <= 0:
        return None

    if mac > 0 and fac > 0:
        if mac == 2 and fac == 2 and ac < 2:
            return m_gt, f_gt, c_gt, c_gq
        return m_gt, f_gt, c_gt, c_gq

    if mac > 0:
        if ac == 2:
            return m_gt, "0/0", c_gt, c_gq
        return m_gt, "0/0", c_gt, c_gq

    if fac > 0:
        if ac == 2:
            return "0/0", f_gt, c_gt, c_gq
        return "0/0", f_gt, c_gt, c_gq

    return None


def trio_bucket(
    ac: int,
    mac: int,
    fac: int,
) -> str | None:
    """Return ``inherited``, ``mendelian_bad``, or None."""
    result = classify_trio(ac, mac, fac)
    if result is None:
        return None
    return "inherited" if result else "mendelian_bad"
