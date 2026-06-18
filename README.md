# Inherited

Analyze rare variants in family trios from a VCF, precomputed gnomAD allele-frequency JSON, and a family relations file.

## Run

**Without installing** (clone and run from the project directory):

```bash
python run.py analyze \
  --vcf chr22.vcf.gz \
  --af-json gnomad_chr22.json \
  --family-file families.tsv \
  -o results/chr22
```

**Or install** so the `inherited` command is on your PATH:

```bash
pip install -e ".[dev]"
inherited analyze --vcf chr22.vcf.gz --af-json gnomad_chr22.json --family-file families.tsv -o results/chr22
```

## Family file format

Tab-separated file (required via `--family-file`):

```text
spid    family_id    father_id    mother_id
child1  fam1         fa1          ma1
fa1     fam1         0            0
ma1     fam1         0            0
```

Complete trios (`father_id` and `mother_id` both not `0`) are used for analysis.

## gnomAD AF JSON

```json
{
  "var_key": 0.0001,
  "22:12345:A:G": {"AF": 0.001, "AF_EUR": 0.0005}
}
```

Variants with AF above the threshold are skipped. Missing keys default to `0` and are kept.

## Usage

```bash
python run.py analyze \
  --vcf chr22.vcf.gz \
  --af-json gnomad_chr22.json \
  --family-file families.tsv \
  -o results/chr22
```

Disable multiallelic per-ALT processing:

```bash
python run.py analyze --vcf chr22.vcf.gz --af-json gnomad.json --family-file families.tsv -o results/chr22 --no-multiallelic
```

## Output

The output directory contains:

- `inherited.json` — rare inherited variant calls per variant key and child sample
- `mendelian_bad.json` — mendelian-inconsistent calls
- `inherited_per_variant.json` — variant key → number of people with that inherited variant
- `inherited_per_person.json` — person id → number of inherited variants for that person
- `mendelian_bad_per_gt.json` — `mother_gt:father_gt:child_gt` → number of people with that pattern
- `stats.json` — summary counts
- `params.json` — parameters used for this run

Each chromosome output directory (e.g. `results/chr22/`) contains its own `params.json` alongside the result files.

Each entry in `inherited.json` / `mendelian_bad.json` stores `(mother_GT, father_GT, child_GT, child_GQ)`.

## Tests

```bash
pytest
```
