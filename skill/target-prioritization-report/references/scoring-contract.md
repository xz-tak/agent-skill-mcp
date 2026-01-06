# Scoring Contract (Abstract)

This is the documented contract your loaders + scorers should implement so the report remains traceable and consistent across indications.

## Key principles

- Preserve raw values in `data_dict` (`*_raw`) with a `source_file`.
- Normalize to 0–100 only at scoring time (or for display).
- Missing sources degrade gracefully to `50.0` (neutral), and should be visible in traceability output.

## Overall score

Default weighted sum (CI is folded into Opportunity; do not add CI as a standalone overall component):

`Overall = 0.30·Clinical + 0.30·Disease + 0.10·Safety + 0.20·Opportunity + 0.10·Novelty`

## Clinical Validation (0–100)

- Input: raw clinical score from a clinical report table/section.
- Normalize: `Clinical = clamp((raw / MAX_RAW)·100, 0, 100)`
- Combo: average member-gene raw, then normalize.

## Disease Association (0–100)

Components (example set; keep consistent in code + report):

- DEG evidence (**0.40**)
- BioBridge percentile (**0.25**)
- ULTRA percentile (**0.25**)
- PrimeKG connectivity (path-based score) (**0.10**)

Recommended disease score:

`Disease = 0.40·DEG + 0.25·BioBridge + 0.25·ULTRA + 0.10·PrimeKG`

Typical combo rule:
- Average member genes for DEG/ULTRA/PrimeKG after normalization.
- Use combo-level BioBridge percentile when available; otherwise fallback to a member gene.

### PrimeKG connectivity normalization

Preferred raw inputs from PrimeKG (from `results/kgpred_<indication>/primekg/...`):
- `primekg_connectivity_score_raw`: expected in `[0, 1]` where **higher = more connected** (shorter/stronger paths)
- optionally `primekg_path_length_raw`: integer path length (1 = direct)

Transform to 0–100 connectivity for scoring:

- If `primekg_connectivity_score_raw` is present:
  - `Connectivity = clamp(primekg_connectivity_score_raw * 100, 0, 100)`
- Else if only `primekg_path_length_raw` is present:
  - `score = 0.9^(path_length - 1)` (with `path_length >= 1`)
  - `Connectivity = clamp(score * 100, 0, 100)`
- Else:
  - `Connectivity = 50.0` (neutral)

This replaces older “connection-count / max-in-run” normalizations.

## Safety (0–100)

Compute a severity-weighted average from OFF‑X breakdown counts:

- Base scores: `very_high=0`, `high=10`, `medium=20`, `low=60`, `very_low=80`, `not_assoc=100`
- NA rows are excluded from the denominator (unknown).

`Safety = Σ(count[class]·score[class]) / Σ(count[class])`

Combo rule:
- Prefer combo-level OFF‑X breakdown if you have it; otherwise average gene safety scores.

## Competitive Intelligence (CI) (0–100)

Source: CI dashboard HTML with `<script id="data">` JSON.

Interpretation: higher = more whitespace / less competition.

Phase-weighted CI definition (gene-level):

- `Weighted = 1.0·Marketed + 0.7·PhaseIII + 0.4·PhaseII + 0.2·PhaseI + 0.1·Preclinical`
- `CI = 100 - 100·(Weighted / TotalPrograms)` (clamped 0–100)
- If CI missing: default `50.0`

Family-level CI definition (crowding):

- Each entry has a `targetFamilyPrimary` label (e.g., JAK, WNT, IL-23).
- For each entry: `w(entry) = weight(entry.ibdPhase)` using the same phase weights above.
- `FamilyWeighted = Σ w(entry)` for entries in that family (IBD-tagged entries only).
- `AllWeighted = Σ w(entry)` across all families (IBD-tagged entries only).
- `CI_family = 100 - 100·(FamilyWeighted / AllWeighted)` (clamped 0–100)

Blended CI (recommended default):

- `CI_blend = 0.5·CI_gene + 0.5·CI_family`

Combo rule (simple + stable):
- `CI_combo = mean(CI_blend)` across member genes with available CI.

## Opportunity (0–100)

Inputs:
- Disease (0–100)
- Clinical Novelty = `100 - Clinical`
- CI (0–100)

Individual:

`Opp = 0.40·Disease + 0.30·ClinicalNovelty + 0.30·CI_blend`

Combo (example pattern used in GI2):
- `MeanOpp = 0.40·Disease + 0.30·ClinicalNovelty + 0.30·CI_blend`
- `NovelMechanismBonus = 100 - 100·Synergy`
- `Opp_combo = 0.60·MeanOpp + 0.40·NovelMechanismBonus`

## Novelty (0–100 display)

Definition:
- Clinical Novelty = `100 - Clinical`
- Literature Novelty proxy (PrimeKG): **higher connectivity ⇒ lower novelty**

Compute:
- `LiteratureNovelty = 100 - Connectivity`
- `Novelty = 0.70·ClinicalNovelty + 0.30·LiteratureNovelty`

Combo rule:
- Compute `Connectivity` per gene, take the mean across genes, then compute `LiteratureNovelty`.

If your internal scoring adds a combo premium, the report may optionally display the base 0–100 equivalent for comparability.
