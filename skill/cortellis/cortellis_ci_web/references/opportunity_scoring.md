# Opportunity scoring & interpretation (template)

This dashboard uses **heuristic, dataset-derived** scores to map “crowded” vs “blue-ocean” spaces across the selected grouping (e.g., target family, target type, target, modality, route, MoA class).

## Core rates (observed maturity)

Rates are computed **within the currently filtered slice** of assets:

- **Clinical entry rate** = fraction of assets at or beyond Phase I.
- **Late-stage rate** = fraction of *clinical* assets that are Phase III / Registered / Marketed.
- **Marketed rate** = fraction of assets that are Marketed.

These are **not** probability-of-success (PoS) estimates; they are snapshot maturity signals.

## Crowdedness and indices

For each group value `g`:

- **Crowdedness** is based on `log1p(N_g)` normalized to `[0,1]`.
- **Maturity composite** is a blend of late-stage + marketed signals.

Two indices are computed on a `0–100` scale:

- **Competition index**: higher when **crowded + mature** (more incumbents, more de-risked).
- **Blue-ocean index**: higher when **less crowded + mature** (fewer competitors but meaningful late/market signals).

## Role labels (quadrants)

Using medians (within the current view) for `N` and the selected y-metric:

- **Crowded + mature**: incumbent-heavy, de-risked; differentiation needs clear advantage (mechanistic, delivery, safety, access).
- **Crowded + early**: busy pipeline, limited maturity; whitespace often in specificity, biomarkers, differentiated delivery/route.
- **Blue-ocean + mature**: “quiet winners”; investigate adoption constraints (safety, delivery, patient selection) and why sparse.
- **Blue-ocean + early**: emerging, high uncertainty; focus on target validation and translational endpoints.

## How to use it

- Use the **Opportunity Lens bubble plot** to toggle global filters quickly.
- Use the **Opportunity tables** to compare groups and click rows for deeper breakdowns (stage mix, top routes/MoA/modality/targets).
