# Commercial Deep-Research Execution Protocol

When `/target-prioritization-report` is invoked, execute commercial scoring on-the-fly BEFORE report generation. This protocol runs within the active Claude session using Agent() calls — no external subprocess.

## References

- Prompt template: `references/commercial-research-prompt-template.md`
- Scoring contract: `references/commercial-scoring-contract.md`
- Preparation script: `scripts/run_commercial_scoring.py`

## Step 1: Load prompts

Read `<WORKDIR>/test/commercial/commercial_prompts.json`.

If the file doesn't exist, generate it:

```bash
python ~/.claude/skills/target-prioritization-report/scripts/run_commercial_scoring.py \
  --workdir <WORKDIR> --indication <IND> --mode prepare
```

Then copy from `results/commercial_<indication>/` to `test/commercial/`.

## Step 2: Identify existing results (skip-existing)

For each combo in the prompts JSON, check if a valid result already exists:

```
<WORKDIR>/test/commercial/<combo_slug>/commercial_result.json
```

Where `<combo_slug>` = genes joined by `_` (e.g., `TYK2_JAK1`, `TNFRSF25_GREM1`, `CDKN2D_ITGA4_ITGB7`).

**Decision logic:**
- If the file EXISTS and contains a valid `commercial_score` field (number 0-100) → **SKIP** this combo
- If the file is MISSING or does not contain a valid score → **ADD** to the fan-out list

Report to user: `"Found N existing results, M combos need research."`

## Step 3: Fan-out deep-research (3 runs per combo, parallel)

For EACH combo in the fan-out list, spawn **3 independent foreground Agent() calls** (3 runs × N combos total agents). Send ALL Agent tool calls in a **single message** so they run in parallel (concurrency cap ~16 handles throttling).

**Why 3 runs:** Each run independently researches and scores the combo. Final result = mean of scores + union of non-overlapping evidence. This reduces single-agent bias and captures broader evidence.

**Concrete Agent tool call pattern** (emit 3 per combo, vary the description suffix):

```
Agent({
  description: "Commercial: <COMBO_NAME> (run <1|2|3>)",
  prompt: "<combo prompt from commercial_prompts.json>\n\n
IMPORTANT — YOU MUST DO ALL OF THE FOLLOWING:

1. Use WebSearch to find: current market size for this indication, competing
   drugs/pipeline for these targets, strategic landscape in GI/immunology.
2. Use WebFetch to pull specific pages (ClinicalTrials.gov, company pipelines,
   analyst reports) that inform scoring.
3. Score each dimension 0-100 based on BOTH the local CI context above AND
   your web research findings.
4. Compute: commercial_score = 0.40*market_opportunity + 0.40*competitive_profile + 0.20*strategic_fit
5. Run this Bash command: mkdir -p <WORKDIR>/test/commercial/<combo_slug>
6. Use the Write tool to write your JSON result to:
   <WORKDIR>/test/commercial/<combo_slug>/commercial_result_run<1|2|3>.json
7. Return the JSON you wrote as your final message.

Output JSON schema (REQUIRED):
{
  \"combo\": \"GENE1 + GENE2\",
  \"indication\": \"<IND>\",
  \"commercial_score\": <float 0-100>,
  \"confidence\": \"<High|Medium|Low>\",
  \"components\": {
    \"market_opportunity\": {\"score\": <0-100>, \"rationale\": \"...\", \"key_evidence\": [...]},
    \"competitive_profile\": {\"score\": <0-100>, \"rationale\": \"...\", \"key_evidence\": [...]},
    \"strategic_fit\": {\"score\": <0-100>, \"rationale\": \"...\", \"key_evidence\": [...]}
  },
  \"citations\": [{\"source\": \"...\", \"relevance\": \"...\"}],
  \"narrative_summary\": \"<one paragraph>\"
}"
})
```

**Key rules:**
- Each agent handles exactly ONE combo, ONE run (simple prompt = better performance)
- 3 agents per combo with different run IDs (run1, run2, run3)
- Each writes to `commercial_result_run1.json`, `commercial_result_run2.json`, `commercial_result_run3.json`
- The combo's full `prompt` field from `commercial_prompts.json` is inserted at the top
- The agent MUST use WebSearch + WebFetch for real market data (not just local context)
- All agents are foreground (NOT `run_in_background`) — results return directly

## Step 3.1: Merge per-combo runs into final result

After all 3 runs complete for a combo, merge them into `commercial_result.json`:

```python
# Merge logic (executed by parent session):
runs = [run1, run2, run3]  # loaded from commercial_result_run{1,2,3}.json

final = {
  "combo": runs[0]["combo"],
  "indication": runs[0]["indication"],
  "commercial_score": mean([r["commercial_score"] for r in runs]),
  "confidence": majority_vote([r["confidence"] for r in runs]),
  "components": {
    dim: {
      "score": mean([r["components"][dim]["score"] for r in runs]),
      "rationale": runs[0]["components"][dim]["rationale"],  # from first run
      "key_evidence": union_dedup([e for r in runs for e in r["components"][dim]["key_evidence"]])
    }
    for dim in ["market_opportunity", "competitive_profile", "strategic_fit"]
  },
  "citations": union_dedup([c for r in runs for c in r["citations"]], key="source"),
  "narrative_summary": runs[0]["narrative_summary"],  # from highest-scoring run
  "runs": len(runs),
  "score_std": std([r["commercial_score"] for r in runs]),
}
```

**Union dedup rule:** For `key_evidence` and `citations`, collect all items across runs and deduplicate by exact string match (evidence) or by `source` field (citations). This maximizes evidence coverage while avoiding repetition.

**Score:** Mean of all 3 runs' `commercial_score` (and per-component scores).

**Narrative:** Use the narrative from the run with the highest commercial_score (most confident assessment).

Write merged result to: `<WORKDIR>/test/commercial/<combo_slug>/commercial_result.json`

**Root cause fixed:** Background agents go idle without returning results. Foreground agents (3 per combo, all in one message) return directly and write files as double-safety. Concurrency cap (~16) handles throttling.

## Step 3b: Fan-out per-gene commercial research (3 runs per gene, parallel)

In ADDITION to combo scoring, spawn **3 independent foreground Agent() calls per UNIQUE GENE** in the target set. Same 3-run pattern as combos.

**Skip-existing:** Check `<WORKDIR>/test/commercial/genes/<GENE>/commercial_result.json` — skip if exists and valid.

**Agent call pattern** (3 per gene, all in one message alongside combo agents):

```
Agent({
  description: "Commercial Gene: <GENE> (run <1|2|3>)",
  prompt: "Research the commercial landscape for <GENE> as an individual drug target in IBD.

## Assessment Criteria (same 3 dimensions, same formula)
Score the SINGLE TARGET (not a combination):
- Market Opportunity (40%): IBD market relevance for this mechanism
- Competitive Profile (40%): Existing/pipeline drugs targeting THIS gene specifically
- Strategic Fit (20%): GI/immunology alignment

## Local Pipeline Context
<GENE>: <N> active drugs (from CI dashboard)
[list programs if any]

## Instructions
1. Use WebSearch + WebFetch for real competitive data
2. Score each dimension 0-100
3. Compute: commercial_score = 0.40*market + 0.40*competitive + 0.20*strategic
4. Run: mkdir -p <WORKDIR>/test/commercial/genes/<GENE>
5. Write JSON to: <WORKDIR>/test/commercial/genes/<GENE>/commercial_result_run<1|2|3>.json
6. Return the JSON as your final message.

Output schema: identical to combo schema but 'combo' field = '<GENE>'"
})
```

**Output:** `test/commercial/genes/<GENE>/commercial_result_run{1,2,3}.json`

**Merge (Step 3.1 same logic):** After all 3 gene runs complete, merge into `test/commercial/genes/<GENE>/commercial_result.json` using the same mean-scores + union-evidence rule as combos.

The loader (`load_commercial_gene`) checks per-gene merged result FIRST, falls back to mean-of-combos.

## Step 4: Verify per-combo AND per-gene result files

After all foreground agents return, verify each combo has its result file:

```bash
ls <WORKDIR>/test/commercial/*/commercial_result.json | wc -l
```

Compare count against expected (N combos in fan-out list).

For any combo where the file is missing (agent failed or returned without writing):
- Re-send that single Agent() call with the same prompt (foreground, blocking)
- If it fails again, write the result manually from the agent's returned JSON text

All per-combo results sit alongside each other in `test/commercial/<combo_slug>/`.

## Step 5: Combine all results + update summary

Read ALL per-combo `commercial_result.json` files from `test/commercial/*/`.

Merge existing + newly generated results into:

1. **`<WORKDIR>/test/commercial/commercial_scores.json`** — array of all combo score objects
2. **`<WORKDIR>/test/commercial/commercial_summary.md`** — human-readable summary:

```markdown
# Commercial Scoring Summary

Generated: <timestamp>

## Score Table

| Combo | Score | Confidence | Market | Competitive | Strategic |
|-------|-------|------------|--------|-------------|-----------|
| TYK2 + JAK1 | 45.2 | High | 72 | 28 | 85 |
| TNFRSF25 + GREM1 | 78.6 | Medium | 70 | 88 | 82 |
| ... | ... | ... | ... | ... | ... |

## Per-Combo Narratives

### TYK2 + JAK1
<narrative_summary from result>

### TNFRSF25 + GREM1
<narrative_summary from result>

...
```

## Step 6: Validate

```bash
python ~/.claude/skills/target-prioritization-report/scripts/run_commercial_scoring.py \
  --workdir <WORKDIR> --indication <IND> --mode validate \
  --scores-file <WORKDIR>/test/commercial/commercial_scores.json
```

## Step 7: Proceed to report generation

```bash
python tools/generate_ibd_report.py --output reports/<IND>_Combo_Prioritization_Report_offline.html
```

`load_all_data()` reads `test/commercial/commercial_scores.json` → commercial scores flow into the Overall formula (0.10 weight) + combo cards in the HTML report.

---

## Expected Directory Structure

```
test/commercial/
  commercial_prompts.json          # input prompts (all combos)
  commercial_scores.json           # combined scores (all combos)
  commercial_summary.md            # human-readable summary table + narratives
  TYK2_JAK1/
    commercial_result.json         # individual combo result
  TNFRSF25_GREM1/
    commercial_result.json
  TNFRSF25_PCOLCE/
    commercial_result.json
  CDKN2D_ITGA4_ITGB7/
    commercial_result.json
  ... (one subdir per combo)
```

## Incremental Execution

Because of skip-existing logic, this protocol supports incremental runs:
- First invocation: researches all 25 combos, saves results
- Subsequent invocations: finds existing results, skips them, only researches new combos (if any added to config)
- Force re-research: delete the relevant `<combo_slug>/commercial_result.json` file(s) before invoking
