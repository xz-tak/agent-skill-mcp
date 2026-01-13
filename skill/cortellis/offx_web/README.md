# OFF-X Web Automation - Quick Reference

Browser automation for OFF-X (TargetSafety.info) to evaluate target and drug safety and download Master view data.

## Purpose

Evaluate target and drug safety and identify adverse events using OFF-X database via web automation.

## Quick Start

```bash
# 1. One-time authentication setup (uses centralized okta-sso skill)
~/ai-sci-claude-skills/ai-sci/okta-sso/run-okta-login.sh

# 2. Run safety assessment
/home/sagemaker-user/.claude/skills/cortellis/offx_web/run-offx-download.sh itga4
```

## Common Commands

### Single Entity Assessment

```bash
# Default field (Targets), all matches
./run-offx-download.sh itga4

# Specific field
./run-offx-download.sh imatinib --field "Drugs and biologics"

# Exact match only
./run-offx-download.sh egfr --exact
```

### Multiple Entities

```bash
# Multiple targets
./run-offx-download.sh --queries "tyk2,jak1,itga4"

# Multiple drugs
./run-offx-download.sh --queries "imatinib,gefitinib" --field "Drugs and biologics"
```

### Config File

```bash
./run-offx-download.sh --config example/example_targets.json
```

## Available Fields

- **Targets** (default) - Target safety profiles
- **Drugs and biologics** - Drug safety data
- **Drug combinations** - Combination therapy safety
- **Adverse events** - Adverse event analysis

## Output

Results saved to `offx_playwright_result/` in working directory:

```
offx_playwright_result/
└── {entity}_{field}_{datetime}/
    ├── *.xlsx (Master view Excel files)
    ├── *.png (Screenshots)
    └── metadata.json
```

## Example Configs

Located in `example/` directory:
- `example_targets.json` - Target safety
- `example_drugs_biologics.json` - Drug safety
- `example_adverse_events.json` - Adverse events
- `example_comprehensive.json` - Mixed analysis

## Troubleshooting

### Authentication Expired
```bash
# Check status and re-authenticate
~/ai-sci-claude-skills/ai-sci/okta-sso/run-okta-login.sh --status
~/ai-sci-claude-skills/ai-sci/okta-sso/run-okta-login.sh
```

### Debug Issues
Check screenshots in output directory:
- `dropdown_matches.png` - Search results
- `{match}_page.png` - Entity page
- `{match}_master_view.png` - Master view table

## Documentation

**Detailed Guide:** [`OFFX_AUTOMATION.md`](OFFX_AUTOMATION.md)

**Parent Skill:** [`../SKILL.md`](../SKILL.md)
