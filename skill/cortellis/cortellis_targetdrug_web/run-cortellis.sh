#!/bin/bash
# Cortellis search automation with export functionality
# This script will work from any directory

# Save the original working directory (where user called the script from)
ORIGINAL_DIR="$(pwd)"

# Get the directory where the script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the Playwright skill directory (adjust if needed for your environment)
PLAYWRIGHT_SKILL_DIR="${HOME}/.claude/plugins/cache/playwright-skill/playwright-skill/4.1.0/skills/playwright-skill"

# Check if playwright skill exists
if [ ! -d "$PLAYWRIGHT_SKILL_DIR" ]; then
    echo "Error: Playwright skill directory not found at: $PLAYWRIGHT_SKILL_DIR"
    echo "Please update the PLAYWRIGHT_SKILL_DIR variable in this script to point to your playwright-skill directory"
    exit 1
fi

# Export the original directory so the JS script can use it
export CORTELLIS_WORK_DIR="$ORIGINAL_DIR"

# Run the automation from playwright skill directory
cd "$PLAYWRIGHT_SKILL_DIR"
node run.js "$SCRIPT_DIR/cortellis-automation.js" "$@"
