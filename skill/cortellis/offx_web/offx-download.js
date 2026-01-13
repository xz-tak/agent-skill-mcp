const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');
const os = require('os');

// Configuration
const TARGET_URL = 'https://takeda.okta.com/';
// Use environment variable if set (from shell script), otherwise use current working directory
const WORK_DIR = process.env.OFFX_WORK_DIR || process.cwd();
// Centralized Okta auth state (managed by okta-sso skill)
const AUTH_STATE_PATH = path.join(os.homedir(), '.okta', 'auth_state.json');
const OUTPUT_DIR = path.join(WORK_DIR, 'offx_playwright_result');
const DEFAULT_FIELD = 'Targets';

// Validate auth state exists before proceeding
function validateAuthState() {
  if (!fs.existsSync(AUTH_STATE_PATH)) {
    console.error('');
    console.error('❌ Okta session not found at ' + AUTH_STATE_PATH);
    console.error('');
    console.error('To authenticate, run:');
    console.error('  ~/ai-sci-claude-skills/ai-sci/okta-sso/run-okta-login.sh');
    console.error('');
    console.error('Then retry this command.');
    console.error('');
    process.exit(1);
  }
}

// Available fields in off-x
const AVAILABLE_FIELDS = [
  'Targets',
  'Drugs and biologics',
  'Drug combinations',
  'Adverse events'
];

// Parse command line arguments
function parseArgs() {
  let args = process.argv.slice(2);
  args = args.filter(arg => !arg.endsWith('.js'));

  if (args.length === 0) {
    console.log('Usage:');
    console.log('  Option 1: Use config file');
    console.log('    ./run-offx-download.sh --config <config.json>');
    console.log('');
    console.log('  Option 2: Direct arguments (single query)');
    console.log('    ./run-offx-download.sh <entity_name> [--field "Field Name"] [--exact]');
    console.log('');
    console.log('  Option 3: Multiple queries with command line');
    console.log('    ./run-offx-download.sh --queries "entity1,entity2" [--field "Field Name"] [--exact]');
    console.log('');
    console.log('Config file format (JSON):');
    console.log('[');
    console.log('  {');
    console.log('    "entity": "tyk2",');
    console.log('    "field": "Targets",');
    console.log('    "exactMatch": false  // Optional, defaults to false (iterate all matched)');
    console.log('  },');
    console.log('  {');
    console.log('    "entity": "jak1",');
    console.log('    "field": "Drugs and biologics"');
    console.log('  }');
    console.log(']');
    console.log('');
    console.log('Examples:');
    console.log('  # Config file');
    console.log('  ./run-offx-download.sh --config my_searches.json');
    console.log('');
    console.log('  # Single query, default field (Targets), all matches');
    console.log('  ./run-offx-download.sh tyk2');
    console.log('');
    console.log('  # Single query, specific field, all matches');
    console.log('  ./run-offx-download.sh jak1 --field "Drugs and biologics"');
    console.log('');
    console.log('  # Single query, exact match only');
    console.log('  ./run-offx-download.sh grem1 --exact');
    console.log('');
    console.log('  # Multiple queries, default field');
    console.log('  ./run-offx-download.sh --queries "tyk2,jak1,grem1"');
    console.log('');
    console.log('  # Multiple queries, specific field');
    console.log('  ./run-offx-download.sh --queries "tyk2,jak1" --field "Targets"');
    console.log('');
    console.log('Available fields:');
    AVAILABLE_FIELDS.forEach(field => console.log(`  - ${field}`));
    process.exit(1);
  }

  // Check for --config
  if (args[0] === '--config' || args[0] === '-c') {
    if (args.length < 2) {
      console.error('Error: --config option requires a file path');
      process.exit(1);
    }
    const configPath = args[1];
    if (!fs.existsSync(configPath)) {
      console.error(`Error: Config file not found: ${configPath}`);
      process.exit(1);
    }
    try {
      const config = JSON.parse(fs.readFileSync(configPath, 'utf-8'));
      if (!Array.isArray(config)) {
        console.error('Error: Config must be an array');
        process.exit(1);
      }
      // Apply defaults to config entries with field-specific exactMatch defaults
      return config.map(entry => {
        const field = entry.field || DEFAULT_FIELD;
        // Default exactMatch based on field:
        // - Targets and Drugs and biologics: false (get all matches)
        // - Drug combinations and Adverse events: true (get only first match)
        let defaultExactMatch = false;
        if (field === 'Drug combinations' || field === 'Adverse events') {
          defaultExactMatch = true;
        }

        return {
          entity: entry.entity,
          field: field,
          exactMatch: entry.exactMatch !== undefined ? entry.exactMatch : defaultExactMatch
        };
      });
    } catch (e) {
      console.error('Error parsing config file:', e.message);
      process.exit(1);
    }
  }

  // Parse command line arguments
  let entities = [];
  let field = DEFAULT_FIELD;
  let exactMatch = false;

  // Check for --queries
  const queriesIndex = args.findIndex(arg => arg === '--queries' || arg === '-q');
  if (queriesIndex !== -1) {
    if (queriesIndex + 1 >= args.length) {
      console.error('Error: --queries option requires a comma-separated list');
      process.exit(1);
    }
    entities = args[queriesIndex + 1].split(',').map(s => s.trim()).filter(s => s);
    args.splice(queriesIndex, 2);
  } else if (args.length > 0 && !args[0].startsWith('--')) {
    // First argument is a single entity
    entities = [args[0]];
    args.splice(0, 1);
  }

  // Check for --field
  const fieldIndex = args.findIndex(arg => arg === '--field' || arg === '-f');
  if (fieldIndex !== -1) {
    if (fieldIndex + 1 >= args.length) {
      console.error('Error: --field option requires a field name');
      process.exit(1);
    }
    field = args[fieldIndex + 1];
    if (!AVAILABLE_FIELDS.includes(field)) {
      console.error(`Error: Invalid field "${field}". Available fields: ${AVAILABLE_FIELDS.join(', ')}`);
      process.exit(1);
    }
  }

  // Check for --exact (can override defaults)
  let exactMatchOverride = null;
  if (args.includes('--exact') || args.includes('-e')) {
    exactMatchOverride = true;
  }

  if (entities.length === 0) {
    console.error('Error: No entities provided');
    process.exit(1);
  }

  // Default exactMatch based on field:
  // - Targets and Drugs and biologics: false (get all matches)
  // - Drug combinations and Adverse events: true (get only first match)
  let defaultExactMatch = false;
  if (field === 'Drug combinations' || field === 'Adverse events') {
    defaultExactMatch = true;
  }

  const finalExactMatch = exactMatchOverride !== null ? exactMatchOverride : defaultExactMatch;

  // Create config format: each entity gets the same field
  return entities.map(entity => ({
    entity,
    field,
    exactMatch: finalExactMatch
  }));
}

// Setup directories
function setupDirectories() {
  if (!fs.existsSync(OUTPUT_DIR)) {
    fs.mkdirSync(OUTPUT_DIR, { recursive: true });
  }
}

// Select field in off-x
async function selectField(page, fieldName) {
  console.log(`\n🎯 Selecting field: "${fieldName}"`);

  try {
    // Look for the field selector - could be tabs, buttons, or links
    const fieldSelectors = [
      `[role="tab"]:has-text("${fieldName}")`,
      `button:has-text("${fieldName}")`,
      `a:has-text("${fieldName}")`,
      `[role="button"]:has-text("${fieldName}")`,
      `label:has-text("${fieldName}")`,
      `div[class*="tab"]:has-text("${fieldName}")`,
      `div:has-text("${fieldName}")`
    ];

    for (const selector of fieldSelectors) {
      try {
        const elements = await page.locator(selector).all();
        for (const element of elements) {
          if (await element.isVisible({ timeout: 1000 })) {
            const text = await element.innerText().catch(() => '');
            // Make sure it's an exact match (or close enough)
            if (text.trim() === fieldName || text.includes(fieldName)) {
              console.log(`✅ Found field selector, clicking "${fieldName}"...`);
              await element.click();
              await page.waitForTimeout(1000);
              return true;
            }
          }
        }
      } catch (e) {
        // Try next selector
      }
    }

    // If not found as button/tab, might be in a dropdown
    console.log('🔍 Looking for field dropdown...');
    const dropdownSelectors = [
      'select',
      '[role="combobox"]',
      'button[aria-haspopup="listbox"]'
    ];

    for (const selector of dropdownSelectors) {
      try {
        const dropdown = await page.$(selector);
        if (dropdown && await dropdown.isVisible()) {
          await dropdown.click();
          await page.waitForTimeout(500);

          const optionSelectors = [
            `option:has-text("${fieldName}")`,
            `[role="option"]:has-text("${fieldName}")`,
            `li:has-text("${fieldName}")`
          ];

          for (const optionSelector of optionSelectors) {
            try {
              const option = await page.locator(optionSelector).first();
              if (await option.isVisible({ timeout: 1000 })) {
                await option.click();
                await page.waitForTimeout(1000);
                console.log(`✅ Selected field: "${fieldName}"`);
                return true;
              }
            } catch (e) {
              // Try next option selector
            }
          }
        }
      } catch (e) {
        // Try next dropdown selector
      }
    }

    console.log(`⚠️  Could not find field selector for "${fieldName}"`);
    return false;
  } catch (e) {
    console.error(`❌ Error selecting field:`, e.message);
    return false;
  }
}

// Get dropdown matches for entity search
async function getDropdownMatches(page, entity) {
  console.log(`\n🔍 Getting dropdown matches for: "${entity}"`);

  try {
    // Wait longer for dropdown to fully load with items
    await page.waitForTimeout(8000);

    // Wait specifically for the listbox with options to appear
    let dropdownMatches = [];
    let optionsFound = false;

    try {
      await page.waitForSelector('[role="listbox"] [role="option"]', { timeout: 5000 });
      console.log('✅ Dropdown with options detected');

      // Directly get all [role="option"] elements from the page
      // This works regardless of the listbox visibility
      const items = await page.$$('[role="option"]');

      if (items.length > 0) {
        console.log(`✅ Found dropdown: [role="listbox"]`);
        console.log(`📋 Found ${items.length} items using selector: [role="option"]`);

        for (const item of items) {
          const text = await item.innerText().catch(() => '');
          const trimmedText = text.trim();

          // Skip empty items
          if (!trimmedText) continue;

          // Ignore "select multiple targets" button
          if (trimmedText.toLowerCase().includes('select multiple')) {
            console.log(`⏭️  Skipping: "${trimmedText}" (select multiple button)`);
            continue;
          }

          // Ignore anything from/below "Machine learning / Pathway maps results" label
          if (
            trimmedText.toLowerCase().includes('machine learning') ||
            trimmedText.toLowerCase().includes('pathway maps')
          ) {
            console.log(`⏭️  Reached ML/Pathway section, stopping: "${trimmedText}"`);
            break;
          }

          // Add valid match
          dropdownMatches.push({
            text: trimmedText,
            element: item
          });
        }

        optionsFound = true;
      }
    } catch (e) {
      console.log('⏳ Dropdown not detected with role selectors, trying other selectors...');
    }

    // If options were found, filter for exact matches and return
    if (optionsFound && dropdownMatches.length > 0) {
      console.log(`📋 Found ${dropdownMatches.length} dropdown items before filtering`);

      // Filter: only keep matches where text before "[" exactly matches entity (ignoring spaces, hyphens, underscores)
      const entityNormalized = entity.replace(/[\s\-_]+/g, '').toLowerCase();
      const exactMatches = dropdownMatches.filter(match => {
        // Extract text before "["
        const textBeforeBracket = match.text.split('[')[0].trim();
        const textNormalized = textBeforeBracket.replace(/[\s\-_]+/g, '').toLowerCase();

        const isExactMatch = textNormalized === entityNormalized;
        if (!isExactMatch) {
          console.log(`⏭️  Skipping: "${match.text}" (before bracket: "${textBeforeBracket}" doesn't match "${entity}")`);
        }
        return isExactMatch;
      });

      if (exactMatches.length > 0) {
        console.log(`✅ Found ${exactMatches.length} exact match(es):`);
        exactMatches.forEach((match, i) => {
          console.log(`   ${i + 1}. ${match.text}`);
        });
        return exactMatches;
      } else {
        console.log(`⚠️  No exact matches found for "${entity}"`);
        return [];
      }
    }

    // Fallback: try other selectors if listbox didn't work
    console.log('⏳ Trying fallback dropdown selectors...');
    const dropdownSelectors = [
      '[role="menu"]',
      '.autocomplete-results',
      '.dropdown-menu',
      '.search-results',
      '[class*="dropdown"]',
      '[class*="autocomplete"]'
    ];

    for (const selector of dropdownSelectors) {
      try {
        const dropdown = await page.$(selector);
        if (dropdown && await dropdown.isVisible({ timeout: 2000 })) {
          console.log(`✅ Found dropdown: ${selector}`);

          // Get all options/items in the dropdown
          const itemSelectors = [
            '[role="option"]',
            '[role="menuitem"]',
            'li',
            '.result-item',
            '[class*="option"]',
            '[class*="item"]'
          ];

          for (const itemSelector of itemSelectors) {
            try {
              const items = await dropdown.$$(itemSelector);
              if (items.length > 0) {
                console.log(`📋 Found ${items.length} items using selector: ${itemSelector}`);

                for (const item of items) {
                  const text = await item.innerText().catch(() => '');
                  const trimmedText = text.trim();

                  // Skip empty items
                  if (!trimmedText) continue;

                  // Ignore "select multiple targets" button
                  if (trimmedText.toLowerCase().includes('select multiple')) {
                    console.log(`⏭️  Skipping: "${trimmedText}" (select multiple button)`);
                    continue;
                  }

                  // Ignore anything from/below "Machine learning / Pathway maps results" label
                  if (
                    trimmedText.toLowerCase().includes('machine learning') ||
                    trimmedText.toLowerCase().includes('pathway maps')
                  ) {
                    console.log(`⏭️  Reached ML/Pathway section, stopping: "${trimmedText}"`);
                    break;
                  }

                  // Add valid match
                  dropdownMatches.push({
                    text: trimmedText,
                    element: item
                  });
                }

                if (dropdownMatches.length > 0) {
                  break; // Found items, stop searching
                }
              } else {
                console.log(`   No items found with selector: ${itemSelector}`);
              }
            } catch (e) {
              // Try next item selector
            }
          }

          if (dropdownMatches.length > 0) {
            break; // Found dropdown with items, stop searching
          } else {
            console.log(`   Dropdown ${selector} has no valid items, trying next selector...`);
          }
        }
      } catch (e) {
        console.log(`   Selector ${selector} not found or not visible`);
      }
    }

    if (dropdownMatches.length === 0) {
      console.log('⚠️  No dropdown matches found');
      return [];
    }

    // Filter fallback matches for exact match as well
    console.log(`📋 Found ${dropdownMatches.length} dropdown items before filtering (fallback)`);
    const entityNormalized = entity.replace(/[\s\-_]+/g, '').toLowerCase();
    const exactMatches = dropdownMatches.filter(match => {
      const textBeforeBracket = match.text.split('[')[0].trim();
      const textNormalized = textBeforeBracket.replace(/[\s\-_]+/g, '').toLowerCase();
      const isExactMatch = textNormalized === entityNormalized;
      if (!isExactMatch) {
        console.log(`⏭️  Skipping: "${match.text}" (before bracket: "${textBeforeBracket}" doesn't match "${entity}")`);
      }
      return isExactMatch;
    });

    if (exactMatches.length > 0) {
      console.log(`✅ Found ${exactMatches.length} exact match(es) (fallback):`);
      exactMatches.forEach((match, i) => {
        console.log(`   ${i + 1}. ${match.text}`);
      });
      return exactMatches;
    } else {
      console.log(`⚠️  No exact matches found for "${entity}" (fallback)`);
      return [];
    }
  } catch (e) {
    console.error(`❌ Error getting dropdown matches:`, e.message);
    return [];
  }
}

// Helper function to export a single safety profile
async function exportSingleProfile(page, matchText, profileType, resultDir, filenameSuffix = '') {
  console.log(`\n📊 Exporting ${profileType} safety profile for: "${matchText}"`);

  try {
    const matchSlug = matchText.replace(/[^a-z0-9]/gi, '_').toLowerCase();

    // Step 1: Click the appropriate safety profile based on profileType
    let profileName = '';
    let profileSelectors = [];

    if (profileType === 'target') {
      profileName = 'Target safety profile';
      profileSelectors = [
        'button:has-text("Target safety profile")',
        'a:has-text("Target safety profile")',
        '[role="button"]:has-text("Target safety profile")',
        'div:has-text("Target safety profile")',
        '*:has-text("Target safety profile")'
      ];
    } else if (profileType === 'drug') {
      profileName = 'Drug safety profile';
      profileSelectors = [
        'button:has-text("Drug safety profile")',
        'a:has-text("Drug safety profile")',
        '[role="tab"]:has-text("Drug safety profile")',
        '*:has-text("Drug safety profile")'
      ];
    } else if (profileType === 'combination') {
      profileName = 'Combination safety profile';
      profileSelectors = [
        'button:has-text("Combination safety profile")',
        'a:has-text("Combination safety profile")',
        '[role="tab"]:has-text("Combination safety profile")',
        '*:has-text("Combination safety profile")'
      ];
    } else {
      return { success: false, reason: `Unknown profile type: ${profileType}` };
    }

    console.log(`🔍 Looking for "${profileName}" button...`);
    let profileClicked = false;
    for (const selector of profileSelectors) {
      try {
        const element = await page.locator(selector).first();
        if (await element.isVisible({ timeout: 2000 })) {
          console.log(`✅ Found "${profileName}", clicking...`);
          await element.click();
          await page.waitForTimeout(3000);
          profileClicked = true;
          break;
        }
      } catch (e) {
        // Try next selector
      }
    }

    if (!profileClicked) {
      console.log(`⚠️  Could not find "${profileName}" button`);
      return { success: false, reason: `${profileName} button not found` };
    }

    // Take screenshot after clicking safety profile
    const afterSafetyScreenshot = path.join(resultDir, `${matchSlug}${filenameSuffix}_after_safety_profile.png`);
    await page.screenshot({ path: afterSafetyScreenshot, fullPage: true });
    console.log(`📸 After safety profile: ${afterSafetyScreenshot}`);

    // Step 2: Select "Master view" from dropdown
    console.log('🔍 Looking for "Master view" in dropdown below safety profile...');

    // Wait for dropdown to appear below the safety profile button
    await page.waitForTimeout(2000);

    // Try more specific selectors for Master view option
    const masterViewSelectors = [
      '[role="option"]:has-text("Master view")',
      '[role="menuitem"]:has-text("Master view")',
      'li:has-text("Master view")',
      'button:has-text("Master view")',
      'a:has-text("Master view")',
      'div[role="button"]:has-text("Master view")',
      '*:has-text("Master view")'
    ];

    let masterViewClicked = false;
    for (const selector of masterViewSelectors) {
      try {
        console.log(`   Trying selector: ${selector}`);
        const elements = await page.locator(selector).all();
        console.log(`   Found ${elements.length} elements`);

        for (const element of elements) {
          try {
            if (await element.isVisible({ timeout: 1000 })) {
              const text = await element.innerText().catch(() => '');
              console.log(`   Visible element text: "${text}"`);

              // Check if this is exactly "Master view" (case insensitive)
              if (text.toLowerCase().trim() === 'master view') {
                console.log('✅ Found "Master view", clicking...');
                await element.click();
                await page.waitForTimeout(3000);
                masterViewClicked = true;
                break;
              }
            }
          } catch (e) {
            // Try next element
          }
        }

        if (masterViewClicked) break;
      } catch (e) {
        // Try next selector
      }
    }

    if (!masterViewClicked) {
      console.log('⚠️  Could not find "Master view" option');
      return { success: false, reason: 'Master view option not found' };
    }

    // Take screenshot after selecting master view
    const afterMasterViewScreenshot = path.join(resultDir, `${matchSlug}${filenameSuffix}_master_view.png`);
    await page.screenshot({ path: afterMasterViewScreenshot, fullPage: true });
    console.log(`📸 Master view: ${afterMasterViewScreenshot}`);

    // Step 3: Click Export icon (far right of table)
    console.log('🔍 Looking for Export icon...');
    const exportSelectors = [
      'button[aria-label="Export" i]',
      'button[title="Export" i]',
      '[aria-label*="Export" i]',
      '[title*="Export" i]',
      'button:has-text("Export")',
      'svg[aria-label="Export" i]',
      // Look for icons on the far right
      'button svg',
      '[class*="export"]',
      '[class*="Export"]'
    ];

    let exportInitiated = false;

    for (const selector of exportSelectors) {
      try {
        const elements = await page.locator(selector).all();

        for (const element of elements) {
          try {
            if (await element.isVisible({ timeout: 500 })) {
              const ariaLabel = await element.getAttribute('aria-label').catch(() => null);
              const title = await element.getAttribute('title').catch(() => null);

              const buttonInfo = ariaLabel || title || selector;

              if (buttonInfo.toLowerCase().includes('export')) {
                console.log(`✅ Found Export button: "${buttonInfo}"`);

                // Set up download listener before clicking
                const downloadPromise = page.waitForEvent('download', { timeout: 30000 }).catch(() => null);

                await element.click();
                await page.waitForTimeout(2000);

                // Check if download started
                const download = await downloadPromise;
                if (download) {
                  const downloadFileName = download.suggestedFilename() || `${matchSlug}${filenameSuffix}_export.xlsx`;
                  const downloadPath = path.join(resultDir, downloadFileName);
                  await download.saveAs(downloadPath);
                  console.log(`✅ Master view exported: ${downloadPath}`);

                  const stats = fs.statSync(downloadPath);
                  const fileSizeInMB = (stats.size / (1024 * 1024)).toFixed(2);
                  console.log(`   File size: ${fileSizeInMB} MB`);

                  exportInitiated = true;
                  return {
                    success: true,
                    filePath: downloadPath,
                    fileName: downloadFileName,
                    fileSize: stats.size,
                    fileSizeMB: fileSizeInMB
                  };
                }
              }
            }
          } catch (e) {
            // Continue to next element
          }
        }

        if (exportInitiated) break;
      } catch (e) {
        // Continue to next selector
      }
    }

    if (!exportInitiated) {
      console.log('⚠️  Could not initiate export - button may not be available');
      return { success: false, reason: 'Export button not found' };
    }

  } catch (e) {
    console.error(`❌ Error during ${profileType} profile export:`, e.message);
    return { success: false, reason: e.message };
  }
}

// Helper function to determine profile type based on field
function getProfileType(field) {
  if (field === 'Targets') {
    return 'target';
  } else if (field === 'Drugs and biologics') {
    return 'drug';
  } else if (field === 'Drug combinations') {
    return 'combination';
  } else if (field === 'Adverse events') {
    return 'both'; // Special case: export both target and drug reports
  }
  return 'unknown';
}

// Download Master view table export
async function downloadMasterView(page, entity, matchText, field, resultDir) {
  console.log(`\n💾 Downloading Master view for: "${matchText}" (Field: ${field})`);

  try {
    // Wait for page to load
    await page.waitForTimeout(3000);

    // Handle any cookie consent popup that might appear
    try {
      const acceptButton = await page.$('#onetrust-accept-btn-handler');
      if (acceptButton && await acceptButton.isVisible({ timeout: 1000 })) {
        await acceptButton.click({ force: true });
        await page.waitForTimeout(1000);
      }
    } catch (e) {
      // No cookie banner
    }

    // Take screenshot of the page
    const matchSlug = matchText.replace(/[^a-z0-9]/gi, '_').toLowerCase();
    const pageScreenshot = path.join(resultDir, `${matchSlug}_page.png`);
    await page.screenshot({ path: pageScreenshot, fullPage: true });
    console.log(`📸 Page screenshot: ${pageScreenshot}`);

    // Field-specific logic
    const profileType = getProfileType(field);

    // Special handling for Adverse events: export both target and drug profiles
    if (profileType === 'both') {
      console.log('ℹ️  Adverse events field detected - will export both Target and Drug safety profiles');

      const results = [];

      // Export 1: Target safety profile
      console.log('\n📊 Export 1/2: Target safety profile');
      const targetResult = await exportSingleProfile(page, matchText, 'target', resultDir, '_target');
      results.push({ profileType: 'target', ...targetResult });

      // Export 2: Drug safety profile
      console.log('\n📊 Export 2/2: Drug safety profile');
      const drugResult = await exportSingleProfile(page, matchText, 'drug', resultDir, '_drug');
      results.push({ profileType: 'drug', ...drugResult });

      // Return combined results
      const allSuccess = results.every(r => r.success);
      if (allSuccess) {
        return {
          success: true,
          multiProfile: true,
          results: results,
          filePaths: results.map(r => r.filePath).filter(f => f),
          fileNames: results.map(r => r.fileName).filter(f => f)
        };
      } else {
        return {
          success: false,
          multiProfile: true,
          results: results,
          reason: 'One or more profiles failed to export'
        };
      }
    }

    if (profileType === 'target') {
      // Step 1: Click "Target safety profile" on the left side
      console.log('🔍 Looking for "Target safety profile" button...');
      const safetyProfileSelectors = [
        'button:has-text("Target safety profile")',
        'a:has-text("Target safety profile")',
        '[role="button"]:has-text("Target safety profile")',
        'div:has-text("Target safety profile")',
        '*:has-text("Target safety profile")'
      ];

      let safetyProfileClicked = false;
      for (const selector of safetyProfileSelectors) {
        try {
          const element = await page.locator(selector).first();
          if (await element.isVisible({ timeout: 2000 })) {
            console.log('✅ Found "Target safety profile", clicking...');
            await element.click();
            await page.waitForTimeout(2000);
            safetyProfileClicked = true;
            break;
          }
        } catch (e) {
          // Try next selector
        }
      }

      if (!safetyProfileClicked) {
        console.log('⚠️  Could not find "Target safety profile" button');
        return { success: false, reason: 'Target safety profile button not found' };
      }
    } else if (profileType === 'drug') {
      // Step 1: Click "Drug safety profile" tab
      console.log('🔍 Looking for "Drug safety profile" tab...');
      const drugSafetySelectors = [
        'button:has-text("Drug safety profile")',
        'a:has-text("Drug safety profile")',
        '[role="tab"]:has-text("Drug safety profile")',
        '*:has-text("Drug safety profile")'
      ];

      let drugSafetyClicked = false;
      for (const selector of drugSafetySelectors) {
        try {
          const element = await page.locator(selector).first();
          if (await element.isVisible({ timeout: 2000 })) {
            console.log('✅ Found "Drug safety profile", clicking...');
            await element.click();
            await page.waitForTimeout(2000);
            drugSafetyClicked = true;
            break;
          }
        } catch (e) {
          // Try next selector
        }
      }

      if (!drugSafetyClicked) {
        console.log('⚠️  Could not find "Drug safety profile" tab');
        return { success: false, reason: 'Drug safety profile tab not found' };
      }
    } else if (profileType === 'combination') {
      // Step 1: Click "Combination safety profile" tab
      console.log('🔍 Looking for "Combination safety profile" tab...');
      const comboSafetySelectors = [
        'button:has-text("Combination safety profile")',
        'a:has-text("Combination safety profile")',
        '[role="tab"]:has-text("Combination safety profile")',
        '*:has-text("Combination safety profile")'
      ];

      let comboSafetyClicked = false;
      for (const selector of comboSafetySelectors) {
        try {
          const element = await page.locator(selector).first();
          if (await element.isVisible({ timeout: 2000 })) {
            console.log('✅ Found "Combination safety profile", clicking...');
            await element.click();
            await page.waitForTimeout(2000);
            comboSafetyClicked = true;
            break;
          }
        } catch (e) {
          // Try next selector
        }
      }

      if (!comboSafetyClicked) {
        console.log('⚠️  Could not find "Combination safety profile" tab');
        return { success: false, reason: 'Combination safety profile tab not found' };
      }
    } else {
      console.log(`⚠️  Unknown profile type for field "${field}"`);
      return { success: false, reason: 'Unknown field type' };
    }

    // Take screenshot after clicking safety profile
    const afterSafetyScreenshot = path.join(resultDir, `${matchSlug}_after_safety_profile.png`);
    await page.screenshot({ path: afterSafetyScreenshot, fullPage: true });
    console.log(`📸 After safety profile: ${afterSafetyScreenshot}`);

    // Step 2: Select "Master view" from dropdown
    console.log('🔍 Looking for "Master view" in dropdown...');
    const masterViewSelectors = [
      'button:has-text("Master view")',
      'a:has-text("Master view")',
      '[role="option"]:has-text("Master view")',
      '[role="menuitem"]:has-text("Master view")',
      'li:has-text("Master view")',
      'div:has-text("Master view")'
    ];

    let masterViewClicked = false;
    for (const selector of masterViewSelectors) {
      try {
        const element = await page.locator(selector).first();
        if (await element.isVisible({ timeout: 2000 })) {
          console.log('✅ Found "Master view", clicking...');
          await element.click();
          await page.waitForTimeout(3000);
          masterViewClicked = true;
          break;
        }
      } catch (e) {
        // Try next selector
      }
    }

    if (!masterViewClicked) {
      console.log('⚠️  Could not find "Master view" option');
      return { success: false, reason: 'Master view option not found' };
    }

    // Take screenshot after selecting master view
    const afterMasterViewScreenshot = path.join(resultDir, `${matchSlug}_master_view.png`);
    await page.screenshot({ path: afterMasterViewScreenshot, fullPage: true });
    console.log(`📸 Master view: ${afterMasterViewScreenshot}`);

    // Step 3: Click Export icon (far right of table)
    console.log('🔍 Looking for Export icon...');
    const exportSelectors = [
      'button[aria-label="Export" i]',
      'button[title="Export" i]',
      '[aria-label*="Export" i]',
      '[title*="Export" i]',
      'button:has-text("Export")',
      'svg[aria-label="Export" i]',
      // Look for icons on the far right
      'button svg',
      '[class*="export"]',
      '[class*="Export"]'
    ];

    let exportInitiated = false;

    for (const selector of exportSelectors) {
      try {
        const elements = await page.locator(selector).all();

        for (const element of elements) {
          try {
            if (await element.isVisible({ timeout: 500 })) {
              const ariaLabel = await element.getAttribute('aria-label').catch(() => null);
              const title = await element.getAttribute('title').catch(() => null);

              const buttonInfo = ariaLabel || title || selector;

              if (buttonInfo.toLowerCase().includes('export')) {
                console.log(`✅ Found Export button: "${buttonInfo}"`);

                // Set up download listener before clicking
                const downloadPromise = page.waitForEvent('download', { timeout: 30000 }).catch(() => null);

                await element.click();
                await page.waitForTimeout(2000);

                // Check if download started
                const download = await downloadPromise;
                if (download) {
                  const downloadFileName = download.suggestedFilename() || `${matchSlug}_export.xlsx`;
                  const downloadPath = path.join(resultDir, downloadFileName);
                  await download.saveAs(downloadPath);
                  console.log(`✅ Master view exported: ${downloadPath}`);

                  const stats = fs.statSync(downloadPath);
                  const fileSizeInMB = (stats.size / (1024 * 1024)).toFixed(2);
                  console.log(`   File size: ${fileSizeInMB} MB`);

                  exportInitiated = true;
                  return {
                    success: true,
                    filePath: downloadPath,
                    fileName: downloadFileName,
                    fileSize: stats.size,
                    fileSizeMB: fileSizeInMB
                  };
                }
              }
            }
          } catch (e) {
            // Continue to next element
          }
        }

        if (exportInitiated) break;
      } catch (e) {
        // Continue to next selector
      }
    }

    if (!exportInitiated) {
      console.log('⚠️  Could not initiate export - button may not be available');
      return { success: false, reason: 'Export button not found' };
    }

  } catch (e) {
    console.error(`❌ Error during download:`, e.message);
    return { success: false, reason: e.message };
  }
}

// Search for next match using top right search box
async function searchNextMatch(page, entity, matchText) {
  console.log(`\n🔎 Searching for next match: "${matchText}"`);

  try {
    // Look for search box in top right area
    const searchSelectors = [
      'input[type="search"]',
      'input[placeholder*="Search" i]',
      'input[name="search"]',
      '[role="searchbox"]'
    ];

    let searchInput = null;
    for (const selector of searchSelectors) {
      const elements = await page.$$(selector);
      // Try to find the top-right search box (usually the last one or second one)
      if (elements.length > 0) {
        searchInput = elements[elements.length - 1]; // Try last search box (likely top right)
        if (await searchInput.isVisible()) {
          break;
        }
      }
    }

    if (!searchInput) {
      console.log('⚠️  Could not find top-right search box');
      return false;
    }

    console.log('✅ Found search box, entering entity...');
    await searchInput.click();
    await searchInput.fill('');
    await page.waitForTimeout(500);
    await searchInput.fill(entity);
    await page.waitForTimeout(5000);

    // Look for the specific match in dropdown
    console.log(`🔍 Looking for "${matchText}" in dropdown...`);
    const dropdownMatches = await getDropdownMatches(page, entity);

    for (const match of dropdownMatches) {
      if (match.text === matchText) {
        console.log(`✅ Found "${matchText}", clicking...`);
        await match.element.click();
        await page.waitForTimeout(3000);
        return true;
      }
    }

    console.log(`⚠️  Could not find "${matchText}" in dropdown`);
    return false;
  } catch (e) {
    console.error(`❌ Error searching for next match:`, e.message);
    return false;
  }
}

// Helper function to dismiss cookie consent overlays
async function dismissCookieOverlay(page) {
  try {
    // Try to close the OneTrust privacy center overlay if it exists
    const closeSelectors = [
      '.onetrust-close-btn-handler',
      '#close-pc-btn-handler',
      'button[aria-label="Close"]',
      '.ot-pc-refuse-all-handler',
      '#onetrust-pc-btn-handler'
    ];

    for (const selector of closeSelectors) {
      try {
        const closeBtn = await page.$(selector);
        if (closeBtn && await closeBtn.isVisible({ timeout: 500 })) {
          await closeBtn.click({ force: true });
          console.log(`✅ Dismissed cookie overlay using: ${selector}`);
          await page.waitForTimeout(1000);
          return true;
        }
      } catch (e) {
        // Continue to next selector
      }
    }

    // Try pressing ESC key to dismiss
    await page.keyboard.press('Escape');
    await page.waitForTimeout(500);

    // Check if the overlay is gone
    const overlay = await page.$('.onetrust-pc-dark-filter');
    if (!overlay || !await overlay.isVisible({ timeout: 500 })) {
      console.log('✅ Cookie overlay dismissed with ESC key');
      return true;
    }
  } catch (e) {
    // No overlay to dismiss
  }
  return false;
}

// Process one entity with its field
async function processEntity(offxPage, entityConfig, index, total) {
  const { entity, field, exactMatch } = entityConfig;

  console.log(`\n${'='.repeat(70)}`);
  console.log(`🔍 Processing ${index}/${total}: "${entity}"`);
  console.log(`   Field: ${field}`);
  console.log(`   Exact match only: ${exactMatch}`);
  console.log('='.repeat(70));

  const entitySlug = entity.replace(/[^a-z0-9]/gi, '_').toLowerCase();
  const fieldSlug = field.replace(/[^a-z0-9]/gi, '_').toLowerCase();
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19).replace('T', '_');
  const resultDir = path.join(OUTPUT_DIR, `${entitySlug}_${fieldSlug}_${timestamp}`);
  fs.mkdirSync(resultDir, { recursive: true });

  const results = {
    entity,
    field,
    exactMatch,
    downloads: [],
    timestamp: new Date().toISOString()
  };

  // Select field
  const fieldSelected = await selectField(offxPage, field);
  if (!fieldSelected) {
    console.log('⚠️  Could not explicitly select field, but it may already be selected');
    // Don't fail here - field might already be selected or tabs work differently
  }

  // Find search box (below the field selector, center area)
  console.log(`\n🔎 Searching for: "${entity}"`);

  // Use more specific selectors for the off-x search box
  const searchSelectors = [
    'input[placeholder*="Type target" i]',
    'input[placeholder*="gene" i]',
    'input[placeholder*="UniProt" i]',
    'input[type="text"][placeholder]',
    'input[type="search"]',
    'input[placeholder*="Search" i]'
  ];

  let searchInput = null;
  for (const selector of searchSelectors) {
    try {
      const element = await offxPage.$(selector);
      if (element && await element.isVisible()) {
        searchInput = element;
        console.log(`✅ Found search box with selector: ${selector}`);
        break;
      }
    } catch (e) {
      // Try next selector
    }
  }

  if (!searchInput) {
    console.log('❌ Could not find search box');
    return { success: false, ...results, reason: 'Search box not found' };
  }

  // Perform search
  await searchInput.click();
  await searchInput.fill('');
  await offxPage.waitForTimeout(500);
  await searchInput.fill(entity);
  await offxPage.waitForTimeout(5000);

  console.log('⏳ Waiting for dropdown results...');

  // Get dropdown matches
  const dropdownMatches = await getDropdownMatches(offxPage, entity);

  if (dropdownMatches.length === 0) {
    console.log('❌ No matches found in dropdown');
    return { success: false, ...results, reason: 'No matches found' };
  }

  // If exactMatch is true, find the match that begins with "entity [entity]" pattern
  let matchesToProcess;
  if (exactMatch) {
    // Look for pattern: "entity [entity]" (case-insensitive)
    const entityLower = entity.toLowerCase();
    const exactPattern = `${entityLower} [${entityLower}]`;

    const exactMatchFound = dropdownMatches.find(match =>
      match.text.toLowerCase().startsWith(exactPattern)
    );

    if (exactMatchFound) {
      console.log(`✅ Found exact match pattern: "${exactMatchFound.text}"`);
      matchesToProcess = [exactMatchFound];
    } else {
      console.log(`⚠️  No exact match pattern found for "${entity} [${entity}]", using first match`);
      matchesToProcess = [dropdownMatches[0]];
    }
  } else {
    matchesToProcess = dropdownMatches;
  }

  console.log(`\n📋 Will process ${matchesToProcess.length} match(es)`);

  // Take screenshot of dropdown
  const dropdownScreenshot = path.join(resultDir, 'dropdown_matches.png');
  await offxPage.screenshot({ path: dropdownScreenshot, fullPage: true });
  console.log(`📸 Dropdown screenshot: ${dropdownScreenshot}`);

  // Process first match by clicking in dropdown
  console.log(`\n--- Match 1/${matchesToProcess.length}: ${matchesToProcess[0].text} ---`);
  // Use force click to bypass any cookie overlays
  await matchesToProcess[0].element.click({ force: true });
  await offxPage.waitForTimeout(3000);

  console.log(`✅ Navigated to page for: "${matchesToProcess[0].text}"`);
  console.log(`   URL: ${offxPage.url()}`);

  // Download Master view for first match
  const downloadResult = await downloadMasterView(offxPage, entity, matchesToProcess[0].text, field, resultDir);
  results.downloads.push({
    match: matchesToProcess[0].text,
    ...downloadResult
  });

  // Process remaining matches by going back to home and searching again
  for (let i = 1; i < matchesToProcess.length; i++) {
    const match = matchesToProcess[i];
    console.log(`\n--- Match ${i + 1}/${matchesToProcess.length}: ${match.text} ---`);

    // Navigate back to home page
    console.log('🏠 Navigating back to home page...');
    await offxPage.goto('https://www.targetsafety.info/home', { waitUntil: 'networkidle', timeout: 30000 });
    await offxPage.waitForTimeout(2000);

    // Select field again
    await selectField(offxPage, field);

    // Search for entity again
    console.log(`🔎 Searching for "${entity}" again...`);
    const searchInput2 = await offxPage.$('input[placeholder*="Type target" i]');
    if (searchInput2) {
      await searchInput2.click();
      await searchInput2.fill('');
      await offxPage.waitForTimeout(500);
      await searchInput2.fill(entity);
      await offxPage.waitForTimeout(5000);
    }

    // Get dropdown matches again and click on the specific one
    const dropdownMatches2 = await getDropdownMatches(offxPage, entity);
    let matchFound = false;
    let targetMatch = null;
    for (const dropdownMatch of dropdownMatches2) {
      if (dropdownMatch.text === match.text) {
        targetMatch = dropdownMatch;
        matchFound = true;
        break;
      }
    }

    if (matchFound && targetMatch) {
      console.log(`✅ Found "${match.text}" in dropdown, clicking...`);
      // Use force click to bypass any cookie overlays
      await targetMatch.element.click({ force: true });
      await offxPage.waitForTimeout(3000);
    }

    if (!matchFound) {
      console.log(`⚠️  Could not find "${match.text}" in dropdown`);
      results.downloads.push({
        match: match.text,
        success: false,
        reason: 'Could not find match in dropdown'
      });
      continue;
    }

    console.log(`✅ Navigated to page for: "${match.text}"`);
    console.log(`   URL: ${offxPage.url()}`);

    // Download Master view for this match
    const downloadResult = await downloadMasterView(offxPage, entity, match.text, field, resultDir);
    results.downloads.push({
      match: match.text,
      ...downloadResult
    });
  }

  // Save results metadata
  const metadataPath = path.join(resultDir, 'metadata.json');
  fs.writeFileSync(metadataPath, JSON.stringify(results, null, 2));
  console.log(`\n💾 Metadata saved: ${metadataPath}`);

  // Print summary
  console.log(`\n✅ Completed processing: "${entity}"`);
  console.log(`   Successful downloads: ${results.downloads.filter(d => d.success).length}/${matchesToProcess.length}`);

  return { success: true, ...results };
}

// Main execution
(async () => {
  // Validate Okta session before starting
  validateAuthState();

  const entityConfigs = parseArgs();
  console.log(`\n🎯 Will process ${entityConfigs.length} entit(ies):`);
  entityConfigs.forEach((config, i) => {
    console.log(`   ${i + 1}. ${config.entity} (${config.field}, exactMatch: ${config.exactMatch})`);
  });

  setupDirectories();

  const browser = await chromium.launch({
    headless: true, // Set to true for headless mode
    slowMo: 500,
    args: [
      '--disable-blink-features=AutomationControlled',
      '--disable-dev-shm-usage',
      '--no-sandbox'
    ]
  });

  try {
    console.log('\n🔐 Loading authentication...');
    const context = await browser.newContext({
      storageState: AUTH_STATE_PATH,
      acceptDownloads: true,
      userAgent: 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
      viewport: { width: 1920, height: 1080 },
      locale: 'en-US',
      timezoneId: 'America/New_York',
      extraHTTPHeaders: {
        'Accept-Language': 'en-US,en;q=0.9',
        'Accept-Encoding': 'gzip, deflate, br',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8'
      }
    });

    // Add stealth scripts
    await context.addInitScript(() => {
      Object.defineProperty(navigator, 'webdriver', { get: () => false });
      window.chrome = { runtime: {} };
      const originalQuery = window.navigator.permissions.query;
      window.navigator.permissions.query = (parameters) => (
        parameters.name === 'notifications' ?
          Promise.resolve({ state: Notification.permission }) :
          originalQuery(parameters)
      );
      Object.defineProperty(navigator, 'plugins', { get: () => [1, 2, 3, 4, 5] });
      Object.defineProperty(navigator, 'languages', { get: () => ['en-US', 'en'] });
    });

    const page = await context.newPage();

    // Navigate to Okta and launch off-x
    console.log('\n🌐 Navigating to Okta portal...');
    await page.goto(TARGET_URL, { waitUntil: 'networkidle', timeout: 60000 });
    await page.waitForTimeout(3000);

    // Find app search in top left
    console.log('🔍 Looking for app search box...');
    const appSearchInput = await page.$('input[type="search"]');
    if (!appSearchInput) throw new Error('Could not find app search');

    // Search for off-x
    console.log('🔎 Searching for "off-x"...');
    await appSearchInput.click();
    await appSearchInput.fill('off-x');
    await page.waitForTimeout(2000);

    // Launch off-x
    console.log('🚀 Launching off-x...');
    const offxLinkSelectors = [
      'a:has-text("off-x")',
      'a:has-text("Off-X")',
      'a:has-text("OFF-X")',
      '[class*="app"]:has-text("off-x")'
    ];

    let offxLink = null;
    for (const selector of offxLinkSelectors) {
      try {
        const link = await page.locator(selector).first();
        if (await link.isVisible({ timeout: 2000 })) {
          offxLink = link;
          break;
        }
      } catch (e) {
        // Try next selector
      }
    }

    if (!offxLink) {
      throw new Error('Could not find off-x application link');
    }

    const newPagePromise = context.waitForEvent('page', { timeout: 15000 }).catch(() => null);
    await offxLink.click();

    let offxPage = await newPagePromise;
    if (!offxPage) {
      await page.waitForNavigation({ timeout: 10000 }).catch(() => {});
      offxPage = page;
    }

    console.log('⏳ Loading off-x...');
    await offxPage.waitForTimeout(5000);

    // Handle Cloudflare if present
    try {
      const pageTitle = await offxPage.title();
      if (pageTitle.includes('Just a moment') || pageTitle.includes('Verify')) {
        console.log('🔐 Cloudflare detected, waiting...');
        for (let i = 0; i < 15; i++) {
          await offxPage.waitForTimeout(3000);
          const currentTitle = await offxPage.title().catch(() => '');
          if (!currentTitle.includes('Just a moment') && !currentTitle.includes('Verify')) {
            break;
          }
        }
      }
    } catch (e) {
      await offxPage.waitForTimeout(2000);
    }

    await offxPage.waitForTimeout(3000);
    console.log('✅ off-x loaded');
    console.log(`   URL: ${offxPage.url()}`);

    // Handle cookie consent
    try {
      const acceptButton = await offxPage.$('#onetrust-accept-btn-handler');
      if (acceptButton && await acceptButton.isVisible({ timeout: 2000 })) {
        await acceptButton.click({ force: true });
        await offxPage.waitForTimeout(2000);
      }
    } catch (e) {
      // No cookie banner
    }

    // Dismiss any cookie privacy center overlay
    await dismissCookieOverlay(offxPage);

    // Process each entity
    const allResults = [];
    for (let i = 0; i < entityConfigs.length; i++) {
      const result = await processEntity(offxPage, entityConfigs[i], i + 1, entityConfigs.length);
      allResults.push(result);

      if (i < entityConfigs.length - 1) {
        console.log('\n⏳ Waiting 3 seconds before next search...');
        await offxPage.waitForTimeout(3000);
      }
    }

    // Generate final summary
    console.log('\n' + '='.repeat(70));
    console.log('📊 FINAL SUMMARY');
    console.log('='.repeat(70));

    const summary = {
      timestamp: new Date().toISOString(),
      totalEntities: entityConfigs.length,
      totalDownloadsSuccessful: allResults.reduce((sum, r) => sum + (r.downloads ? r.downloads.filter(d => d.success).length : 0), 0),
      results: allResults
    };

    const timestamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19).replace('T', '_');
    const summaryPath = path.join(OUTPUT_DIR, `summary_${timestamp}.json`);
    fs.writeFileSync(summaryPath, JSON.stringify(summary, null, 2));

    console.log(`\n💾 Summary: ${summaryPath}`);
    console.log(`📁 Results: ${OUTPUT_DIR}`);
    console.log(`\n✅ Completed ${summary.totalDownloadsSuccessful} downloads`);

    allResults.forEach((r, i) => {
      const successful = r.downloads ? r.downloads.filter(d => d.success).length : 0;
      const total = r.downloads ? r.downloads.length : 0;
      console.log(`   ${i + 1}. ${r.entity}: ${successful}/${total} downloads`);
    });

  } catch (error) {
    console.error('\n❌ Error:', error.message);
    console.error(error.stack);
    process.exit(1);
  } finally {
    await browser.close();
    console.log('\n✅ Browser closed');
  }
})();
