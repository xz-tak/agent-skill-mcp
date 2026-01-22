const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');

// Configuration
const TARGET_URL = 'https://takeda.okta.com/';
// Use environment variable if set (from shell script), otherwise use current working directory
const WORK_DIR = process.env.CORTELLIS_WORK_DIR || process.cwd();
const AUTH_STATE_PATH = path.join(WORK_DIR, 'okta_auth_state.json');
const OUTPUT_DIR = path.join(WORK_DIR, 'cortellis_playwright_result');
const DEFAULT_CATEGORY = 'Drugs & Biologics';

// Parse command line arguments
function parseArgs() {
  let args = process.argv.slice(2);
  args = args.filter(arg => !arg.endsWith('.js'));

  if (args.length === 0) {
    console.log('Usage:');
    console.log('  Option 1: Use config file');
    console.log('    ./run-cortellis-download.sh --config <config.json>');
    console.log('');
    console.log('  Option 2: Direct arguments (single query)');
    console.log('    ./run-cortellis-download.sh <search_term> [--categories "Cat1,Cat2"]');
    console.log('');
    console.log('  Option 3: Multiple queries with command line');
    console.log('    ./run-cortellis-download.sh --queries "drug1,drug2" [--categories "Cat1,Cat2"]');
    console.log('');
    console.log('Config file format (JSON):');
    console.log('[');
    console.log('  {');
    console.log('    "searchTerm": "imatinib",');
    console.log('    "categories": ["Drugs & Biologics", "Clinical Studies"]');
    console.log('  },');
    console.log('  {');
    console.log('    "searchTerm": "gefitinib",');
    console.log('    "categories": ["Genes & Targets"]  // Optional, defaults to "Drugs & Biologics"');
    console.log('  }');
    console.log(']');
    console.log('');
    console.log('Examples:');
    console.log('  # Config file');
    console.log('  ./run-cortellis-download.sh --config my_searches.json');
    console.log('');
    console.log('  # Single query, default category');
    console.log('  ./run-cortellis-download.sh imatinib');
    console.log('');
    console.log('  # Single query, specific categories');
    console.log('  ./run-cortellis-download.sh imatinib --categories "Clinical Studies,Patents"');
    console.log('');
    console.log('  # Multiple queries, default category');
    console.log('  ./run-cortellis-download.sh --queries "imatinib,gefitinib,dasatinib"');
    console.log('');
    console.log('  # Multiple queries, specific categories');
    console.log('  ./run-cortellis-download.sh --queries "imatinib,gefitinib" --categories "Clinical Studies,Patents"');
    console.log('');
    console.log('Available categories:');
    console.log('  - Drugs & Biologics (default)');
    console.log('  - Genes & Targets');
    console.log('  - Organic Synthesis');
    console.log('  - Experimental Pharmacology');
    console.log('  - Experimental Models');
    console.log('  - Pharmacokinetics');
    console.log('  - Drug Metabolism');
    console.log('  - Drug-Drug Interactions');
    console.log('  - Clinical Studies');
    console.log('  - Organizations');
    console.log('  - Literature');
    console.log('  - Patents');
    console.log('  - Disease Briefings');
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
      // Apply defaults to config entries
      return config.map(entry => ({
        searchTerm: entry.searchTerm,
        categories: entry.categories && entry.categories.length > 0 ? entry.categories : [DEFAULT_CATEGORY]
      }));
    } catch (e) {
      console.error('Error parsing config file:', e.message);
      process.exit(1);
    }
  }

  // Parse command line arguments
  let searchTerms = [];
  let categories = [DEFAULT_CATEGORY];

  // Check for --queries
  const queriesIndex = args.findIndex(arg => arg === '--queries' || arg === '-q');
  if (queriesIndex !== -1) {
    if (queriesIndex + 1 >= args.length) {
      console.error('Error: --queries option requires a comma-separated list');
      process.exit(1);
    }
    searchTerms = args[queriesIndex + 1].split(',').map(s => s.trim()).filter(s => s);
    args.splice(queriesIndex, 2);
  } else if (args.length > 0 && !args[0].startsWith('--')) {
    // First argument is a single search term
    searchTerms = [args[0]];
    args.splice(0, 1);
  }

  // Check for --categories
  const categoriesIndex = args.findIndex(arg => arg === '--categories' || arg === '-c');
  if (categoriesIndex !== -1) {
    if (categoriesIndex + 1 >= args.length) {
      console.error('Error: --categories option requires a comma-separated list');
      process.exit(1);
    }
    categories = args[categoriesIndex + 1].split(',').map(s => s.trim()).filter(s => s);
  }

  if (searchTerms.length === 0) {
    console.error('Error: No search terms provided');
    process.exit(1);
  }

  // Create config format: each search term gets all categories
  return searchTerms.map(searchTerm => ({
    searchTerm,
    categories
  }));
}

// Setup directories
function setupDirectories() {
  if (!fs.existsSync(OUTPUT_DIR)) {
    fs.mkdirSync(OUTPUT_DIR, { recursive: true });
  }
}

// Click "View results" for a specific category
async function clickViewResults(page, categoryName) {
  console.log(`\n📊 Looking for category: "${categoryName}"`);

  // Wait for results to be visible
  await page.waitForTimeout(2000);

  // Try to find the category card and its "View results" button
  try {
    // Find all result cards
    const cards = await page.$$('[class*="card"], [class*="result"], [class*="tile"]');

    for (const card of cards) {
      const text = await card.innerText().catch(() => '');

      // Check if this card matches the category
      if (text.includes(categoryName)) {
        console.log(`✅ Found category card: "${categoryName}"`);

        // Look for "View results" button within or near this card
        const viewButtonSelectors = [
          'button:has-text("View results")',
          'a:has-text("View results")',
          '[aria-label*="View results"]',
          'button:has-text("View")',
          'a:has-text("View")'
        ];

        for (const selector of viewButtonSelectors) {
          try {
            const button = await card.$(selector);
            if (button) {
              console.log(`🖱️  Clicking "View results" for ${categoryName}...`);
              await button.click();
              await page.waitForTimeout(3000);
              return true;
            }
          } catch (e) {
            // Try next selector
          }
        }

        // If no button found in card, try clicking the card itself
        console.log('⚠️  No "View results" button found, trying to click card...');
        await card.click();
        await page.waitForTimeout(3000);
        return true;
      }
    }

    console.log(`⚠️  Could not find category: "${categoryName}"`);
    return false;
  } catch (e) {
    console.error(`❌ Error finding category ${categoryName}:`, e.message);
    return false;
  }
}

// Click "..." menu and export (same as cortellis-automation.js)
async function downloadResults(page, searchTerm, categoryName, resultDir) {
  console.log(`\n💾 Looking for export menu...`);

  try {
    // Wait for page to load
    await page.waitForTimeout(2000);

    // Take screenshot of the results page
    const categorySlug = categoryName.replace(/[^a-z0-9]/gi, '_');
    const screenshotPath = path.join(resultDir, `${categorySlug}_page.png`);
    await page.screenshot({ path: screenshotPath, fullPage: true });
    console.log(`📸 Results page screenshot: ${screenshotPath}`);

    console.log('🔍 Scanning page for export/menu buttons...');

    // Look for "..." menu button or export triggers
    const possibleExportTriggers = [
      // Icon buttons
      'button svg',
      'button[class*="icon"]',
      '[class*="IconButton"]',
      '[class*="icon-button"]',
      // Specific aria labels
      'button[aria-label*="export" i]',
      'button[aria-label*="download" i]',
      'button[aria-label*="more" i]',
      'button[aria-label*="menu" i]',
      'button[aria-label*="options" i]',
      '[aria-label*="export" i]',
      '[aria-label*="download" i]',
      // Title attributes
      'button[title*="export" i]',
      'button[title*="download" i]',
      'button[title*="more" i]',
      'button[title*="menu" i]',
      // Text-based
      'button:has-text("Export")',
      'button:has-text("Download")',
      'button:has-text("...")',
      'a:has-text("Export")',
      'a:has-text("Download")',
      // Common toolbar/action buttons
      '[role="toolbar"] button',
      '[class*="toolbar"] button',
      '[class*="actions"] button',
      '[class*="action-bar"] button',
      // Data test IDs
      '[data-testid*="export"]',
      '[data-testid*="download"]',
      '[data-testid*="menu"]',
      '[data-testid*="more"]'
    ];

    let exportInitiated = false;

    // Try each selector
    for (const selector of possibleExportTriggers) {
      try {
        const elements = await page.locator(selector).all();

        for (const element of elements) {
          try {
            if (await element.isVisible({ timeout: 500 })) {
              // Try to get button text or aria-label to understand what it is
              const ariaLabel = await element.getAttribute('aria-label').catch(() => null);
              const title = await element.getAttribute('title').catch(() => null);
              const text = await element.innerText().catch(() => '');

              const buttonInfo = ariaLabel || title || text || selector;

              // Check if this looks like an export button
              if (
                buttonInfo.toLowerCase().includes('export') ||
                buttonInfo.toLowerCase().includes('download') ||
                buttonInfo.toLowerCase().includes('...') ||
                (buttonInfo.toLowerCase().includes('more') && !buttonInfo.toLowerCase().includes('learn more'))
              ) {
                console.log(`✅ Found potential export trigger: "${buttonInfo}" (${selector})`);

                // Set up download listener before clicking
                const downloadPromise = page.waitForEvent('download', { timeout: 30000 }).catch(() => null);

                await element.click();
                await page.waitForTimeout(2000);

                // Check if download started immediately
                const download = await downloadPromise;
                if (download) {
                  const downloadFileName = download.suggestedFilename() || `export_${searchTerm.replace(/[^a-z0-9]/gi, '_')}.xlsx`;
                  const downloadPath = path.join(resultDir, `${categorySlug}_${downloadFileName}`);
                  await download.saveAs(downloadPath);
                  console.log(`✅ Export downloaded directly: ${downloadPath}`);

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

                // Take screenshot after clicking to see what opened
                const afterClickPath = path.join(resultDir, `${categorySlug}_after_button_click.png`);
                await page.screenshot({ path: afterClickPath, fullPage: true });
                console.log(`📸 Screenshot after button click: ${afterClickPath}`);

                // Now look for export option in any opened menu/dialog
                console.log('💾 Looking for "Export" option in menu...');
                const exportMenuSelectors = [
                  'button:has-text("Export")',
                  'a:has-text("Export")',
                  '[role="menuitem"]:has-text("Export")',
                  'li:has-text("Export")',
                  'div[role="menuitem"]:has-text("Export")',
                  '[class*="menu"] button:has-text("Export")',
                  '[class*="dropdown"] button:has-text("Export")',
                  '[class*="menu-item"]:has-text("Export")'
                ];

                for (const menuSelector of exportMenuSelectors) {
                  try {
                    const menuElement = page.locator(menuSelector).first();
                    if (await menuElement.isVisible({ timeout: 2000 })) {
                      console.log('✅ Found "Export" in menu, clicking...');

                      const downloadPromise2 = page.waitForEvent('download', { timeout: 30000 }).catch(() => null);

                      await menuElement.click();
                      await page.waitForTimeout(2000);

                      // Check if download started
                      const download2 = await downloadPromise2;
                      if (download2) {
                        const downloadFileName = download2.suggestedFilename() || `export_${searchTerm.replace(/[^a-z0-9]/gi, '_')}.xlsx`;
                        const downloadPath = path.join(resultDir, `${categorySlug}_${downloadFileName}`);
                        await download2.saveAs(downloadPath);
                        console.log(`✅ Export downloaded: ${downloadPath}`);

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

                      // If no download yet, might need to confirm in a dialog
                      console.log('⬇️ Looking for confirmation button...');
                      const confirmSelectors = [
                        'button:has-text("Export")',
                        'button:has-text("Download")',
                        'button:has-text("Confirm")',
                        'button:has-text("OK")',
                        'button[type="submit"]',
                        '[role="dialog"] button:has-text("Export")',
                        '.modal button:has-text("Export")'
                      ];

                      for (const confirmSelector of confirmSelectors) {
                        try {
                          const confirmElements = await page.locator(confirmSelector).all();
                          for (const confirmElement of confirmElements) {
                            if (await confirmElement.isVisible({ timeout: 1000 })) {
                              console.log(`✅ Found confirmation button, clicking...`);

                              const downloadPromise3 = page.waitForEvent('download', { timeout: 30000 }).catch(() => null);

                              await confirmElement.click();
                              await page.waitForTimeout(3000);

                              const download3 = await downloadPromise3;
                              if (download3) {
                                const downloadFileName = download3.suggestedFilename() || `export_${searchTerm.replace(/[^a-z0-9]/gi, '_')}.xlsx`;
                                const downloadPath = path.join(resultDir, `${categorySlug}_${downloadFileName}`);
                                await download3.saveAs(downloadPath);
                                console.log(`✅ Export downloaded: ${downloadPath}`);

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
                          if (exportInitiated) break;
                        } catch (e) {
                          // Continue
                        }
                      }

                      if (exportInitiated) break;
                    }
                  } catch (e) {
                    // Continue
                  }
                }

                if (exportInitiated) break;
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
      console.log('⚠️  Could not initiate export - button may not be available or export may require different steps');
      console.log('💡 Check the screenshots to see the actual UI');
      return { success: false, reason: 'Export button/menu not found' };
    }

  } catch (e) {
    console.error(`❌ Error during export:`, e.message);
    return { success: false, reason: e.message };
  }
}

// Process one search term with its categories
async function processSearch(cortellisPage, searchConfig, index, total) {
  const { searchTerm, categories } = searchConfig;

  console.log(`\n${'='.repeat(70)}`);
  console.log(`🔍 Processing ${index}/${total}: "${searchTerm}"`);
  console.log(`   Categories to download: ${categories.join(', ')}`);
  console.log('='.repeat(70));

  const searchTermSlug = searchTerm.replace(/[^a-z0-9]/gi, '_').toLowerCase();
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19).replace('T', '_');
  const categoriesSlug = categories.map(c => c.replace(/[^a-z0-9]/gi, '_')).join('-');
  const resultDir = path.join(OUTPUT_DIR, `${searchTermSlug}_${categoriesSlug}_${timestamp}`);
  fs.mkdirSync(resultDir, { recursive: true });

  const results = {
    searchTerm,
    requestedCategories: categories,
    downloads: [],
    timestamp: new Date().toISOString()
  };

  // Find and use search input
  console.log(`\n🔎 Searching for: "${searchTerm}"`);

  const searchSelectors = [
    'input[type="search"]',
    'input[placeholder*="Search" i]',
    'input[name="search"]',
    '[role="searchbox"]'
  ];

  let searchInput = null;
  for (const selector of searchSelectors) {
    const element = await cortellisPage.$(selector);
    if (element && await element.isVisible()) {
      searchInput = element;
      break;
    }
  }

  if (!searchInput) {
    console.log('❌ Could not find search box');
    return { success: false, ...results };
  }

  // Perform search
  await searchInput.click();
  await searchInput.fill('');
  await cortellisPage.waitForTimeout(500);
  await searchInput.fill(searchTerm);
  await cortellisPage.waitForTimeout(1000);
  await searchInput.press('Enter');

  console.log('⏳ Waiting for search results...');
  await cortellisPage.waitForTimeout(5000);

  // Take screenshot of search results overview
  const overviewScreenshot = path.join(resultDir, 'search_overview.png');
  await cortellisPage.screenshot({ path: overviewScreenshot, fullPage: true });
  console.log(`📸 Overview screenshot: ${overviewScreenshot}`);

  // Process each requested category
  for (let i = 0; i < categories.length; i++) {
    const category = categories[i];

    console.log(`\n--- Category ${i + 1}/${categories.length}: ${category} ---`);

    // Click "View results" for this category
    const found = await clickViewResults(cortellisPage, category);

    if (!found) {
      results.downloads.push({
        category,
        success: false,
        reason: 'Category not found or could not click'
      });

      // Go back to search results for next category
      await cortellisPage.goBack();
      await cortellisPage.waitForTimeout(2000);
      continue;
    }

    // We're now on the category detail page
    console.log(`✅ Navigated to ${category} results page`);
    console.log(`   URL: ${cortellisPage.url()}`);

    // Download the results
    const downloadResult = await downloadResults(cortellisPage, searchTerm, category, resultDir);
    results.downloads.push({
      category,
      ...downloadResult
    });

    // Go back to search results for next category
    console.log('🔙 Returning to search results...');
    await cortellisPage.goBack();
    await cortellisPage.waitForTimeout(3000);
  }

  // Save results metadata
  const metadataPath = path.join(resultDir, 'metadata.json');
  fs.writeFileSync(metadataPath, JSON.stringify(results, null, 2));
  console.log(`\n💾 Metadata saved: ${metadataPath}`);

  // Print summary
  console.log(`\n✅ Completed processing: "${searchTerm}"`);
  console.log(`   Successful downloads: ${results.downloads.filter(d => d.success).length}/${categories.length}`);

  return { success: true, ...results };
}

// Main execution
(async () => {
  const searchConfigs = parseArgs();
  console.log(`\n🎯 Will process ${searchConfigs.length} search(es):`);
  searchConfigs.forEach((config, i) => {
    console.log(`   ${i + 1}. ${config.searchTerm} (${config.categories.length} categories)`);
  });

  setupDirectories();

  const browser = await chromium.launch({
    headless: true,
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
      acceptDownloads: true, // Enable downloads
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

    // Navigate to Okta and launch Cortellis
    console.log('\n🌐 Navigating to Okta portal...');
    await page.goto(TARGET_URL, { waitUntil: 'networkidle', timeout: 60000 });
    await page.waitForTimeout(3000);

    // Find app search
    const appSearchInput = await page.$('input[type="search"]');
    if (!appSearchInput) throw new Error('Could not find app search');

    await appSearchInput.click();
    await appSearchInput.fill('cortellis drug discovery');
    await page.waitForTimeout(2000);

    // Launch Cortellis
    const cortellisLink = await page.locator('a:has-text("Cortellis Drug Discovery")').first();
    const newPagePromise = context.waitForEvent('page', { timeout: 15000 }).catch(() => null);
    await cortellisLink.click();

    let cortellisPage = await newPagePromise;
    if (!cortellisPage) {
      await page.waitForNavigation({ timeout: 10000 }).catch(() => {});
      cortellisPage = page;
    }

    console.log('⏳ Loading Cortellis...');
    await cortellisPage.waitForTimeout(5000);

    // Handle Cloudflare
    try {
      const pageTitle = await cortellisPage.title();
      if (pageTitle.includes('Just a moment') || pageTitle.includes('Verify')) {
        console.log('🔐 Cloudflare detected, waiting...');
        for (let i = 0; i < 15; i++) {
          await cortellisPage.waitForTimeout(3000);
          const currentTitle = await cortellisPage.title().catch(() => '');
          if (!currentTitle.includes('Just a moment') && !currentTitle.includes('Verify')) {
            break;
          }
        }
      }
    } catch (e) {
      await cortellisPage.waitForTimeout(2000);
    }

    await cortellisPage.waitForTimeout(3000);
    console.log('✅ Cortellis loaded');

    // Handle cookie consent
    try {
      const acceptButton = await cortellisPage.$('#onetrust-accept-btn-handler');
      if (acceptButton && await acceptButton.isVisible({ timeout: 2000 })) {
        await acceptButton.click({ force: true });
        await cortellisPage.waitForTimeout(2000);
      }
    } catch (e) {
      // No cookie banner
    }

    // Process each search
    const allResults = [];
    for (let i = 0; i < searchConfigs.length; i++) {
      const result = await processSearch(cortellisPage, searchConfigs[i], i + 1, searchConfigs.length);
      allResults.push(result);

      if (i < searchConfigs.length - 1) {
        console.log('\n⏳ Waiting 3 seconds before next search...');
        await cortellisPage.waitForTimeout(3000);
      }
    }

    // Generate final summary
    console.log('\n' + '='.repeat(70));
    console.log('📊 FINAL SUMMARY');
    console.log('='.repeat(70));

    const summary = {
      timestamp: new Date().toISOString(),
      totalSearches: searchConfigs.length,
      totalCategoriesRequested: searchConfigs.reduce((sum, c) => sum + c.categories.length, 0),
      totalDownloadsSuccessful: allResults.reduce((sum, r) => sum + r.downloads.filter(d => d.success).length, 0),
      results: allResults
    };

    const timestamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19).replace('T', '_');
    const summaryPath = path.join(OUTPUT_DIR, `summary_${timestamp}.json`);
    fs.writeFileSync(summaryPath, JSON.stringify(summary, null, 2));

    console.log(`\n💾 Summary: ${summaryPath}`);
    console.log(`📁 Results: ${OUTPUT_DIR}`);
    console.log(`\n✅ Completed ${summary.totalDownloadsSuccessful}/${summary.totalCategoriesRequested} downloads`);

    allResults.forEach((r, i) => {
      const successful = r.downloads.filter(d => d.success).length;
      console.log(`   ${i + 1}. ${r.searchTerm}: ${successful}/${r.requestedCategories.length} downloads`);
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
