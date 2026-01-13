const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');
const os = require('os');

// Configuration
const TARGET_URL = 'https://takeda.okta.com/';
// Use environment variable if set (from shell script), otherwise use current working directory
const WORK_DIR = process.env.CORTELLIS_WORK_DIR || process.cwd();
// Centralized Okta auth state (managed by okta-sso skill)
const AUTH_STATE_PATH = path.join(os.homedir(), '.okta', 'auth_state.json');
const OUTPUT_DIR = path.join(WORK_DIR, 'cortellis_playwright_result');

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

// Parse command line arguments
function parseArgs() {
  let args = process.argv.slice(2);

  // Filter out .js file paths (from run.js wrapper)
  args = args.filter(arg => !arg.endsWith('.js'));

  if (args.length === 0) {
    console.log('Usage:');
    console.log('  ./run-cortellis.sh <search_term1> [search_term2] ... [--category <category_name>]');
    console.log('  ./run-cortellis.sh --file <path_to_file> [--category <category_name>]');
    console.log('');
    console.log('Options:');
    console.log('  --category <name>   Specify which category to export (default: first available)');
    console.log('                      Examples: "Drugs & Biologics", "Clinical Studies", "Patents"');
    console.log('');
    console.log('Examples:');
    console.log('  ./run-cortellis.sh imatinib dasatinib');
    console.log('  ./run-cortellis.sh imatinib --category "Clinical Studies"');
    console.log('  ./run-cortellis.sh --file drugs.txt --category "Patents"');
    console.log('');
    console.log('Available categories:');
    console.log('  - Drugs & Biologics');
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
    console.log('');
    console.log('File format: One search term per line');
    process.exit(1);
  }

  let searchTerms = [];
  let category = null;

  // Check for --category flag
  const categoryIndex = args.findIndex(arg => arg === '--category' || arg === '-c');
  if (categoryIndex !== -1) {
    if (categoryIndex + 1 >= args.length) {
      console.error('Error: --category option requires a category name');
      process.exit(1);
    }
    category = args[categoryIndex + 1];
    // Remove category flag and value from args
    args.splice(categoryIndex, 2);
  }

  if (args[0] === '--file' || args[0] === '-f') {
    if (args.length < 2) {
      console.error('Error: --file option requires a file path');
      process.exit(1);
    }
    const filePath = args[1];
    if (!fs.existsSync(filePath)) {
      console.error(`Error: File not found: ${filePath}`);
      process.exit(1);
    }
    const content = fs.readFileSync(filePath, 'utf-8');
    searchTerms = content.split('\n')
      .map(line => line.trim())
      .filter(line => line.length > 0 && !line.startsWith('#'));
  } else {
    searchTerms = args;
  }

  return { searchTerms, category };
}

// Create output directory
function setupOutputDir() {
  if (!fs.existsSync(OUTPUT_DIR)) {
    fs.mkdirSync(OUTPUT_DIR, { recursive: true });
    console.log(`📁 Created output directory: ${OUTPUT_DIR}`);
  }
}

// Extract search results data
async function extractSearchResults(page) {
  try {
    await page.waitForTimeout(2000);

    const results = await page.evaluate(() => {
      const data = {
        timestamp: new Date().toISOString(),
        categories: []
      };

      // Try to find result cards/tiles
      const resultCards = document.querySelectorAll('[class*="card"], [class*="result"], [class*="tile"]');

      resultCards.forEach(card => {
        const text = card.innerText || '';
        const lines = text.split('\n').map(l => l.trim()).filter(l => l);

        if (lines.length >= 2) {
          // Try to extract count and category name
          const countMatch = lines[0].match(/\d+[,\d]*/);
          const category = lines[1] || lines[0];

          if (countMatch) {
            data.categories.push({
              count: countMatch[0].replace(/,/g, ''),
              category: category
            });
          }
        }
      });

      return data;
    });

    return results;
  } catch (e) {
    console.log('⚠️  Could not extract structured data:', e.message);
    return null;
  }
}

// Main search function
async function searchCortellis(cortellisPage, searchTerm, index, total, category = null) {
  console.log(`\n${'='.repeat(60)}`);
  console.log(`🔍 Searching ${index}/${total}: "${searchTerm}"`);
  if (category) {
    console.log(`📂 Target category: "${category}"`);
  }
  console.log('='.repeat(60));

  const searchTermSlug = searchTerm.replace(/[^a-z0-9]/gi, '_').toLowerCase();
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19).replace('T', '_');
  const categorySlug = category ? category.replace(/[^a-z0-9]/gi, '_') : 'default';
  const resultDir = path.join(OUTPUT_DIR, `${searchTermSlug}_${categorySlug}_${timestamp}`);

  // Create directory for this search
  fs.mkdirSync(resultDir, { recursive: true });

  // Find search input
  const cortellisSearchSelectors = [
    'input[type="search"]',
    'input[placeholder*="Search" i]',
    'input[name="search"]',
    'input[id*="search" i]',
    '[role="searchbox"]'
  ];

  let cortellisSearchInput = null;
  for (const selector of cortellisSearchSelectors) {
    try {
      const element = await cortellisPage.$(selector);
      if (element && await element.isVisible()) {
        cortellisSearchInput = element;
        break;
      }
    } catch (e) {
      // Continue
    }
  }

  if (!cortellisSearchInput) {
    console.log('⚠️  Could not find Cortellis search box');
    return { success: false, term: searchTerm };
  }

  // Clear previous search if any
  await cortellisSearchInput.click();
  await cortellisSearchInput.fill('');
  await cortellisPage.waitForTimeout(500);

  // Enter search term
  console.log(`📝 Entering search term: "${searchTerm}"`);
  await cortellisSearchInput.fill(searchTerm);
  await cortellisPage.waitForTimeout(1000);

  // Submit search
  await cortellisSearchInput.press('Enter');
  console.log('⏳ Waiting for search results...');
  await cortellisPage.waitForTimeout(5000);

  // Take screenshot of search results
  const searchScreenshotPath = path.join(resultDir, 'search_results.png');
  await cortellisPage.screenshot({ path: searchScreenshotPath, fullPage: true });
  console.log(`📸 Search screenshot saved: ${searchScreenshotPath}`);

  // Extract structured data
  console.log('📊 Extracting search results data...');
  const resultsData = await extractSearchResults(cortellisPage);

  // Click "View Results" to navigate to detailed results page
  console.log('🔍 Looking for "View Results" button...');

  let viewResultsClicked = false;

  if (category) {
    // Try to find the specific category card and click its "View results" button
    console.log(`📂 Searching for category: "${category}"`);

    try {
      // Find all text elements that might contain the category name
      const categoryElement = cortellisPage.locator(`text="${category}"`).first();

      if (await categoryElement.isVisible({ timeout: 5000 })) {
        console.log(`✅ Found category: "${category}"`);

        // Find the parent card/container
        const card = categoryElement.locator('xpath=ancestor::div[contains(@class, "card") or contains(@class, "tile") or contains(@class, "result")]').first();

        // Try to find "View results" button within this card
        const viewResultsInCard = card.locator('button:has-text("View results"), a:has-text("View results")').first();

        if (await viewResultsInCard.isVisible({ timeout: 3000 })) {
          console.log(`✅ Found "View Results" button for "${category}", clicking...`);
          await viewResultsInCard.click();
          viewResultsClicked = true;
          await cortellisPage.waitForTimeout(5000);
        } else {
          console.log(`⚠️  Could not find "View Results" button for "${category}"`);
        }
      } else {
        console.log(`⚠️  Category "${category}" not found in search results`);
      }
    } catch (e) {
      console.log(`⚠️  Error finding category "${category}": ${e.message}`);
    }
  }

  // If category not specified or not found, use the first "View Results" button
  if (!viewResultsClicked) {
    if (category) {
      console.log('⚠️  Falling back to first "View Results" button...');
    }

    const viewResultsSelectors = [
      'button:has-text("View Results")',
      'a:has-text("View Results")',
      'button:has-text("View results")',
      'a:has-text("View results")',
      '[role="button"]:has-text("View")',
      'button:has-text("View")'
    ];

    for (const selector of viewResultsSelectors) {
      try {
        const element = cortellisPage.locator(selector).first();
        if (await element.isVisible({ timeout: 3000 })) {
          console.log('✅ Found "View Results" button, clicking...');
          await element.click();
          viewResultsClicked = true;
          await cortellisPage.waitForTimeout(5000);
          break;
        }
      } catch (e) {
        // Continue to next selector
      }
    }
  }

  if (!viewResultsClicked) {
    console.log('⚠️  Could not find "View Results" button, checking if already on results page...');
  } else {
    console.log('✅ Navigated to results page');
  }

  // Take screenshot of results page
  const resultsScreenshotPath = path.join(resultDir, 'results_page.png');
  await cortellisPage.screenshot({ path: resultsScreenshotPath, fullPage: true });
  console.log(`📸 Results page screenshot saved: ${resultsScreenshotPath}`);

  // Wait a bit for any dynamic content to load
  await cortellisPage.waitForTimeout(2000);

  // Try to find all clickable buttons in the page to help debug
  console.log('🔍 Scanning page for export/menu buttons...');

  // Look for any buttons or icons that might trigger export
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
      const elements = await cortellisPage.locator(selector).all();

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
              const downloadPromise = cortellisPage.waitForEvent('download', { timeout: 30000 }).catch(() => null);

              await element.click();
              await cortellisPage.waitForTimeout(2000);

              // Check if download started immediately
              const download = await downloadPromise;
              if (download) {
                const downloadFileName = download.suggestedFilename() || `export_${searchTerm.replace(/[^a-z0-9]/gi, '_')}.xlsx`;
                const downloadPath = path.join(resultDir, downloadFileName);
                await download.saveAs(downloadPath);
                console.log(`✅ Export downloaded directly: ${downloadPath}`);
                exportInitiated = true;
                break;
              }

              // Take screenshot after clicking to see what opened
              const afterClickPath = path.join(resultDir, 'after_button_click.png');
              await cortellisPage.screenshot({ path: afterClickPath, fullPage: true });
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
                  const menuElement = cortellisPage.locator(menuSelector).first();
                  if (await menuElement.isVisible({ timeout: 2000 })) {
                    console.log('✅ Found "Export" in menu, clicking...');

                    const downloadPromise2 = cortellisPage.waitForEvent('download', { timeout: 30000 }).catch(() => null);

                    await menuElement.click();
                    await cortellisPage.waitForTimeout(2000);

                    // Check if download started
                    const download2 = await downloadPromise2;
                    if (download2) {
                      const downloadFileName = download2.suggestedFilename() || `export_${searchTerm.replace(/[^a-z0-9]/gi, '_')}.xlsx`;
                      const downloadPath = path.join(resultDir, downloadFileName);
                      await download2.saveAs(downloadPath);
                      console.log(`✅ Export downloaded: ${downloadPath}`);
                      exportInitiated = true;
                      break;
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
                        const confirmElements = await cortellisPage.locator(confirmSelector).all();
                        for (const confirmElement of confirmElements) {
                          if (await confirmElement.isVisible({ timeout: 1000 })) {
                            console.log(`✅ Found confirmation button, clicking...`);

                            const downloadPromise3 = cortellisPage.waitForEvent('download', { timeout: 30000 }).catch(() => null);

                            await confirmElement.click();
                            await cortellisPage.waitForTimeout(3000);

                            const download3 = await downloadPromise3;
                            if (download3) {
                              const downloadFileName = download3.suggestedFilename() || `export_${searchTerm.replace(/[^a-z0-9]/gi, '_')}.xlsx`;
                              const downloadPath = path.join(resultDir, downloadFileName);
                              await download3.saveAs(downloadPath);
                              console.log(`✅ Export downloaded: ${downloadPath}`);
                              exportInitiated = true;
                              break;
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
  }

  // Take final screenshot
  const screenshotPath = path.join(resultDir, 'final_state.png');
  await cortellisPage.screenshot({ path: screenshotPath, fullPage: true });
  console.log(`📸 Final screenshot saved: ${screenshotPath}`);

  if (resultsData) {
    resultsData.searchTerm = searchTerm;
    resultsData.url = cortellisPage.url();

    const jsonPath = path.join(resultDir, 'results.json');
    fs.writeFileSync(jsonPath, JSON.stringify(resultsData, null, 2));
    console.log(`💾 Data saved: ${jsonPath}`);

    // Print summary
    console.log('\n📋 Results Summary:');
    console.log(`   Search term: ${searchTerm}`);
    console.log(`   Categories found: ${resultsData.categories.length}`);
    resultsData.categories.forEach(cat => {
      console.log(`   - ${cat.category}: ${cat.count}`);
    });
  }

  // Save HTML content
  const htmlPath = path.join(resultDir, 'results.html');
  const htmlContent = await cortellisPage.content();
  fs.writeFileSync(htmlPath, htmlContent);
  console.log(`💾 HTML saved: ${htmlPath}`);

  console.log(`✅ Completed search for: "${searchTerm}"`);
  return { success: true, term: searchTerm, data: resultsData, dir: resultDir };
}

// Main execution
(async () => {
  // Validate Okta session before starting
  validateAuthState();

  const { searchTerms, category } = parseArgs();
  console.log(`\n🎯 Will search for ${searchTerms.length} term(s):`);
  searchTerms.forEach((term, i) => console.log(`   ${i + 1}. ${term}`));
  if (category) {
    console.log(`📂 Target category: "${category}"`);
  }

  setupOutputDir();

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
    console.log('\n🔐 Loading authentication state from:', AUTH_STATE_PATH);
    const context = await browser.newContext({
      storageState: AUTH_STATE_PATH,
      userAgent: 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
      viewport: { width: 1920, height: 1080 },
      locale: 'en-US',
      timezoneId: 'America/New_York',
      permissions: [],
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

    // Navigate to Okta portal
    console.log('\n🌐 Navigating to Takeda Okta portal...');
    await page.goto(TARGET_URL, { waitUntil: 'networkidle', timeout: 60000 });
    console.log('✅ Portal loaded');
    await page.waitForTimeout(3000);

    // Find and use app search
    console.log('\n🔍 Looking for app launcher...');
    const appSearchSelectors = [
      'input[type="search"]',
      'input[placeholder*="Search" i]',
      'input[aria-label*="Search" i]',
      '[role="searchbox"]'
    ];

    let appSearchInput = null;
    for (const selector of appSearchSelectors) {
      const element = await page.$(selector);
      if (element && await element.isVisible()) {
        appSearchInput = element;
        break;
      }
    }

    if (!appSearchInput) {
      throw new Error('Could not find app search box');
    }

    console.log('📝 Searching for Cortellis Drug Discovery...');
    await appSearchInput.click();
    await appSearchInput.fill('cortellis drug discovery');
    await page.waitForTimeout(2000);

    // Click on Cortellis app
    const cortellisLink = await page.locator('a:has-text("Cortellis Drug Discovery")').first();
    if (!await cortellisLink.isVisible({ timeout: 5000 })) {
      throw new Error('Could not find Cortellis app link');
    }

    console.log('🖱️  Launching Cortellis Drug Discovery...');
    const newPagePromise = context.waitForEvent('page', { timeout: 15000 }).catch(() => null);
    await cortellisLink.click();

    let cortellisPage = await newPagePromise;
    if (!cortellisPage) {
      await page.waitForNavigation({ timeout: 10000 }).catch(() => {});
      cortellisPage = page;
    }

    console.log('⏳ Waiting for Cortellis to load...');
    await cortellisPage.waitForTimeout(5000);

    // Handle Cloudflare if needed
    let pageTitle;
    try {
      pageTitle = await cortellisPage.title();
    } catch (e) {
      // Page may have navigated, wait and retry
      await cortellisPage.waitForTimeout(2000);
      pageTitle = await cortellisPage.title();
    }

    if (pageTitle.includes('Just a moment') || pageTitle.includes('Verify')) {
      console.log('🔐 Cloudflare verification detected, waiting...');

      for (let i = 0; i < 15; i++) {
        await cortellisPage.waitForTimeout(3000);

        try {
          const currentTitle = await cortellisPage.title();

          if (!currentTitle.includes('Just a moment') && !currentTitle.includes('Verify')) {
            console.log('✅ Cloudflare verification passed!');
            break;
          }
        } catch (e) {
          // Page may be navigating, continue waiting
          console.log(`   Waiting... (${i + 1}/15)`);
        }

        if (i === 14) {
          throw new Error('Unable to bypass Cloudflare verification');
        }
      }
    }

    // Wait for page to be stable
    await cortellisPage.waitForTimeout(3000);

    console.log('✅ Cortellis loaded successfully');
    console.log('URL:', cortellisPage.url());
    await cortellisPage.waitForTimeout(2000);

    // Handle cookie consent banner
    console.log('🍪 Checking for cookie consent banner...');
    try {
      const acceptButtonSelectors = [
        'button:has-text("Accept")',
        'button:has-text("Accept All")',
        '#onetrust-accept-btn-handler',
        '.onetrust-accept-btn-handler',
        'button[id*="accept"]',
        'button[class*="accept"]'
      ];

      for (const selector of acceptButtonSelectors) {
        try {
          const button = await cortellisPage.$(selector);
          if (button && await button.isVisible({ timeout: 2000 })) {
            console.log('✅ Found cookie consent button, clicking...');
            await button.click({ force: true });
            await cortellisPage.waitForTimeout(2000);
            break;
          }
        } catch (e) {
          // Continue to next selector
        }
      }
    } catch (e) {
      console.log('⚠️  No cookie consent banner found (or already handled)');
    }

    // Perform searches iteratively
    const results = [];
    for (let i = 0; i < searchTerms.length; i++) {
      const result = await searchCortellis(cortellisPage, searchTerms[i], i + 1, searchTerms.length, category);
      results.push(result);

      // Wait between searches to avoid rate limiting
      if (i < searchTerms.length - 1) {
        console.log('\n⏳ Waiting 3 seconds before next search...');
        await cortellisPage.waitForTimeout(3000);
      }
    }

    // Generate summary report
    console.log('\n' + '='.repeat(60));
    console.log('📊 FINAL SUMMARY');
    console.log('='.repeat(60));

    const summaryPath = path.join(OUTPUT_DIR, `summary_${new Date().toISOString().replace(/[:.]/g, '-')}.json`);
    const summary = {
      timestamp: new Date().toISOString(),
      totalSearches: searchTerms.length,
      successful: results.filter(r => r.success).length,
      failed: results.filter(r => !r.success).length,
      results: results
    };

    fs.writeFileSync(summaryPath, JSON.stringify(summary, null, 2));
    console.log(`\n💾 Summary report saved: ${summaryPath}`);
    console.log(`📁 All results saved to: ${OUTPUT_DIR}`);

    console.log('\n✅ Searches completed:');
    results.forEach((r, i) => {
      console.log(`   ${i + 1}. ${r.term}: ${r.success ? '✅ Success' : '❌ Failed'}`);
    });

    // Close Cortellis page
    if (cortellisPage !== page) {
      await cortellisPage.close();
    }

  } catch (error) {
    console.error('\n❌ Error:', error.message);
    console.error(error.stack);

    try {
      const allPages = await browser.contexts()[0]?.pages();
      if (allPages && allPages.length > 0) {
        for (let i = 0; i < allPages.length; i++) {
          await allPages[i].screenshot({ path: `/tmp/cortellis-error-page${i}.png`, fullPage: true });
          console.log(`📸 Error screenshot saved: /tmp/cortellis-error-page${i}.png`);
        }
      }
    } catch (e) {
      // Ignore screenshot errors
    }

    process.exit(1);
  } finally {
    await browser.close();
    console.log('\n✅ Browser closed');
  }
})();
