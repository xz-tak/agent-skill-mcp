"""
NCBI Paper Query - Paper Downloader Module

Downloads papers from various sources with subscription-only mode support.
Supports institutional access via OpenAthens/Shibboleth SSO authentication.
"""

import os
import re
import json
import base64
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup

from ..config import logger


class PaperDownloader:
    """Download papers from various sources with subscription-only mode support.

    Supports institutional access via credentials loaded from .env file.
    Uses OpenAthens/Shibboleth SSO authentication for major publishers.
    """

    PMC_OA_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi"
    PMC_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/pmc/"

    # Publisher-specific institutional access URLs and patterns
    PUBLISHER_AUTH_PATTERNS = {
        "elsevier": {
            "domains": ["sciencedirect.com", "cell.com", "thelancet.com", "linkinghub.elsevier.com"],
            "institutional_link_selectors": [
                'a:has-text("Access through your institution")',
                'a:has-text("Institutional access")',
                'button:has-text("Access through your institution")',
                '#institution-access',
                '.institution-login',
            ],
            "openathens_url": "https://auth.elsevier.com/ShibAuth/institutionLogin.asp",
        },
        "springer": {
            "domains": ["springer.com", "nature.com", "link.springer.com"],
            "institutional_link_selectors": [
                'a:has-text("Access through your institution")',
                'a:has-text("Log in via an institution")',
                'button:has-text("institutional")',
                '.c-login-institutions__link',
            ],
            "openathens_url": "https://wayf.springernature.com/",
        },
        "wiley": {
            "domains": ["wiley.com", "onlinelibrary.wiley.com"],
            "institutional_link_selectors": [
                'a:has-text("Institutional Login")',
                'a:has-text("Shibboleth")',
                'button:has-text("institutional")',
                '.institutional-login',
            ],
            "openathens_url": "https://www.onlinelibrary.wiley.com/action/ssostart",
        },
        "oxford": {
            "domains": ["academic.oup.com", "oup.com"],
            "institutional_link_selectors": [
                'a:has-text("Sign in via your Institution")',
                'a:has-text("Institutional access")',
                '.institutional-signin',
            ],
            "openathens_url": "https://academic.oup.com/my-account/oauth/shibboleth",
        },
        "taylor_francis": {
            "domains": ["tandfonline.com"],
            "institutional_link_selectors": [
                'a:has-text("Login with your institution")',
                '.shibboleth-login',
            ],
            "openathens_url": "https://www.tandfonline.com/action/ssostart",
        },
        "acs": {
            "domains": ["pubs.acs.org"],
            "institutional_link_selectors": [
                'a:has-text("Find my institution")',
                '.institutional-access',
            ],
            "openathens_url": "https://pubs.acs.org/action/ssostart",
        },
        "aaas": {
            "domains": ["science.org", "sciencemag.org", "stm.sciencemag.org"],
            "institutional_link_selectors": [
                'a:has-text("Check Access")',
                'button:has-text("Check Access")',
                'a:has-text("Institution")',
                'a:has-text("Institutional")',
                '.institution-access',
                'a[href*="institution"]',
                '.article-controls a[href*="pdf"]',
                'a:has-text("Full Text PDF")',
                'a:has-text("View PDF")',
            ],
            "openathens_url": "https://www.science.org/action/ssostart",
        },
        "biorxiv": {
            "domains": ["biorxiv.org", "medrxiv.org"],
            "institutional_link_selectors": [],  # Free access - no auth needed
            "pdf_selectors": [
                'a:has-text("Download PDF")',
                'a[href*=".full.pdf"]',
                'a[href*="/pdf/"]',
            ],
            "openathens_url": None,  # Free access
            "free_access": True,
        },
        "wolters_kluwer": {
            "domains": ["journals.lww.com", "lww.com"],
            "institutional_link_selectors": [
                'a:has-text("Log In")',
                'a:has-text("Sign in")',
                'a:has-text("Institutional Access")',
                'button:has-text("Log In")',
                '.login-link',
            ],
            "pdf_selectors": [
                'a.ejp-article-pdf-link',
                'a[href*="/Fulltext.pdf"]',
                'a:has-text("Full Text PDF")',
                'a:has-text("Download PDF")',
                '.article-tools a[href*="pdf"]',
            ],
            "openathens_url": "https://journals.lww.com/pages/login.aspx",
        },
        "bmj": {
            "domains": ["bmj.com", "gut.bmj.com", "thorax.bmj.com"],
            "institutional_link_selectors": [
                'a:has-text("Log in")',
                'a:has-text("Institutional access")',
                'a:has-text("via your institution")',
                '.login-link',
            ],
            "pdf_selectors": [
                'a.article-pdf-download',
                'a[data-track-action="download pdf"]',
                'a:has-text("PDF")',
                'a[href*=".full.pdf"]',
                '.article-tools a[href*="pdf"]',
            ],
            "openathens_url": "https://sso.bmj.com/",
        },
        "portland_press": {
            "domains": ["portlandpress.com"],
            "institutional_link_selectors": [
                'a:has-text("Sign in via your Institution")',
                'a:has-text("Institutional access")',
            ],
            "pdf_selectors": [
                'a.article-pdfLink',
                'a[href*="/pdf/"]',
                'a:has-text("Full Text PDF")',
                'a:has-text("Download PDF")',
                '.article-tools a[href*="pdf"]',
            ],
            "openathens_url": None,
        },
    }

    # Institution-specific identifiers for OpenAthens/Shibboleth
    INSTITUTION_IDENTIFIERS = [
        "Takeda",
        "takeda.com",
        "Takeda Pharmaceutical",
        "Takeda Pharmaceuticals",
    ]

    def __init__(self, download_dir: str = "downloads", subscription_only: bool = False):
        self.download_dir = Path(download_dir)
        self.download_dir.mkdir(parents=True, exist_ok=True)
        self.subscription_only = subscription_only
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
        })
        self._playwright_browser = None
        self._playwright_context = None
        self._is_authenticated = False

        # Load credentials from .env
        self.tak_account = os.getenv("TAK_ACCOUNT")
        self.tak_key = os.getenv("TAK_KEY")

        if self.tak_account and self.tak_key:
            logger.info(f"Institutional credentials loaded for: {self.tak_account}")
        else:
            logger.warning("Institutional credentials not found in .env file. Institutional access may be limited.")

    def get_paper_dir(self, pmid: str) -> Path:
        """Get directory for a paper's downloads."""
        paper_dir = self.download_dir / f"PMID_{pmid}"
        paper_dir.mkdir(parents=True, exist_ok=True)
        return paper_dir

    def has_free_access(self, pmc_id: str) -> bool:
        """Check if paper has free access via PMC."""
        if not pmc_id:
            return False

        try:
            # PMC OA API now returns XML (not JSON)
            params = {"id": pmc_id}
            response = requests.get(self.PMC_OA_URL, params=params, timeout=10)

            if response.status_code == 200:
                # Parse XML response
                soup = BeautifulSoup(response.text, 'xml')
                # Check for PDF link in the response
                links = soup.find_all('link')
                for link in links:
                    if link.get('format') == 'pdf':
                        return True

                # Also check if there's a record (indicates OA availability)
                records = soup.find_all('record')
                if records:
                    return True

            return False

        except Exception as e:
            logger.debug(f"Free access check failed: {e}")
            return False

    def check_pmc_access(self, pmc_id: str) -> Optional[str]:
        """Check if paper is available in PMC Open Access."""
        if not pmc_id:
            return None

        try:
            # Query PMC OA service (XML format)
            params = {"id": pmc_id}
            response = requests.get(self.PMC_OA_URL, params=params, timeout=10)

            if response.status_code == 200:
                # Parse XML response
                soup = BeautifulSoup(response.text, 'xml')
                # Look for PDF link
                links = soup.find_all('link')
                for link in links:
                    if link.get('format') == 'pdf':
                        href = link.get('href')
                        if href:
                            # Convert FTP to HTTPS if needed
                            if href.startswith('ftp://'):
                                href = href.replace('ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/',
                                                   'https://www.ncbi.nlm.nih.gov/pmc/articles/')
                            return href

            # Try direct PMC URL as fallback
            pmc_clean = pmc_id.replace("PMC", "")
            return f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_clean}/pdf/"

        except Exception as e:
            logger.debug(f"PMC access check failed: {e}")
            return None

    def get_doi_url(self, doi: str) -> Optional[str]:
        """Resolve DOI to publisher URL."""
        if not doi:
            return None

        try:
            response = self.session.head(
                f"https://doi.org/{doi}",
                allow_redirects=True,
                timeout=10
            )
            return response.url
        except:
            return f"https://doi.org/{doi}"

    def download_pdf(self, url: str, save_path: Path) -> bool:
        """Download PDF from URL."""
        try:
            response = self.session.get(url, timeout=30, allow_redirects=True)

            if response.status_code == 200:
                content_type = response.headers.get('content-type', '')

                if 'pdf' in content_type.lower() or url.endswith('.pdf'):
                    with open(save_path, 'wb') as f:
                        f.write(response.content)

                    # Verify it's a valid PDF
                    if save_path.stat().st_size > 1000:  # At least 1KB
                        return True
                    else:
                        save_path.unlink()

            return False

        except Exception as e:
            logger.debug(f"PDF download failed: {e}")
            return False

    def _detect_publisher(self, url: str) -> Optional[str]:
        """Detect the publisher from URL."""
        url_lower = url.lower()
        for publisher, config in self.PUBLISHER_AUTH_PATTERNS.items():
            for domain in config["domains"]:
                if domain in url_lower:
                    return publisher
        return None

    def _get_persistent_context_dir(self) -> Path:
        """Get directory for persistent browser context to retain auth cookies."""
        context_dir = self.download_dir / ".browser_context"
        context_dir.mkdir(parents=True, exist_ok=True)
        return context_dir

    def _handle_sso_login(self, page, publisher: str = None) -> bool:
        """
        Handle SSO authentication via OpenAthens/Shibboleth.

        Returns True if authentication succeeded.
        """
        if not self.tak_account or not self.tak_key:
            logger.warning("Institutional credentials not available for SSO login")
            return False

        try:
            # Wait for any redirects to complete
            page.wait_for_load_state("networkidle", timeout=10000)
            current_url = page.url.lower()

            # Check if we're on a Shibboleth/OpenAthens login page
            shibboleth_indicators = [
                "shibboleth", "openathens", "idp", "wayf",
                "login", "signin", "authenticate", "sso"
            ]

            is_auth_page = any(ind in current_url for ind in shibboleth_indicators)

            if not is_auth_page:
                # Check page content for auth indicators
                page_text = page.content().lower()
                is_auth_page = any(ind in page_text for ind in [
                    "sign in", "log in", "authenticate", "institution"
                ])

            if is_auth_page:
                logger.info("Detected authentication page, attempting institutional SSO login...")

                # Try to find and click institution search/select
                institution_selectors = [
                    'input[placeholder*="institution"]',
                    'input[placeholder*="search"]',
                    'input[name*="institution"]',
                    'input[id*="institution"]',
                    '#institution-search',
                    '.institution-search',
                    'input[type="search"]',
                ]

                for selector in institution_selectors:
                    try:
                        search_input = page.locator(selector).first
                        if search_input.is_visible(timeout=2000):
                            search_input.fill("Takeda")
                            page.wait_for_timeout(1000)

                            # Look for Takeda in dropdown/results
                            for identifier in self.INSTITUTION_IDENTIFIERS:
                                try:
                                    takeda_option = page.locator(f'text="{identifier}"').first
                                    if takeda_option.is_visible(timeout=2000):
                                        takeda_option.click()
                                        page.wait_for_load_state("networkidle", timeout=15000)
                                        break
                                except:
                                    continue
                            break
                    except:
                        continue

                # Now look for username/email and password fields
                page.wait_for_timeout(1000)

                # Common username/email field selectors
                username_selectors = [
                    'input[type="email"]',
                    'input[name="username"]',
                    'input[name="email"]',
                    'input[name="j_username"]',
                    'input[id*="username"]',
                    'input[id*="email"]',
                    'input[placeholder*="email"]',
                    'input[placeholder*="username"]',
                    '#username',
                    '#email',
                ]

                password_selectors = [
                    'input[type="password"]',
                    'input[name="password"]',
                    'input[name="j_password"]',
                    'input[id*="password"]',
                    '#password',
                ]

                # Try to fill username
                for selector in username_selectors:
                    try:
                        username_input = page.locator(selector).first
                        if username_input.is_visible(timeout=2000):
                            username_input.fill(self.tak_account)
                            logger.debug(f"Filled username with {self.tak_account}")
                            break
                    except:
                        continue

                # Try to fill password
                for selector in password_selectors:
                    try:
                        password_input = page.locator(selector).first
                        if password_input.is_visible(timeout=2000):
                            password_input.fill(self.tak_key)
                            logger.debug("Filled password field")
                            break
                    except:
                        continue

                # Look for and click submit/login button
                submit_selectors = [
                    'button[type="submit"]',
                    'input[type="submit"]',
                    'button:has-text("Sign in")',
                    'button:has-text("Log in")',
                    'button:has-text("Login")',
                    'button:has-text("Submit")',
                    'button:has-text("Continue")',
                    '#submit',
                    '.submit-button',
                    '.login-button',
                ]

                for selector in submit_selectors:
                    try:
                        submit_btn = page.locator(selector).first
                        if submit_btn.is_visible(timeout=2000):
                            submit_btn.click()
                            page.wait_for_load_state("networkidle", timeout=30000)
                            logger.info("Submitted login form")
                            break
                    except:
                        continue

                # Check if login was successful by looking for auth cookies or redirects
                page.wait_for_timeout(2000)
                self._is_authenticated = True
                return True

            return False

        except Exception as e:
            logger.error(f"Institutional SSO login failed: {e}")
            return False

    def _navigate_to_institutional_access(self, page, publisher: str) -> bool:
        """Navigate to institutional access link on publisher page."""
        if publisher not in self.PUBLISHER_AUTH_PATTERNS:
            return False

        config = self.PUBLISHER_AUTH_PATTERNS[publisher]
        selectors = config.get("institutional_link_selectors", [])

        for selector in selectors:
            try:
                inst_link = page.locator(selector).first
                if inst_link.is_visible(timeout=3000):
                    inst_link.click()
                    page.wait_for_load_state("networkidle", timeout=15000)
                    logger.info(f"Clicked institutional access link for {publisher}")
                    return True
            except:
                continue

        return False

    def download_via_playwright(self, doi_url: str, save_path: Path, use_institutional_auth: bool = True) -> bool:
        """
        Download paper via Playwright browser automation with institutional access support.

        Uses institutional credentials from .env for OpenAthens/Shibboleth SSO.

        Args:
            doi_url: The DOI URL or publisher URL to download from
            save_path: Path to save the downloaded PDF
            use_institutional_auth: Whether to attempt institutional login if needed

        Returns:
            True if download succeeded, False otherwise
        """
        try:
            from playwright.sync_api import sync_playwright
        except ImportError:
            logger.warning("Playwright not installed. Run: pip install playwright && playwright install")
            return False

        try:
            publisher = self._detect_publisher(doi_url)
            logger.info(f"Attempting download from {publisher or 'unknown publisher'}: {doi_url}")

            with sync_playwright() as p:
                # Use persistent context to retain authentication cookies across downloads
                context_dir = self._get_persistent_context_dir()

                browser = p.chromium.launch(headless=True)
                context = browser.new_context(
                    accept_downloads=True,
                    user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36",
                    viewport={"width": 1920, "height": 1080},
                    java_script_enabled=True,
                    locale="en-US",
                    timezone_id="America/New_York",
                )

                # Load cookies from persistent storage if available
                cookies_file = context_dir / "cookies.json"
                if cookies_file.exists():
                    try:
                        with open(cookies_file, 'r') as f:
                            cookies = json.load(f)
                        context.add_cookies(cookies)
                        logger.debug("Loaded saved authentication cookies")
                    except Exception as e:
                        logger.debug(f"Failed to load cookies: {e}")

                page = context.new_page()

                # Navigate to DOI URL
                page.goto(doi_url, timeout=60000)
                page.wait_for_load_state("networkidle")

                # Check if we can access the PDF directly (already authenticated or open access)
                success = self._try_download_pdf(page, save_path)
                if success:
                    self._save_cookies(context, cookies_file)
                    browser.close()
                    return True

                # If not, attempt institutional authentication
                if use_institutional_auth and self.tak_account and self.tak_key:
                    logger.info("PDF not directly available, attempting institutional authentication...")

                    # Try to find and click institutional access link
                    if publisher:
                        self._navigate_to_institutional_access(page, publisher)

                    # Handle SSO login
                    self._handle_sso_login(page, publisher)

                    # After authentication, navigate back to the paper and try download again
                    page.goto(doi_url, timeout=60000)
                    page.wait_for_load_state("networkidle")

                    success = self._try_download_pdf(page, save_path)
                    if success:
                        self._save_cookies(context, cookies_file)
                        browser.close()
                        return True

                # Try alternative download methods
                success = self._try_alternative_download_methods(page, doi_url, save_path)
                if success:
                    self._save_cookies(context, cookies_file)
                    browser.close()
                    return True

                browser.close()
                return False

        except Exception as e:
            logger.debug(f"Playwright download failed: {e}")
            return False

    def _is_supplementary_link(self, href: str, text: str) -> bool:
        """Check if a link points to supplementary material instead of main paper."""
        if not href and not text:
            return False

        href_lower = (href or "").lower()
        text_lower = (text or "").lower()

        # Patterns that indicate supplementary/supporting material
        supplementary_patterns = [
            "supplement", "supporting", "appendix", "extended",
            "additional", "extra", "si_", "s1_", "s2_", "s3_",
            "supp_", "suppl", "esm", "online_resource",
            "table_s", "figure_s", "data_s", "file_s",
            "mmc", "multimedia", "source_data", "source-data",
        ]

        for pattern in supplementary_patterns:
            if pattern in href_lower or pattern in text_lower:
                return True

        return False

    def _try_download_pdf(self, page, save_path: Path) -> bool:
        """Try to find and download the MAIN PAPER PDF from the current page.

        This method specifically targets the main article PDF and filters out
        supplementary materials, supporting information, and other auxiliary files.
        """
        # PRIORITY 1: Publisher-specific main article PDF selectors (most reliable)
        main_paper_selectors = [
            # Nature/Springer - main article PDF
            'a[data-track-action="download pdf"]',
            '.c-pdf-download__link',
            'a[href*="/content/pdf/"][href$=".pdf"]:not([href*="supplement"]):not([href*="esm"])',

            # Elsevier/ScienceDirect - main article PDF
            '.PdfDropDownMenu a[href*="pdf"]',
            'a.pdf-download-btn',
            '#pdfLink',
            'a[href*="/pdfdirect/"]:not([href*="mmc"])',
            '.article-tools a[href*="pdf"]',

            # Wiley - main article PDF
            '.epub-section__item a[href*="epdf"]',
            '.article-tool-pdf a',
            'a[href*="/doi/epdf/"]:not([href*="supp"])',
            'a[href*="/doi/pdf/"]:not([href*="supp"])',

            # Oxford Academic - main article PDF
            '.article-pdfLink',
            '.accessbar__link--pdf',
            'a.al-link[href*="pdf"]',

            # Cell Press
            'a.pdfLink[href*="pdf"]',
            '.article-tools a.pdf',

            # PNAS/Science
            'a.article-dl-pdf-link',
            '.toolbar-link-pdf a',

            # Wolters Kluwer (LWW)
            'a.ejp-article-pdf-link',
            'a[href*="/Fulltext.pdf"]',

            # BMJ/Gut
            'a.article-pdf-download',

            # Portland Press
            'a.article-pdfLink',

            # bioRxiv/medRxiv (preprints - free access)
            'a[href*=".full.pdf"]',

            # Generic main article selectors (text-based, more specific)
            'a:has-text("Download PDF"):not(:has-text("Supplement"))',
            'a:has-text("Full Text PDF")',
            'a:has-text("Article PDF")',
            'a:has-text("PDF ("):not(:has-text("Supplement"))',  # "PDF (1.2 MB)" format
            'a:has-text("View PDF"):not(:has-text("Supplement"))',
            'a:has-text("Get PDF"):not(:has-text("Supplement"))',
        ]

        # PRIORITY 2: Generic PDF selectors (less reliable, may catch supplements)
        generic_selectors = [
            'a[href*="/pdf/"][href$=".pdf"]',
            'a[href*=".pdf"]:not([href*="supplement"]):not([href*="supp"]):not([href*="mmc"])',
            'button:has-text("PDF"):not(:has-text("Supplement"))',
        ]

        # Try main paper selectors first
        for selector in main_paper_selectors:
            if self._try_selector_download(page, selector, save_path, check_supplementary=True):
                return True

        # Fall back to generic selectors with strict filtering
        for selector in generic_selectors:
            if self._try_selector_download(page, selector, save_path, check_supplementary=True):
                return True

        return False

    def _try_selector_download(self, page, selector: str, save_path: Path, check_supplementary: bool = True) -> bool:
        """Try to download PDF using a specific selector."""
        try:
            # Get all matching elements
            links = page.locator(selector).all()

            for link in links[:5]:  # Try up to 5 matches
                try:
                    if not link.is_visible(timeout=1000):
                        continue

                    href = link.get_attribute('href') or ""
                    text = link.text_content() or ""

                    # Skip supplementary materials
                    if check_supplementary and self._is_supplementary_link(href, text):
                        logger.debug(f"Skipping supplementary link: {href[:50]}...")
                        continue

                    # Skip very small files (likely icons/thumbnails)
                    # Skip links without PDF indicators
                    if not any(ind in href.lower() for ind in ['pdf', 'download', 'epdf']):
                        if not any(ind in text.lower() for ind in ['pdf', 'download']):
                            continue

                    logger.debug(f"Trying PDF link: {href[:80]}... text: {text[:50]}")

                    # Try to download
                    try:
                        with page.expect_download(timeout=90000) as download_info:
                            link.click()
                        download = download_info.value
                        suggested_filename = download.suggested_filename.lower()

                        # Check if downloaded file is supplementary
                        if check_supplementary and self._is_supplementary_link(suggested_filename, ""):
                            logger.warning(f"Downloaded file appears to be supplementary: {suggested_filename}")
                            download.delete()
                            continue

                        download.save_as(save_path)

                        # Validate it's a real PDF with reasonable size
                        if save_path.exists():
                            size = save_path.stat().st_size
                            if size > 50000:  # Main papers are usually > 50KB
                                with open(save_path, 'rb') as f:
                                    header = f.read(4)
                                    if header == b'%PDF':
                                        logger.info(f"Downloaded main paper PDF ({size/1024:.0f} KB) via {selector[:50]}")
                                        return True
                                    else:
                                        logger.debug(f"Downloaded file is not a valid PDF")
                                        save_path.unlink()
                            else:
                                logger.debug(f"Downloaded file too small ({size} bytes), likely not main paper")
                                save_path.unlink()
                    except Exception as e:
                        logger.debug(f"Download attempt failed: {e}")
                        continue

                except Exception as e:
                    logger.debug(f"Link processing failed: {e}")
                    continue

        except Exception as e:
            logger.debug(f"Selector {selector[:50]} failed: {e}")

        return False

    def _save_cookies(self, context, cookies_file: Path):
        """Save browser cookies to file for persistence."""
        try:
            cookies = context.cookies()
            with open(cookies_file, 'w') as f:
                json.dump(cookies, f)
            logger.debug(f"Saved {len(cookies)} cookies to {cookies_file}")
        except Exception as e:
            logger.debug(f"Failed to save cookies: {e}")

    def _try_alternative_download_methods(self, page, doi_url: str, save_path: Path) -> bool:
        """Try alternative methods to download the MAIN PAPER PDF when standard methods fail.

        Validates downloaded files to ensure they are main papers, not supplements.
        """
        # Method 1: Transform URL patterns for main article PDF
        if "/doi/" in doi_url or "doi.org" in doi_url:
            pdf_url_patterns = [
                # Elsevier/ScienceDirect patterns - main article
                (r"sciencedirect.com/science/article/pii/([A-Z0-9]+)", r"sciencedirect.com/science/article/pii/\1/pdfft"),
                # Nature/Springer patterns - main article
                (r"nature.com/articles/([^/]+)$", r"nature.com/articles/\1.pdf"),
                (r"link.springer.com/article/([^/]+)$", r"link.springer.com/content/pdf/\1.pdf"),
                # Wiley patterns - main article
                (r"onlinelibrary.wiley.com/doi/([^/]+)/([^/]+)", r"onlinelibrary.wiley.com/doi/pdfdirect/\1/\2"),
                # Generic patterns
                (r"/abs/", r"/pdf/"),
                (r"/full/", r"/pdf/"),
                (r"/abstract/", r"/pdf/"),
            ]

            for pattern, replacement in pdf_url_patterns:
                try:
                    pdf_url = re.sub(pattern, replacement, doi_url)
                    if pdf_url != doi_url:
                        # Skip if URL looks like supplementary
                        if self._is_supplementary_link(pdf_url, ""):
                            continue

                        response = page.goto(pdf_url, timeout=30000)
                        if response and response.status == 200:
                            content_type = response.headers.get("content-type", "")
                            if "pdf" in content_type.lower():
                                try:
                                    with page.expect_download(timeout=30000) as download_info:
                                        page.goto(pdf_url)
                                    download = download_info.value
                                    filename = download.suggested_filename.lower()

                                    # Skip if filename suggests supplementary
                                    if self._is_supplementary_link(filename, ""):
                                        logger.debug(f"Skipping supplementary file: {filename}")
                                        download.delete()
                                        continue

                                    download.save_as(save_path)
                                    if self._validate_main_paper_pdf(save_path):
                                        logger.info(f"Downloaded main paper via URL transformation: {pdf_url[:60]}")
                                        return True
                                except:
                                    pass
                except:
                    continue

        # Method 2: Try direct append patterns for main article
        direct_patterns = [
            doi_url + ".pdf",
            doi_url + "/pdf",
            doi_url.rstrip('/') + "/pdf",
        ]

        for pdf_url in direct_patterns:
            # Skip if URL looks like supplementary
            if self._is_supplementary_link(pdf_url, ""):
                continue

            try:
                with page.expect_download(timeout=30000) as download_info:
                    page.goto(pdf_url)
                download = download_info.value
                filename = download.suggested_filename.lower()

                # Skip if filename suggests supplementary
                if self._is_supplementary_link(filename, ""):
                    logger.debug(f"Skipping supplementary file: {filename}")
                    download.delete()
                    continue

                download.save_as(save_path)
                if self._validate_main_paper_pdf(save_path):
                    logger.info(f"Downloaded main paper via direct pattern: {pdf_url[:60]}")
                    return True
            except:
                continue

        # Method 3: Look for iframe containing main PDF (not supplement)
        try:
            iframes = page.locator('iframe[src*=".pdf"], iframe[src*="/pdf"]').all()
            for iframe in iframes[:3]:
                if iframe.is_visible(timeout=2000):
                    iframe_src = iframe.get_attribute('src') or ""

                    # Skip supplementary iframes
                    if self._is_supplementary_link(iframe_src, ""):
                        continue

                    try:
                        with page.expect_download(timeout=30000) as download_info:
                            page.goto(iframe_src)
                        download = download_info.value
                        filename = download.suggested_filename.lower()

                        if self._is_supplementary_link(filename, ""):
                            download.delete()
                            continue

                        download.save_as(save_path)
                        if self._validate_main_paper_pdf(save_path):
                            logger.info(f"Downloaded main paper from iframe: {iframe_src[:60]}")
                            return True
                    except:
                        continue
        except:
            pass

        return False

    def _try_print_to_pdf_fallback(self, page, save_path: Path, article_url: str = None) -> bool:
        """
        Fallback method: Use Chrome DevTools Protocol to print page as PDF.

        This works when:
        - Direct PDF download is blocked by Cloudflare or other bot protection
        - Article content is viewable in browser but PDF download requires authentication
        - Publisher allows viewing but restricts downloads

        The method captures the rendered article page as a PDF using Chrome's
        built-in print-to-PDF functionality via CDP (Chrome DevTools Protocol).

        Args:
            page: Playwright page object
            save_path: Path to save the PDF
            article_url: Optional article URL to navigate to first (not PDF URL)

        Returns:
            True if PDF was successfully created, False otherwise
        """
        try:
            # If article_url provided, navigate to article page (not PDF URL)
            if article_url:
                logger.debug(f"Navigating to article page for print-to-PDF: {article_url}")
                page.goto(article_url, timeout=60000)
                page.wait_for_load_state("networkidle")

                # Wait for Cloudflare if present
                for _ in range(10):
                    title = page.title()
                    if "Just a moment" not in title and "Attention" not in title:
                        break
                    page.wait_for_timeout(2000)

            # Check if we have actual content (not Cloudflare or error page)
            title = page.title()
            if "Just a moment" in title or "Attention Required" in title:
                logger.debug("Page still showing Cloudflare challenge")
                return False

            if "404" in title or "Error" in title or "Access Denied" in title:
                logger.debug(f"Error page detected: {title}")
                return False

            # Wait for article content to fully render
            page.wait_for_timeout(3000)

            # Use CDP to print page as PDF
            logger.info("Attempting print-to-PDF fallback via CDP...")

            # Create CDP session
            cdp_session = page.context.new_cdp_session(page)

            # Print to PDF with optimal settings for articles
            pdf_data = cdp_session.send("Page.printToPDF", {
                "printBackground": True,
                "preferCSSPageSize": True,
                "scale": 0.9,
                "paperWidth": 8.5,
                "paperHeight": 11,
                "marginTop": 0.4,
                "marginBottom": 0.4,
                "marginLeft": 0.4,
                "marginRight": 0.4,
            })

            # Decode and save PDF
            pdf_bytes = base64.b64decode(pdf_data["data"])

            # Validate size (article PDFs should be substantial)
            if len(pdf_bytes) < 50000:
                logger.debug(f"Print-to-PDF output too small ({len(pdf_bytes)} bytes)")
                return False

            # Save the PDF
            save_path.parent.mkdir(parents=True, exist_ok=True)
            with open(save_path, 'wb') as f:
                f.write(pdf_bytes)

            # Verify it's a valid PDF
            with open(save_path, 'rb') as f:
                header = f.read(4)
                if header != b'%PDF':
                    logger.debug("Print-to-PDF output is not valid PDF")
                    save_path.unlink()
                    return False

            file_size = save_path.stat().st_size
            logger.info(f"Saved article via print-to-PDF ({file_size/1024:.1f} KB)")
            return True

        except Exception as e:
            logger.debug(f"Print-to-PDF fallback failed: {e}")
            return False

    def _validate_main_paper_pdf(self, save_path: Path) -> bool:
        """Validate that downloaded PDF is likely the main paper, not supplementary."""
        if not save_path.exists():
            return False

        size = save_path.stat().st_size

        # Main papers are typically > 100KB, supplements can be smaller
        if size < 50000:
            logger.debug(f"File too small ({size} bytes) - likely not main paper")
            save_path.unlink()
            return False

        # Check it's a valid PDF
        try:
            with open(save_path, 'rb') as f:
                header = f.read(4)
                if header != b'%PDF':
                    logger.debug("File is not a valid PDF")
                    save_path.unlink()
                    return False
        except:
            save_path.unlink()
            return False

        return True

    def download_paper(
        self,
        pub: Dict,
        subscription_only: bool = None,
        use_institutional_auth: bool = True
    ) -> Tuple[bool, Optional[str], str, bool]:
        """
        Download a paper and return (success, path, access_url, requires_subscription).

        Downloads papers using multiple strategies:
        1. PMC Open Access (direct HTTP)
        2. PMC via Playwright (if direct fails)
        3. DOI/Publisher (direct HTTP)
        4. DOI/Publisher via Playwright with institutional auth

        Args:
            pub: Publication dictionary with pmid, pmc_id, doi
            subscription_only: Override instance setting. If True, skip free papers
            use_institutional_auth: If True, attempt institutional SSO for subscription papers

        Returns:
            Tuple of (download_success, file_path, access_url, requires_subscription)
        """
        pmid = pub.get("pmid", "")
        pmc_id = pub.get("pmc_id")
        doi = pub.get("doi")

        paper_dir = self.get_paper_dir(pmid)
        pdf_path = paper_dir / "paper.pdf"

        # Determine if subscription-only mode
        use_subscription_only = subscription_only if subscription_only is not None else self.subscription_only

        # Check if paper has free access via PMC
        has_free = self.has_free_access(pmc_id)

        # If subscription_only mode and paper has free access, skip download
        if use_subscription_only and has_free:
            logger.info(f"Skipping PMID {pmid}: has free access (subscription-only mode)")
            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            return False, None, pubmed_url, False

        # Already downloaded?
        if pdf_path.exists() and pdf_path.stat().st_size > 1000:
            logger.info(f"PMID {pmid}: Already downloaded ({pdf_path.stat().st_size} bytes)")
            return True, str(pdf_path), str(pdf_path), not has_free

        # =============================================================
        # Strategy 1: PMC Open Access (free papers)
        # =============================================================
        if not use_subscription_only and pmc_id:
            pmc_url = self.check_pmc_access(pmc_id)
            if pmc_url:
                logger.info(f"PMID {pmid}: Trying PMC direct download...")

                # Try direct HTTP download first
                if self.download_pdf(pmc_url, pdf_path):
                    logger.info(f"PMID {pmid}: Downloaded via PMC direct")
                    return True, str(pdf_path), pmc_url, False

                # If direct fails, try PMC article page via Playwright
                pmc_article_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/"
                logger.info(f"PMID {pmid}: PMC direct failed, trying Playwright...")
                if self.download_via_playwright(pmc_article_url, pdf_path, use_institutional_auth=False):
                    logger.info(f"PMID {pmid}: Downloaded via PMC Playwright")
                    return True, str(pdf_path), pmc_article_url, False

        # =============================================================
        # Strategy 2: DOI/Publisher (may require institutional access)
        # =============================================================
        if doi:
            doi_url = self.get_doi_url(doi)
            if doi_url:
                publisher = self._detect_publisher(doi_url)
                logger.info(f"PMID {pmid}: Trying DOI from {publisher or 'unknown publisher'}...")

                # Try standard HTTP download first
                try:
                    response = self.session.get(doi_url, timeout=15)
                    soup = BeautifulSoup(response.text, 'html.parser')

                    # Look for direct PDF link
                    pdf_links = soup.find_all('a', href=re.compile(r'\.pdf', re.I))
                    for link in pdf_links[:3]:
                        pdf_url = link.get('href')
                        if pdf_url:
                            if not pdf_url.startswith('http'):
                                pdf_url = urljoin(doi_url, pdf_url)

                            if self.download_pdf(pdf_url, pdf_path):
                                logger.info(f"PMID {pmid}: Downloaded via DOI direct")
                                return True, str(pdf_path), doi_url, not has_free

                except Exception as e:
                    logger.debug(f"PMID {pmid}: Publisher page fetch failed: {e}")

                # =============================================================
                # Strategy 3: Playwright with institutional authentication
                # =============================================================
                # For papers without free access OR in subscription-only mode
                if use_institutional_auth and (not has_free or use_subscription_only):
                    logger.info(f"PMID {pmid}: Attempting institutional access download...")

                    # Check if we have credentials
                    if self.tak_account and self.tak_key:
                        if self.download_via_playwright(doi_url, pdf_path, use_institutional_auth=True):
                            logger.info(f"PMID {pmid}: Downloaded via institutional access")
                            return True, str(pdf_path), doi_url, True
                        else:
                            logger.warning(f"PMID {pmid}: Institutional access download failed")
                    else:
                        logger.warning(f"PMID {pmid}: No credentials for institutional access")
                        # Still try Playwright without auth (might work for some papers)
                        if self.download_via_playwright(doi_url, pdf_path, use_institutional_auth=False):
                            logger.info(f"PMID {pmid}: Downloaded via Playwright (no auth)")
                            return True, str(pdf_path), doi_url, not has_free

                return False, None, doi_url, not has_free

        # =============================================================
        # Fallback: Return PubMed URL (download failed)
        # =============================================================
        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        logger.warning(f"PMID {pmid}: Download failed - no PMC or DOI available")
        return False, None, pubmed_url, not has_free

    # =============================================================
    # Web-Based Full Text Extraction (without PDF download)
    # =============================================================

    def fetch_article_html(
        self,
        pmc_id: Optional[str],
        doi: Optional[str],
        use_institutional_auth: bool = True
    ) -> Optional[str]:
        """
        Fetch full article HTML content via Playwright without downloading PDF.

        Priority:
        1. PMC HTML (open access, consistent structure)
        2. Publisher DOI page (may require institutional auth)

        Args:
            pmc_id: PubMed Central ID (e.g., "PMC1234567")
            doi: Digital Object Identifier
            use_institutional_auth: Whether to attempt institutional login

        Returns:
            Plain text extracted from HTML, or None if access fails.
        """
        try:
            from playwright.sync_api import sync_playwright
        except ImportError:
            logger.warning("Playwright not installed for web scraping")
            return None

        # Strategy 1: Try PMC first (open access, consistent HTML)
        if pmc_id:
            pmc_clean = pmc_id.replace("PMC", "") if pmc_id.startswith("PMC") else pmc_id
            pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_clean}/"
            logger.debug(f"Trying PMC HTML: {pmc_url}")
            text = self._fetch_html_from_url(pmc_url, use_auth=False)
            if text and len(text) > 1000:
                logger.info(f"Fetched article from PMC ({len(text)} chars)")
                return text

        # Strategy 2: Try publisher DOI (may need auth)
        if doi:
            doi_url = self.get_doi_url(doi)
            if doi_url:
                logger.debug(f"Trying publisher HTML: {doi_url}")
                text = self._fetch_html_from_url(doi_url, use_auth=use_institutional_auth)
                if text and len(text) > 1000:
                    logger.info(f"Fetched article from publisher ({len(text)} chars)")
                    return text

        return None

    def _fetch_html_from_url(self, url: str, use_auth: bool = True) -> Optional[str]:
        """
        Fetch and parse HTML from a single URL using Playwright.

        Args:
            url: URL to fetch
            use_auth: Whether to attempt institutional authentication

        Returns:
            Extracted text content, or None if fetch fails.
        """
        try:
            from playwright.sync_api import sync_playwright

            with sync_playwright() as p:
                browser = p.chromium.launch(headless=True)
                context = browser.new_context(
                    user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36",
                    viewport={"width": 1920, "height": 1080},
                    java_script_enabled=True,
                    locale="en-US",
                    timezone_id="America/New_York",
                )

                # Load cookies if available
                context_dir = self._get_persistent_context_dir()
                cookies_file = context_dir / "cookies.json"
                if cookies_file.exists():
                    try:
                        with open(cookies_file, 'r') as f:
                            cookies = json.load(f)
                        context.add_cookies(cookies)
                    except Exception as e:
                        logger.debug(f"Failed to load cookies: {e}")

                page = context.new_page()

                page.goto(url, timeout=60000)
                page.wait_for_load_state("networkidle")

                # Check if we need institutional auth
                if use_auth and self._needs_authentication(page):
                    publisher = self._detect_publisher(url)
                    if publisher and self.tak_account and self.tak_key:
                        logger.info(f"Page requires authentication, attempting institutional access...")
                        self._navigate_to_institutional_access(page, publisher)
                        self._handle_sso_login(page, publisher)
                        page.goto(url, timeout=60000)
                        page.wait_for_load_state("networkidle")

                html_content = page.content()
                self._save_cookies(context, cookies_file)
                browser.close()

                return self._extract_text_from_html(html_content, url)

        except Exception as e:
            logger.debug(f"Failed to fetch HTML from {url}: {e}")
            return None

    def _needs_authentication(self, page) -> bool:
        """
        Check if page requires authentication to view full content.

        Looks for common paywall/login indicators on publisher pages.
        """
        paywall_indicators = [
            'text="Access through your institution"',
            'text="Sign in"',
            'text="Subscribe"',
            'text="Purchase"',
            'text="Get access"',
            'text="Log in"',
            '.paywall',
            '.login-required',
            '[data-access-type="restricted"]',
            '.access-options',
            '#access-denial',
        ]

        for indicator in paywall_indicators:
            try:
                if page.locator(indicator).count() > 0:
                    logger.debug(f"Found paywall indicator: {indicator}")
                    return True
            except:
                pass

        return False

    def _extract_text_from_html(self, html: str, url: str) -> Optional[str]:
        """
        Extract article text from HTML using publisher-specific selectors.

        Args:
            html: Raw HTML content
            url: Source URL (used to detect publisher)

        Returns:
            Extracted plain text, or None if extraction fails.
        """
        soup = BeautifulSoup(html, 'html.parser')
        publisher = self._detect_publisher(url)

        # Publisher-specific article content selectors
        ARTICLE_SELECTORS = {
            "springer": ["article", ".c-article-body", ".article-body", "#body"],
            "elsevier": ["article", ".article-content", "#body", ".Body", ".article__body"],
            "wiley": [".article-section__content", ".article__body", "article"],
            "oxford": [".article-body", ".widget-ArticleFulltext", "article"],
            "taylor_francis": [".article-body", "#article-content", "article"],
            "acs": [".article_content", ".NLM_sec_level_1", "article"],
            None: ["article", "main", ".content", "#content", "#article-body"],  # Generic/PMC fallback
        }

        selectors = ARTICLE_SELECTORS.get(publisher, ARTICLE_SELECTORS[None])

        for selector in selectors:
            element = soup.select_one(selector)
            if element:
                # Remove script/style/nav elements
                for tag in element.find_all(['script', 'style', 'nav', 'header', 'footer', 'aside']):
                    tag.decompose()

                text = element.get_text(separator='\n', strip=True)
                if len(text) > 500:  # Minimum content threshold
                    logger.debug(f"Extracted {len(text)} chars using selector: {selector}")
                    return text

        # Fallback: extract all body text
        body = soup.find('body')
        if body:
            for tag in body.find_all(['script', 'style', 'nav', 'header', 'footer', 'aside']):
                tag.decompose()
            return body.get_text(separator='\n', strip=True)

        return None
