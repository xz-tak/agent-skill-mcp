#!/usr/bin/env python3
"""
auth_setup.py

One-time Okta login using Playwright (headless) and save authenticated
storage state to a JSON file (okta_auth_state.json).

Usage
=====
Set environment variables:

    export OKTA_EMAIL="your.email@company.com"
    export OKTA_PASSWORD="your-okta-password"
    # Optional overrides:
    # export OKTA_ENTRY_URL="https://yourcompany.okta.com"
    # export APP_LANDING_URL_PATTERN="**/your/app/landing*"

Then run:

    python auth_setup.py

This will:
  - open Okta in a headless browser
  - fill email + password
  - click through MFA factor selection (e.g. Okta Verify / Push) if present
  - wait for you to APPROVE the push on your phone
  - once redirected into the target app, save cookies/localStorage to
    okta_auth_state.json

Later, in your real automation scripts, you can do:

    context = await browser.new_context(storage_state="okta_auth_state.json")

to reuse the logged-in session without repeating the login flow.
"""

import os
import re
import asyncio
from playwright.async_api import async_playwright, TimeoutError as PlaywrightTimeoutError


# ----------------------------------------------------------------------
# Configuration – override via environment if needed
# ----------------------------------------------------------------------

# Okta entry URL –  this can be the Okta org URL or an app-specific link
OKTA_ENTRY_URL = os.environ.get("OKTA_ENTRY_URL", "https://takeda.okta.com")

# Pattern of the URL once you are successfully inside the target app.
# Use a glob-style pattern accepted by Playwright's wait_for_url, e.g.:
#   "**/app/yourapp/**" or "**/dashboard*"
APP_LANDING_URL_PATTERN = os.environ.get(
    "APP_LANDING_URL_PATTERN",
    "**/app/**",
)

# Where to save the storage state
AUTH_STATE_FILE = os.environ.get("OKTA_AUTH_STATE_FILE", "okta_auth_state.json")

# Credentials – MUST be set in env; fail loud if missing
OKTA_EMAIL = os.environ.get("OKTA_EMAIL")
OKTA_PASSWORD = os.environ.get("OKTA_PASSWORD")

if not OKTA_EMAIL or not OKTA_PASSWORD:
    raise SystemExit(
        "ERROR: Please set OKTA_EMAIL and OKTA_PASSWORD environment variables before running.\n"
        "Example:\n"
        "  export OKTA_EMAIL='your.email@company.com'\n"
        "  export OKTA_PASSWORD='your-okta-password'\n"
    )


# ----------------------------------------------------------------------
# Helper: login flow
# ----------------------------------------------------------------------
async def login_okta_and_save_state() -> None:
    async with async_playwright() as p:
        browser = await p.chromium.launch(
            headless=True,
            args=["--no-sandbox", "--disable-dev-shm-usage"],
        )
        context = await browser.new_context()
        page = await context.new_page()

        print(f"🔐 Going to Okta entry URL: {OKTA_ENTRY_URL}")
        await page.goto(OKTA_ENTRY_URL)

        # --------------------------------------------------------------
        # STEP 1 – email / username
        # --------------------------------------------------------------
        email_selector = (
            'input[type="email"], '
            'input[name="identifier" i], '
            'input[name="username" i], '
            '#okta-signin-username'
        )
        print("✉️  Waiting for email/username field...")
        await page.wait_for_selector(email_selector, timeout=60_000)
        await page.locator(email_selector).first.fill(OKTA_EMAIL)

        # Click "Next" / "Sign in"
        try:
            await page.click("#idp-discovery-submit", timeout=5_000)
        except PlaywrightTimeoutError:
            print("⬇️  Clicking email submit button by role/name...")
            await page.get_by_role(
                "button",
                name=re.compile(r"(next|sign in|continue)", re.I),
            ).click(timeout=30_000)

        # --------------------------------------------------------------
        # STEP 2 – first password page
        # --------------------------------------------------------------
        print("🔑 Waiting for password field...")
        await page.wait_for_selector(
            'input[type="password"], #okta-signin-password', timeout=60_000
        )
        await page.fill(
            'input[type="password"], #okta-signin-password', OKTA_PASSWORD
        )

        # Click "Sign in" / "Verify"
        try:
            await page.click("#okta-signin-submit", timeout=5_000)
        except PlaywrightTimeoutError:
            print("⬇️  Clicking password submit button by role/name...")
            await page.get_by_role(
                "button",
                name=re.compile(r"(verify|sign in|continue)", re.I),
            ).click(timeout=30_000)

        # --------------------------------------------------------------
        # STEP 3 – optional "Verify with your password" page
        # --------------------------------------------------------------
        try:
            print("🛡  Checking for 'Verify with your password' page…")
            await page.get_by_text("Verify with your password").wait_for(timeout=30_000)

            shot = "okta_verify_with_password.png"
            await page.screenshot(path=shot, full_page=True)
            print(f"📷 Screenshot saved: {shot}")

            # Ensure password is filled on this screen too
            try:
                await page.fill('input[type="password"]', OKTA_PASSWORD)
            except Exception:
                pass

            # Click the blue Verify button
            print("🔒 Clicking blue 'Verify' button…")
            await page.get_by_role(
                "button",
                name=re.compile(r"^verify$", re.I),
            ).click(timeout=30_000)

        except PlaywrightTimeoutError:
            print("ℹ️  Second 'Verify with your password' page not shown; continuing.")

        # --------------------------------------------------------------
        # STEP 4 – MFA options: choose "Get a push notification"
        # --------------------------------------------------------------
        try:
            print("🛡  Waiting for MFA options screen…")
            await page.get_by_text(
                "Verify it's you with a security method"
            ).wait_for(timeout=60_000)

            mfa_shot = "okta_mfa_options.png"
            await page.screenshot(path=mfa_shot, full_page=True)
            print(f"📷 MFA options screenshot: {mfa_shot}")

            # First try 'button' role
            buttons = page.get_by_role("button", name="Select")
            count = await buttons.count()
            print(f"   Found {count} 'Select' buttons with role=button")

            if count == 0:
                # Some themes use <a> as links instead of buttons
                links = page.get_by_role("link", name="Select")
                count_links = await links.count()
                print(f"   Found {count_links} 'Select' links with role=link")
                buttons = links
                count = count_links

            if count == 0:
                raise PlaywrightTimeoutError(
                    "No 'Select' control found on MFA options screen."
                )

            # 0: "Enter a code", 1: "Get a push notification"
            if count >= 2:
                await buttons.nth(1).click(timeout=30_000)
                print("📲 Clicked second 'Select' (Get a push notification).")
            else:
                await buttons.first().click(timeout=30_000)
                print("📲 Only one 'Select' found; clicked it (fallback).")

            # 🔍 Wait a moment for the code screen to render, then screenshot it
            await page.wait_for_timeout(3_000)  # 3 seconds; adjust if needed
            push_code_shot = "okta_push_code.png"
            await page.screenshot(path=push_code_shot, full_page=True)
            print(f"📷 Push code screen screenshot: {push_code_shot}")

        except PlaywrightTimeoutError as e:
            mfa_timeout_shot = "okta_mfa_step4_timeout.png"
            await page.screenshot(path=mfa_timeout_shot, full_page=True)
            print(f"⚠️  Could not click MFA 'Select' button: {e}")
            print(f"📷 Screenshot: {mfa_timeout_shot}")
            # We still continue: maybe Okta accepted previous step as enough.
            
        # --------------------------------------------------------------
        # STEP 5 – wait until we’re in Okta dashboard / app,
        #          then screenshot the page *after* push verification
        # --------------------------------------------------------------
        print("\n⏳ Waiting for Okta dashboard / app page…")
        try:
            # ✅ At this point Okta has accepted the push and redirected.
            post_login_shot = "okta_after_push_verified_0.png"
            await page.screenshot(path=post_login_shot, full_page=True)
            print(f"📷 Post-login screenshot saved: {post_login_shot}")
            
            # While this is waiting, you approve the push on your phone.
            await page.wait_for_url(APP_LANDING_URL_PATTERN, timeout=180_000)

            # ✅ At this point Okta has accepted the push and redirected.
            post_login_shot = "okta_after_push_verified.png"
            await page.screenshot(path=post_login_shot, full_page=True)
            print(f"📷 Post-login screenshot saved: {post_login_shot}")

        except PlaywrightTimeoutError:
            current_url = page.url
            post_shot = "okta_post_login_timeout.png"
            await page.screenshot(path=post_shot, full_page=True)
            print(f"❌ Timed out waiting for landing page. Current URL: {current_url}")
            print(f"📷 Screenshot (timeout state): {post_shot}")
            await browser.close()
            raise SystemExit(
                "Could not detect successful app login. "
                "You may need to tweak APP_LANDING_URL_PATTERN."
            )

        # --------------------------------------------------------------
        # STEP 6 – save storage state
        # --------------------------------------------------------------
        print(f"✅ Logged in. Saving storage state to {AUTH_STATE_FILE} …")
        await context.storage_state(path=AUTH_STATE_FILE)
        print("✅ Done. You can reuse this state in other Playwright scripts.")
        await browser.close()


if __name__ == "__main__":
    asyncio.run(login_okta_and_save_state())
