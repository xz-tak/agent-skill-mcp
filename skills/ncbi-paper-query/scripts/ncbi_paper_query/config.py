"""
NCBI Paper Query - Configuration Module

Global configuration, logging setup, and NCBI API settings.
"""

import os
import logging
import shutil
from pathlib import Path

from dotenv import load_dotenv
from Bio import Entrez

# Load environment variables
load_dotenv()

# Skill directory - defaults to ~/.claude/skills/ncbi-paper-query
# Can be overridden by SKILL_DIR environment variable
SKILL_DIR = Path(os.environ.get("SKILL_DIR", Path.home() / ".claude" / "skills" / "ncbi-paper-query"))

# Check if claude CLI is available
CLAUDE_CLI_AVAILABLE = shutil.which("claude") is not None

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# NCBI API Configuration
NCBI_EMAIL = os.getenv("TAK_ACCOUNT", "user@example.com")
NCBI_API_KEY = os.getenv("NCBI_API_KEY", None)
Entrez.email = NCBI_EMAIL
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

# PDF parsing availability
try:
    from PyPDF2 import PdfReader
    PDF_READER_AVAILABLE = True
except ImportError:
    PdfReader = None
    PDF_READER_AVAILABLE = False
