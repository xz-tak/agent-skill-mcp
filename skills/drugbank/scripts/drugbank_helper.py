"""
DrugBank Database Helper

Provides utilities for querying the DrugBank PostgreSQL database.
Credentials are automatically loaded from AWS Secrets Manager:
- Secret name: DRUGBANK_RO_PASSWORD
- Profiles: cmp-dev, sci-dev

Falls back to environment variables if AWS Secrets Manager is not available:
- DRUGBANK_USERNAME
- DRUGBANK_PASSWORD
"""

import os
import sys
import json
from typing import Optional, Dict, List, Any
import psycopg2
from psycopg2.extras import RealDictCursor


# Database connection constants
DB_HOST = "usvgarps11158-dev003.cm9aqaugy64i.us-east-1.rds.amazonaws.com"
DB_PORT = 5442
DB_NAME = "drugbank"


def get_credentials_from_secrets_manager(
    secret_name: str = "DRUGBANK_RO_PASSWORD",
    profile: Optional[str] = None
) -> Optional[tuple[str, str]]:
    """
    Get DrugBank credentials from AWS Secrets Manager.

    Args:
        secret_name: Name of the secret in AWS Secrets Manager
        profile: AWS profile to use (tries cmp-dev and sci-dev by default)

    Returns:
        Tuple of (username, password) or None if not found
    """
    try:
        import boto3
        from botocore.exceptions import ClientError, NoCredentialsError

        profiles = [profile] if profile else ["cmp-dev", "sci-dev"]

        for prof in profiles:
            try:
                session = boto3.Session(profile_name=prof)
                client = session.client("secretsmanager")

                response = client.get_secret_value(SecretId=secret_name)

                if "SecretString" in response:
                    secret = json.loads(response["SecretString"])
                    username = secret.get("username")
                    password = secret.get("password")

                    if username and password:
                        return username, password

            except (ClientError, NoCredentialsError):
                continue

    except ImportError:
        pass

    return None


def get_credentials() -> tuple[str, str]:
    """
    Get DrugBank credentials from AWS Secrets Manager or environment variables.

    Tries in order:
    1. AWS Secrets Manager (profiles: cmp-dev, sci-dev)
    2. Environment variables (DRUGBANK_USERNAME, DRUGBANK_PASSWORD)

    Returns:
        Tuple of (username, password)

    Raises:
        ValueError: If credentials are not found
    """
    # Try AWS Secrets Manager first
    credentials = get_credentials_from_secrets_manager()
    if credentials:
        return credentials

    # Fall back to environment variables
    username = os.getenv("DRUGBANK_USERNAME")
    password = os.getenv("DRUGBANK_PASSWORD")

    if username and password:
        return username, password

    raise ValueError(
        "DrugBank credentials not found. Please configure credentials using one of:\n\n"
        "1. AWS Secrets Manager:\n"
        "   - Secret name: DRUGBANK_RO_PASSWORD\n"
        "   - Profiles: cmp-dev or sci-dev\n"
        "   - Format: {\"username\": \"your_username\", \"password\": \"your_password\"}\n\n"
        "2. Environment variables:\n"
        "   export DRUGBANK_USERNAME=your_username\n"
        "   export DRUGBANK_PASSWORD=your_password\n"
    )


def get_connection():
    """
    Create and return a database connection.

    Returns:
        psycopg2 connection object

    Raises:
        ValueError: If credentials are not set
        psycopg2.Error: If connection fails
    """
    username, password = get_credentials()

    return psycopg2.connect(
        host=DB_HOST,
        port=DB_PORT,
        database=DB_NAME,
        user=username,
        password=password
    )


def execute_query(query: str, params: Optional[tuple] = None) -> List[Dict[str, Any]]:
    """
    Execute a SQL query and return results as a list of dictionaries.

    Args:
        query: SQL query string
        params: Optional tuple of query parameters

    Returns:
        List of dictionaries with column names as keys

    Raises:
        ValueError: If credentials are not set
        psycopg2.Error: If query fails
    """
    conn = get_connection()
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute(query, params or ())
            return [dict(row) for row in cursor.fetchall()]
    finally:
        conn.close()


# ============================================================================
# Pre-built Query Functions
# ============================================================================

def search_drugs(search_term: str, limit: int = 50) -> List[Dict[str, Any]]:
    """
    Search drugs by name or DrugBank ID.

    Args:
        search_term: Drug name or ID to search for
        limit: Maximum number of results (default: 50)

    Returns:
        List of drug dictionaries with drugbank_id, name, type, state, description
    """
    query = """
        SELECT d.drugbank_id, d.name, d.type, d.state, d.description
        FROM drugs d
        WHERE d.name ILIKE %s OR d.drugbank_id ILIKE %s
        LIMIT %s
    """
    search_pattern = f"%{search_term}%"
    return execute_query(query, (search_pattern, search_pattern, limit))


def get_drug_by_id(drugbank_id: str) -> Optional[Dict[str, Any]]:
    """
    Get detailed drug information by DrugBank ID.

    Args:
        drugbank_id: DrugBank ID (e.g., "DB00112")

    Returns:
        Dictionary with drug details including molecular properties, or None if not found
    """
    query = """
        SELECT d.*,
               dcp.smiles,
               dcp.molecular_formula,
               dcp.molecular_weight,
               dcp.inchi,
               dcp.inchikey
        FROM drugs d
        LEFT JOIN drug_calculated_properties dcp ON d.id = dcp.drug_id
        WHERE d.drugbank_id = %s
    """
    results = execute_query(query, (drugbank_id,))
    return results[0] if results else None


def get_drug_targets(drug_identifier: str) -> List[Dict[str, Any]]:
    """
    Get all targets for a drug (by name or DrugBank ID).

    Args:
        drug_identifier: Drug name or DrugBank ID

    Returns:
        List of target dictionaries with target name, gene, mechanism info
    """
    query = """
        SELECT DISTINCT
            d.drugbank_id,
            d.name AS drug_name,
            be.name AS target_name,
            p.gene_name,
            p.general_function,
            p.specific_function,
            b.inhibitor,
            b.agonist,
            b.antagonist,
            b.pharmacological_action
        FROM drugs d
        JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
        JOIN bio_entities be ON b.biodb_id = be.biodb_id
        LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
        LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
        WHERE d.name ILIKE %s OR d.drugbank_id ILIKE %s
    """
    search_pattern = f"%{drug_identifier}%"
    return execute_query(query, (search_pattern, search_pattern))


def get_clinical_trials(drug_identifier: str, limit: int = 100) -> List[Dict[str, Any]]:
    """
    Get clinical trials for a drug.

    Args:
        drug_identifier: Drug name or DrugBank ID
        limit: Maximum number of results (default: 100)

    Returns:
        List of clinical trial dictionaries
    """
    query = """
        SELECT DISTINCT
            ct.identifier,
            ct.title,
            ct.status,
            ct.purpose,
            ct.start_date,
            ct.end_date,
            ct.phase,
            d.drugbank_id,
            d.name AS drug_name
        FROM clinical_trials ct
        JOIN clinical_trial_interventions cti ON ct.identifier = cti.trial_id
        JOIN clinical_trial_interventions_drugs ctid ON cti.id = ctid.intervention_id
        JOIN drugs d ON ctid.drug_id = d.id
        WHERE d.name ILIKE %s OR d.drugbank_id ILIKE %s
        ORDER BY ct.start_date DESC NULLS LAST
        LIMIT %s
    """
    search_pattern = f"%{drug_identifier}%"
    return execute_query(query, (search_pattern, search_pattern, limit))


def get_drug_indications(drug_identifier: str) -> List[Dict[str, Any]]:
    """
    Get approved and off-label indications for a drug.

    Args:
        drug_identifier: Drug name or DrugBank ID

    Returns:
        List of indication dictionaries
    """
    query = """
        SELECT DISTINCT
            d.drugbank_id,
            d.name AS drug_name,
            c.title AS indication,
            si.kind,
            si.off_label,
            si.country,
            c.snomed_id,
            c.meddra_id,
            c.icd10_id
        FROM structured_indications si
        JOIN drugs d ON si.drug_id = d.id
        JOIN indication_conditions ic ON si.id = ic.indication_id
            AND ic.relationship = 'for_condition'
        JOIN conditions c ON ic.condition_id = c.id
        WHERE d.name ILIKE %s OR d.drugbank_id ILIKE %s
        ORDER BY si.off_label, c.title
    """
    search_pattern = f"%{drug_identifier}%"
    return execute_query(query, (search_pattern, search_pattern))


def get_adverse_effects(
    drug_identifier: str,
    min_frequency: Optional[float] = None
) -> List[Dict[str, Any]]:
    """
    Get adverse effects for a drug.

    Args:
        drug_identifier: Drug name or DrugBank ID
        min_frequency: Optional minimum frequency threshold (e.g., 0.01 for 1%)

    Returns:
        List of adverse effect dictionaries with frequency data
    """
    query = """
        SELECT
            d.drugbank_id,
            d.name AS drug_name,
            c.title AS adverse_effect,
            aei.percent AS frequency_percent,
            sae.severity,
            c.snomed_id,
            c.meddra_id
        FROM structured_adverse_effects sae
        JOIN drugs d ON sae.drug_id = d.id
        JOIN adverse_effect_conditions aec ON sae.id = aec.adverse_effect_id
            AND aec.relationship = 'effect'
        JOIN conditions c ON aec.condition_id = c.id
        LEFT JOIN adverse_effect_incidences aei ON sae.id = aei.adverse_effect_id
        WHERE (d.name ILIKE %s OR d.drugbank_id ILIKE %s)
    """

    search_pattern = f"%{drug_identifier}%"
    params = [search_pattern, search_pattern]

    if min_frequency is not None:
        query += " AND aei.percent >= %s"
        params.append(min_frequency)

    query += " ORDER BY aei.percent DESC NULLS LAST"

    return execute_query(query, tuple(params))


def search_by_target(target_name: str, limit: int = 50) -> List[Dict[str, Any]]:
    """
    Find drugs targeting a specific protein or gene.

    Args:
        target_name: Protein name or gene name
        limit: Maximum number of results (default: 50)

    Returns:
        List of drug dictionaries with target information
    """
    query = """
        SELECT DISTINCT
            d.drugbank_id,
            d.name AS drug_name,
            d.type,
            d.state,
            be.name AS target_name,
            p.gene_name,
            b.inhibitor,
            b.agonist,
            b.antagonist,
            b.pharmacological_action
        FROM drugs d
        JOIN bonds b ON d.id = b.drug_id AND b.type = 'TargetBond'
        JOIN bio_entities be ON b.biodb_id = be.biodb_id
        LEFT JOIN bio_entity_components bec ON be.biodb_id = bec.biodb_id
        LEFT JOIN polypeptides p ON bec.component_id = p.uniprot_id
        WHERE be.name ILIKE %s OR p.gene_name ILIKE %s
        ORDER BY d.state, d.name
        LIMIT %s
    """
    search_pattern = f"%{target_name}%"
    return execute_query(query, (search_pattern, search_pattern, limit))


def get_drug_interactions(drug_identifier: str) -> List[Dict[str, Any]]:
    """
    Get drug-drug interactions for a drug.

    Args:
        drug_identifier: Drug name or DrugBank ID

    Returns:
        List of drug interaction dictionaries
    """
    query = """
        SELECT DISTINCT
            d1.drugbank_id AS drug1_id,
            d1.name AS drug1_name,
            d2.drugbank_id AS drug2_id,
            d2.name AS drug2_name,
            sdi.description
        FROM structured_drug_interactions sdi
        JOIN drugs d1 ON sdi.drug_id = d1.id
        JOIN drugs d2 ON sdi.interacting_drug_id = d2.id
        WHERE d1.name ILIKE %s OR d1.drugbank_id ILIKE %s
        ORDER BY d2.name
    """
    search_pattern = f"%{drug_identifier}%"
    return execute_query(query, (search_pattern, search_pattern))


# ============================================================================
# CLI Interface
# ============================================================================

def main():
    """Command-line interface for testing queries."""
    if len(sys.argv) < 2:
        print("Usage: python drugbank_helper.py <function> [args...]")
        print("\nAvailable functions:")
        print("  search_drugs <term>")
        print("  get_drug_by_id <drugbank_id>")
        print("  get_drug_targets <drug>")
        print("  get_clinical_trials <drug>")
        print("  get_drug_indications <drug>")
        print("  get_adverse_effects <drug> [min_frequency]")
        print("  search_by_target <target>")
        print("  get_drug_interactions <drug>")
        sys.exit(1)

    func_name = sys.argv[1]
    args = sys.argv[2:]

    try:
        if func_name == "search_drugs":
            results = search_drugs(args[0])
        elif func_name == "get_drug_by_id":
            results = get_drug_by_id(args[0])
        elif func_name == "get_drug_targets":
            results = get_drug_targets(args[0])
        elif func_name == "get_clinical_trials":
            results = get_clinical_trials(args[0])
        elif func_name == "get_drug_indications":
            results = get_drug_indications(args[0])
        elif func_name == "get_adverse_effects":
            min_freq = float(args[1]) if len(args) > 1 else None
            results = get_adverse_effects(args[0], min_freq)
        elif func_name == "search_by_target":
            results = search_by_target(args[0])
        elif func_name == "get_drug_interactions":
            results = get_drug_interactions(args[0])
        else:
            print(f"Unknown function: {func_name}")
            sys.exit(1)

        # Pretty print results
        import json
        print(json.dumps(results, indent=2, default=str))

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
