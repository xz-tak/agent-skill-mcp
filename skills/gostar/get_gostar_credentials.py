#!/usr/bin/env python3
"""
Retrieve GOSTAR database credentials from AWS Secrets Manager or environment.
"""

import os
import json
import boto3
from botocore.exceptions import ClientError


def get_gostar_password():
    """
    Get GOSTAR password from environment or AWS Secrets Manager.

    Priority:
    1. Environment variable GOSTAR_RO_PASSWORD (for GitHub Actions)
    2. AWS Secrets Manager (for EC2/local development)

    Returns:
        str: The GOSTAR read-only password

    Raises:
        ValueError: If password cannot be retrieved
    """
    # First try environment variable (GitHub Actions)
    password = os.environ.get('GOSTAR_RO_PASSWORD')
    if password:
        return password

    # Try AWS Secrets Manager (EC2/local)
    secret_name = "gostar/ro-password"
    region_name = "us-east-1"

    try:
        session = boto3.session.Session()
        client = session.client(
            service_name='secretsmanager',
            region_name=region_name
        )

        get_secret_value_response = client.get_secret_value(
            SecretId=secret_name
        )

        # Secrets Manager returns string or binary
        if 'SecretString' in get_secret_value_response:
            secret = get_secret_value_response['SecretString']
            # Try to parse as JSON first
            try:
                secret_dict = json.loads(secret)
                return secret_dict.get('password', secret)
            except json.JSONDecodeError:
                # It's a plain string
                return secret
        else:
            # Binary secret
            return get_secret_value_response['SecretBinary'].decode('utf-8')

    except ClientError as e:
        error_code = e.response['Error']['Code']
        if error_code == 'ResourceNotFoundException':
            raise ValueError(f"Secret '{secret_name}' not found in AWS Secrets Manager")
        elif error_code == 'AccessDeniedException':
            raise ValueError(f"No permission to access secret '{secret_name}'")
        else:
            raise ValueError(f"Error retrieving secret: {e}")

    raise ValueError("GOSTAR_RO_PASSWORD not found in environment or AWS Secrets Manager")


def get_gostar_connection_string():
    """
    Get PostgreSQL connection string for GOSTAR database.

    Returns:
        str: PostgreSQL connection string
    """
    password = get_gostar_password()
    return f"postgresql://gostar_ro:{password}@usvgarps11158-dev003.cm9aqaugy64i.us-east-1.rds.amazonaws.com:5442/gostar"


def get_gostar_connection_params():
    """
    Get connection parameters for psycopg2.

    Returns:
        dict: Connection parameters
    """
    return {
        'host': 'usvgarps11158-dev003.cm9aqaugy64i.us-east-1.rds.amazonaws.com',
        'port': 5442,
        'database': 'gostar',
        'user': 'gostar_ro',
        'password': get_gostar_password(),
        'connect_timeout': 10
    }


if __name__ == "__main__":
    # Test credential retrieval
    try:
        password = get_gostar_password()
        print(f"✓ Successfully retrieved password (length: {len(password)})")
        print(f"✓ Connection string: postgresql://gostar_ro:***@usvgarps11158-dev003.cm9aqaugy64i.us-east-1.rds.amazonaws.com:5442/gostar")
    except ValueError as e:
        print(f"✗ Failed to retrieve password: {e}")
        exit(1)
