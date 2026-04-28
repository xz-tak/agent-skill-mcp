#!/usr/bin/env python3
"""
Setup DrugBank credentials in AWS Secrets Manager.

This script stores DrugBank database credentials in AWS Secrets Manager
for the specified profiles.
"""

import sys
import json
import argparse
import getpass
import boto3
from botocore.exceptions import ClientError


def create_or_update_secret(
    secret_name: str,
    username: str,
    password: str,
    profile: str,
    region: str = "us-east-1"
) -> bool:
    """
    Create or update a secret in AWS Secrets Manager.

    Args:
        secret_name: Name of the secret
        username: Database username
        password: Database password
        profile: AWS profile to use
        region: AWS region

    Returns:
        True if successful, False otherwise
    """
    try:
        session = boto3.Session(profile_name=profile, region_name=region)
        client = session.client("secretsmanager")

        secret_value = json.dumps({
            "username": username,
            "password": password
        })

        try:
            # Try to create the secret
            response = client.create_secret(
                Name=secret_name,
                Description="DrugBank database credentials",
                SecretString=secret_value
            )
            print(f"✓ Created DrugBank credentials secret in profile '{profile}'")
            print(f"  ARN: {response['ARN']}")
            return True

        except ClientError as e:
            if e.response["Error"]["Code"] == "ResourceExistsException":
                # Secret exists, update it
                response = client.update_secret(
                    SecretId=secret_name,
                    SecretString=secret_value
                )
                print(f"✓ Updated DrugBank credentials secret in profile '{profile}'")
                print(f"  ARN: {response['ARN']}")
                return True
            else:
                raise

    except ClientError as e:
        print(f"✗ Error with profile '{profile}': {e.response['Error']['Message']}")
        return False
    except Exception as e:
        print(f"✗ Error with profile '{profile}': {str(e)}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Store DrugBank credentials in AWS Secrets Manager"
    )
    parser.add_argument(
        "--secret-name",
        default="DRUGBANK_RO_PASSWORD",
        help="Name of the secret (default: DRUGBANK_RO_PASSWORD)"
    )
    parser.add_argument(
        "--username",
        help="Database username (will prompt if not provided)"
    )
    parser.add_argument(
        "--password",
        help="Database password (will prompt if not provided)"
    )
    parser.add_argument(
        "--profiles",
        nargs="+",
        default=["cmp-dev", "sci-dev"],
        help="AWS profiles to use (default: cmp-dev sci-dev)"
    )
    parser.add_argument(
        "--region",
        default="us-east-1",
        help="AWS region (default: us-east-1)"
    )

    args = parser.parse_args()

    # Get credentials
    username = args.username or input("Database username: ")
    password = args.password or getpass.getpass("Database password: ")

    if not username or not password:
        print("Error: Username and password are required")
        sys.exit(1)

    print(f"\nStoring credentials in AWS Secrets Manager...")
    print("Secret name: [REDACTED]")
    print(f"Profiles: {', '.join(args.profiles)}")
    print(f"Region: {args.region}\n")

    success_count = 0
    for profile in args.profiles:
        if create_or_update_secret(
            args.secret_name,
            username,
            password,
            profile,
            args.region
        ):
            success_count += 1

    print(f"\n{'='*60}")
    print(f"Successfully configured {success_count}/{len(args.profiles)} profiles")
    print(f"{'='*60}")

    if success_count < len(args.profiles):
        sys.exit(1)


if __name__ == "__main__":
    main()