#!/usr/bin/env python3
"""
Test connection to GOSTAR database and display basic information.
"""

import os
import psycopg2
import polars as pl
from get_gostar_credentials import get_gostar_connection_params, get_gostar_connection_string

def test_psycopg2_connection():
    """Test connection using psycopg2."""
    print("Testing psycopg2 connection...")

    try:
        conn_params = get_gostar_connection_params()
    except Exception as e:
        print(f"✗ Failed to get credentials: {e}")
        return False

    try:
        conn = psycopg2.connect(**conn_params)
        cursor = conn.cursor()

        # Test query
        cursor.execute("SELECT version();")
        version = cursor.fetchone()
        print(f"✓ Connected successfully!")
        print(f"  PostgreSQL version: {version[0]}")

        # Get database size
        cursor.execute("SELECT pg_size_pretty(pg_database_size('gostar'));")
        size = cursor.fetchone()
        print(f"  Database size: {size[0]}")

        cursor.close()
        conn.close()

        return True

    except Exception as e:
        print(f"✗ Connection failed: {e}")
        return False

def test_sqlalchemy_connection():
    """Test connection using Polars."""
    print("\nTesting Polars connection...")

    try:
        conn_string = get_gostar_connection_string()
    except Exception as e:
        print(f"✗ Failed to get credentials: {e}")
        return False

    try:
        # Test query using Polars
        df = pl.read_database_uri("SELECT current_database(), current_user;", uri=conn_string)
        print(f"✓ Connected successfully!")
        print(f"  Database: {df[0, 0]}")
        print(f"  User: {df[0, 1]}")

        return True

    except Exception as e:
        print(f"✗ Connection failed: {e}")
        return False

def list_tables():
    """List all tables in the database."""
    print("\nListing available tables...")

    try:
        conn_string = get_gostar_connection_string()
    except Exception as e:
        print(f"✗ Failed to get credentials: {e}")
        return False

    try:
        query = """
        SELECT
            schemaname,
            tablename,
            pg_size_pretty(pg_total_relation_size(schemaname||'.'||tablename)) as size
        FROM pg_tables
        WHERE schemaname NOT IN ('pg_catalog', 'information_schema')
        ORDER BY schemaname, tablename
        LIMIT 50
        """

        df = pl.read_database_uri(query, uri=conn_string)

        if len(df) > 0:
            print(f"\n✓ Found {len(df)} tables:")
            print(df)
        else:
            print("  No user tables found")

        return True

    except Exception as e:
        print(f"✗ Failed to list tables: {e}")
        return False

def main():
    """Run all tests."""
    print("=" * 60)
    print("GOSTAR Database Connection Test")
    print("=" * 60)

    # Test connections
    psycopg2_ok = test_psycopg2_connection()
    sqlalchemy_ok = test_sqlalchemy_connection()

    # List tables if connection successful
    if psycopg2_ok or sqlalchemy_ok:
        list_tables()

    print("\n" + "=" * 60)
    print("Test completed")
    print("=" * 60)

if __name__ == "__main__":
    main()
