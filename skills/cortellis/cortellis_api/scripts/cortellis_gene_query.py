#!/usr/bin/env python3
"""
Cortellis Gene Query Script
Query Cortellis API for gene annotations and related data.
"""

import os
import sys
import json
import argparse
import requests
import xml.etree.ElementTree as ET
from typing import List, Dict, Any, Optional
from pathlib import Path
from requests.auth import HTTPDigestAuth
from dotenv import load_dotenv
import pandas as pd
from datetime import datetime
import asyncio
import aiohttp

# Load environment variables from .env file
load_dotenv()

# API Credentials from environment
API_KEY = os.getenv("CORTELLIS_API_KEY")
API_SECRET = os.getenv("CORTELLIS_API_SECRET")

if not API_KEY or not API_SECRET:
    print("ERROR: CORTELLIS_API_KEY and CORTELLIS_API_SECRET must be set in .env file")
    sys.exit(1)

# API URLs
CI_AUTH_URL = "https://api.cortellis.com/api-ws/ws/rs/auth-v2/token"
SCI_AUTH_URL = "https://api.cortellis.com/api-ws/ws/rs/auth-v4/token"
DDI_AUTH_URL = "https://api.cortellis.com/api-ws/ws/rs/auth-v1/token"
TARGET_AUTH_URL = "https://api.cortellis.com/api-ws/ws/rs/auth-v3/token"

CORTELLIS_TARGETS_URL = "https://api.cortellis.com/api-ws/ws/rs/targets-v2"
CORTELLIS_BIOMARKER_URL = "https://api.cortellis.com/api-ws/ws/rs/biomarkers-v3"
CORTELLIS_DRUGS_URL = "https://api.cortellis.com/api-ws/ws/rs/drugs-v2"


def fetch_token(url: str, user: str = API_KEY, pw: str = API_SECRET) -> str:
    """Return token handling both XML and plain-text replies."""
    r = requests.get(url, auth=HTTPDigestAuth(user, pw), timeout=120)
    r.raise_for_status()

    ctype = r.headers.get("Content-Type", "")
    body = r.text.strip()

    if ctype.startswith("application/xml") or body.startswith("<"):
        root = ET.fromstring(body)
        token = root.findtext("token")
    else:
        token = body

    if not token:
        raise RuntimeError("Could not extract token - check credentials/endpoint")
    return token


# Fetch authentication tokens
print("Authenticating with Cortellis API...")
TARGET_HEADER = {"API-Token": fetch_token(TARGET_AUTH_URL)}
BIOMARKER_HEADER = {"API-Token": fetch_token(SCI_AUTH_URL)}
DRUGS_HEADER = {"API-Token": fetch_token(CI_AUTH_URL)}  # Investigational Drugs uses CI auth
print("Authentication successful.\n")


def as_list(x):
    """Convert to list if not already."""
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def join_ids(ids_iterable):
    """Join IDs into comma-separated string."""
    ids = [s for s in {i for i in ids_iterable if i} if isinstance(s, str)]
    return ",".join(ids)


async def fetch_url_async(session: aiohttp.ClientSession, url: str, headers: Dict, timeout_seconds: int = 120) -> Optional[Dict]:
    """Async fetch a single URL and return JSON response."""
    try:
        timeout = aiohttp.ClientTimeout(total=timeout_seconds)
        async with session.get(url, headers=headers, timeout=timeout) as response:
            if response.status == 200:
                return await response.json()
            else:
                return None
    except Exception as e:
        return None


async def fetch_chunks_async(urls: List[str], headers: Dict, max_concurrent: int = 10, timeout_seconds: int = 120) -> List[Optional[Dict]]:
    """Fetch multiple URLs concurrently with concurrency limit."""
    semaphore = asyncio.Semaphore(max_concurrent)

    async def fetch_with_semaphore(url):
        async with semaphore:
            return await fetch_url_async(session, url, headers, timeout_seconds)

    async with aiohttp.ClientSession() as session:
        tasks = [fetch_with_semaphore(url) for url in urls]
        results = await asyncio.gather(*tasks, return_exceptions=True)
        # Convert exceptions to None
        return [r if not isinstance(r, Exception) else None for r in results]


def fetch_with_adaptive_concurrency(urls: List[str], headers: Dict, initial_concurrent: int = 25, fallback_concurrent: int = 20, timeout_seconds: int = 120) -> List[Optional[Dict]]:
    """
    Fetch URLs with adaptive concurrency and exponential backoff.

    Retry strategy:
    - Try initial_concurrent (e.g., 25)
    - If >1% fail: Wait 2s, retry at fallback_concurrent (e.g., 20)
    - If still failures: Wait 5s, retry at 15
    - If still failures: Wait 10s, retry at 10
    - If still failures: Wait 15s, retry at 5

    Args:
        urls: List of URLs to fetch
        headers: HTTP headers
        initial_concurrent: Initial concurrency level (default: 25)
        fallback_concurrent: Fallback concurrency if errors detected (default: 20)
        timeout_seconds: Request timeout

    Returns:
        List of response dictionaries or None for failures
    """
    import time

    # Try with initial concurrency
    results = asyncio.run(fetch_chunks_async(urls, headers, initial_concurrent, timeout_seconds))

    # Define retry cascade: (wait_seconds, concurrency_level)
    retry_cascade = [
        (2, fallback_concurrent),  # Stage 1: Wait 2s, try at 20
        (5, 15),                   # Stage 2: Wait 5s, try at 15
        (10, 10),                  # Stage 3: Wait 10s, try at 10
        (15, 5),                   # Stage 4: Wait 15s, try at 5
    ]

    for stage_num, (wait_seconds, concurrency) in enumerate(retry_cascade, 1):
        # Count current failures
        failure_count = sum(1 for r in results if r is None)
        failure_rate = failure_count / len(results) if results else 0

        # If failure rate > 1%, retry failed requests (very aggressive)
        if failure_rate > 0.01 and failure_count > 0:
            print(f"    Stage {stage_num}: Detected {failure_count} failures ({failure_rate:.1%}), waiting {wait_seconds}s then retrying at {concurrency} concurrent...")

            # Wait before retry (exponential backoff)
            time.sleep(wait_seconds)

            # Find indices of failed requests
            failed_indices = [i for i, r in enumerate(results) if r is None]
            failed_urls = [urls[i] for i in failed_indices]

            # Retry failed requests with lower concurrency
            retry_results = asyncio.run(fetch_chunks_async(failed_urls, headers, concurrency, timeout_seconds))

            # Merge results back
            for idx, retry_result in zip(failed_indices, retry_results):
                results[idx] = retry_result
        else:
            # Success rate acceptable, no need for further retries
            break

    # Final failure check
    final_failures = sum(1 for r in results if r is None)
    if final_failures > 0:
        print(f"    Warning: {final_failures}/{len(results)} requests still failed after all retries")

    return results


def query_cortellis_target_info(
    query_param: str,
    hits: int = 30,
    offset: int = 0,
    sort_by: str = "targetNameMain",
    sort_direction: str = "ascending",
    output_format: str = "json",
    headers: Dict = TARGET_HEADER
) -> Dict:
    """Query Cortellis Targets API."""
    params = {
        "query": query_param,
        "hits": hits,
        "offset": offset,
        "sortBy": sort_by,
        "sortDirection": sort_direction,
        "fmt": output_format
    }

    req = requests.Request('GET', f"{CORTELLIS_TARGETS_URL}/target/search", params=params)
    query_url = req.prepare().url

    print(f"  Querying: {query_param}")

    response = requests.get(query_url, headers=headers, timeout=120)
    response.raise_for_status()

    if output_format == "json":
        return response.json()
    else:
        return response.text


def target_summaryrecord(target: str) -> tuple[Dict, str]:
    """Get target summary record and return (summary_dict, target_id)."""
    target_summary_json = query_cortellis_target_info(
        query_param=f"targetSynonyms:{target}",
        hits=50,
        sort_by="targetNameMain"
    )

    sr = (target_summary_json or {}).get("TargetResultsOutput", {}).get("SearchResults", {})
    results = as_list(sr.get("TargetResult"))

    if not results:
        raise ValueError(f"No TargetResult for query targetSynonyms:{target}")

    lower_t = str(target).lower()

    def score(r):
        # Check main name fields
        name_fields = [
            str(r.get("targetNameMain", "")).lower(),
            str(r.get("@DisplayName", "")).lower(),
            str(r.get("Symbol", "")).lower(),
            str(r.get("NameMain", "")).lower(),
            str(r.get("@namemain", "")).lower(),
        ]

        # Check if target matches any main field exactly
        if any(lower_t == f for f in name_fields if f):
            return True

        # Check synonyms for exact match
        synonyms_data = r.get('Synonyms', {})
        if isinstance(synonyms_data, dict):
            synonyms = as_list(synonyms_data.get('Synonym', []))
        else:
            synonyms = as_list(synonyms_data)

        for syn in synonyms:
            syn_str = str(syn).lower().strip()
            # Exact match
            if syn_str == lower_t:
                return True

        return False

    picked = next((r for r in results if score(r)), results[0])

    # Get target ID
    tid = picked.get("@Id") or picked.get("@id")
    if not tid:
        raise KeyError("No @Id/@id in TargetResult")

    # Get full target record
    target_record_url = requests.Request(
        'GET',
        f"{CORTELLIS_TARGETS_URL}/targets",
        params={'idList': tid, 'fmt': 'json'}
    ).prepare().url

    target_record = requests.get(target_record_url, headers=TARGET_HEADER, timeout=120).json()
    record_target = (target_record or {}).get('TargetRecordsOutput', {}).get('Targets', {}).get('Target', {})

    # Merge selected summary & full record
    target_summary = dict(picked)
    if isinstance(record_target, dict):
        target_summary.update(record_target)

    return target_summary, tid


def target_interaction(target_id: str) -> Dict:
    """Get target interactions."""
    target_interaction_url = requests.Request(
        'GET',
        f"{CORTELLIS_TARGETS_URL}/target/interactions",
        params={'idList': target_id, 'fmt': 'json'}
    ).prepare().url

    target_inter = requests.get(target_interaction_url, headers=TARGET_HEADER, timeout=120).json()
    return target_inter


def target_drug(target_id: str, annotation: Optional[Dict] = None) -> Dict:
    """Get condition drug associations, drug records, and trial records.

    Args:
        target_id: The main target ID
        annotation: Optional annotation data containing RelatedTargets
    """
    # Collect target IDs to query (main target + related targets)
    target_ids_to_query = [target_id]

    # Extract related target IDs if annotation is provided
    if annotation:
        related_targets = annotation.get('RelatedTargets', {})
        if related_targets:
            related_ids = as_list(related_targets.get('Id', []))
            for rel_id in related_ids:
                if isinstance(rel_id, dict):
                    rel_target_id = rel_id.get('$')
                    if rel_target_id:
                        target_ids_to_query.append(rel_target_id)
                elif isinstance(rel_id, str):
                    target_ids_to_query.append(rel_id)

    if len(target_ids_to_query) > 1:
        print(f"  Querying {len(target_ids_to_query)} targets (main + {len(target_ids_to_query)-1} related)")

    # Get ConditionDrugAssociations for all targets
    all_drug_ids = set()
    target_disdrug = None

    for tid in target_ids_to_query:
        target_disdrug_url = requests.Request(
            'GET',
            f"{CORTELLIS_TARGETS_URL}/target/conditionDrugAssociations",
            params={'idList': tid, 'fmt': 'json'}
        ).prepare().url

        try:
            response = requests.get(target_disdrug_url, headers=TARGET_HEADER, timeout=120).json()
            tid_data = response['TargetRecordsOutput']['Targets']

            # Keep the main target's data for return
            if tid == target_id:
                target_disdrug = tid_data

            # Extract drug IDs from this target
            conditions = (((tid_data.get('Target') or {})
                           .get('ConditionDrugAssociations') or {})
                           .get('Condition') or [])

            for condition in as_list(conditions):
                drug_ids = condition.get('DrugId')
                if isinstance(drug_ids, dict):
                    val = drug_ids.get('$')
                    if val:
                        all_drug_ids.add(val)
                else:
                    for d in as_list(drug_ids):
                        if isinstance(d, dict):
                            val = d.get('$')
                            if val:
                                all_drug_ids.add(val)
                        elif isinstance(d, str):
                            all_drug_ids.add(d)
        except Exception as e:
            print(f"    Warning: Error querying related target {tid}: {str(e)[:100]}")
            continue

    # Fallback if main target data wasn't set
    if target_disdrug is None:
        target_disdrug = {'Target': {'ConditionDrugAssociations': {}}}

    # Get Drug Records using async bulk fetch
    target_drug = {'Drug': []}

    if all_drug_ids:
        ids_list = list(all_drug_ids)
        chunk_size = 30
        total_chunks = (len(ids_list) + chunk_size - 1) // chunk_size
        print(f"  Fetching {len(ids_list)} drug records in {total_chunks} chunks (async)...")

        # Prepare all chunk URLs
        chunk_urls = []
        for i in range(0, len(ids_list), chunk_size):
            chunk = ids_list[i:i + chunk_size]
            drug_id_str = ",".join(chunk)
            target_drug_url = requests.Request(
                'GET',
                f"{CORTELLIS_TARGETS_URL}/drugs",
                params={'idList': drug_id_str, 'fmt': 'json'}
            ).prepare().url
            chunk_urls.append(target_drug_url)

        # Fetch all chunks concurrently with adaptive concurrency (25→20)
        results = fetch_with_adaptive_concurrency(chunk_urls, TARGET_HEADER, initial_concurrent=25, fallback_concurrent=20)

        # Process results
        failed_chunks = 0
        for idx, resp in enumerate(results):
            if resp:
                payload = (resp or {}).get('drugRecordsOutput', {}).get('Drug')
                if payload:
                    target_drug['Drug'].extend(payload if isinstance(payload, list) else [payload])
            else:
                failed_chunks += 1

        # Track data loss
        data_loss_tracker['failed_drug_chunks'] += failed_chunks

        if failed_chunks > 0:
            print(f"    Warning: {failed_chunks}/{total_chunks} chunks failed (after retry)")

        print(f"  Retrieved {len(target_drug['Drug'])} drug records")

    # Get Trial Records using async bulk fetch
    target_trail = {'Trial': []}
    all_trail_ids = set()
    for drug in target_drug.get('Drug', []):
        for trial in as_list((drug.get('RelatedTrials') or {}).get('Trial')):
            tid = trial.get('@id')
            if tid:
                all_trail_ids.add(tid)

    if all_trail_ids:
        ids_list = list(all_trail_ids)
        chunk_size = 30
        total_chunks = (len(ids_list) + chunk_size - 1) // chunk_size
        print(f"  Fetching {len(ids_list)} trial records in {total_chunks} chunks (async)...")

        # Prepare all chunk URLs
        chunk_urls = []
        for i in range(0, len(ids_list), chunk_size):
            chunk = ids_list[i:i + chunk_size]
            trail_id_str = ",".join(chunk)
            target_trail_url = requests.Request(
                'GET',
                f"{CORTELLIS_TARGETS_URL}/target/trials",
                params={'idList': trail_id_str, 'fmt': 'json'}
            ).prepare().url
            chunk_urls.append(target_trail_url)

        # Fetch all chunks concurrently with adaptive concurrency (25→20)
        results = fetch_with_adaptive_concurrency(chunk_urls, TARGET_HEADER, initial_concurrent=25, fallback_concurrent=20)

        # Process results
        failed_chunks = 0
        for resp in results:
            if resp:
                payload = (resp or {}).get('TrialRecordsOutput', {}).get('Trial')
                if payload:
                    target_trail['Trial'].extend(payload if isinstance(payload, list) else [payload])
            else:
                failed_chunks += 1

        # Track data loss
        data_loss_tracker['failed_trial_chunks'] += failed_chunks

        if failed_chunks > 0:
            print(f"    Warning: {failed_chunks}/{total_chunks} trial chunks failed (after retry)")

        print(f"  Retrieved {len(target_trail['Trial'])} trial records")

    # Filter out drugs whose @namemain is all numbers (these are typically internal IDs)
    filtered_drugs = []
    for drug in target_drug.get('Drug', []):
        drug_name = drug.get('@namemain', '')
        # Keep drug if name is not all numbers
        if drug_name and not drug_name.isdigit():
            filtered_drugs.append(drug)

    print(f"  Filtered out {len(target_drug.get('Drug', [])) - len(filtered_drugs)} drugs with numeric-only names")

    # Update target_drug to only contain filtered drugs
    target_drug['Drug'] = filtered_drugs

    # Get comprehensive drug records from Investigational Drugs API
    comprehensive_drug_records = get_comprehensive_drug_records(filtered_drugs, max_records=10000)

    return target_disdrug | target_drug | target_trail | {'DrugRecord': comprehensive_drug_records}


def get_comprehensive_drug_records(basic_drug_records: List[Dict], max_records: int = 10000) -> Dict:
    """
    Query Investigational Drugs API v2.0 for comprehensive drug records.

    Uses drug search with @namemain from Targets API to get ALL associated drug IDs,
    then bulk fetches comprehensive records via /drugs?idList= endpoint using async.

    Args:
        basic_drug_records: List of basic drug records from Targets API
        max_records: Maximum number of drug names to search (default: 10000)

    Returns dict with drug ID as key and full drug record as value.
    """
    comprehensive_records = {}

    if not basic_drug_records:
        return comprehensive_records

    # Step 1: Search for all drug IDs associated with each drug name (async)
    records_to_search = basic_drug_records[:max_records]
    print(f"  Searching for drug IDs from {len(records_to_search)} drug names (async)...")

    # Prepare search URLs and metadata
    search_urls = []
    drug_metadata = []  # List of (drug_namemain, drug_name) tuples

    for drug in records_to_search:
        drug_namemain = drug.get('@namemain', '')
        if not drug_namemain:
            continue

        # Extract drug name for reference
        drug_name = ''
        names_data = drug.get('NamesChemicalAndDescriptions', {})
        if isinstance(names_data, dict):
            name_entry = names_data.get('Name', '')
            if isinstance(name_entry, dict):
                drug_name = extract_value(name_entry)
            elif isinstance(name_entry, str):
                drug_name = name_entry
        if not drug_name:
            drug_name = drug.get('displayName', drug_namemain)

        search_url = requests.Request(
            'GET',
            f"{CORTELLIS_DRUGS_URL}/drug/search",
            params={
                'query': f'drugNamesAll:{drug_namemain}',
                'hits': 2000,
                'returnFilterCount': 5000,
                'fmt': 'json'
            }
        ).prepare().url

        search_urls.append(search_url)
        drug_metadata.append((drug_namemain, drug_name))

    # Fetch all drug searches concurrently with adaptive concurrency (25→20)
    search_results = fetch_with_adaptive_concurrency(search_urls, DRUGS_HEADER, initial_concurrent=25, fallback_concurrent=20)

    # Process search results to collect drug IDs
    all_drug_ids = {}  # Maps drug_id -> (source_drug_name, source_display_name)
    failed_drug_names = []  # Track which specific drugs failed to search
    no_results_drug_names = []  # Track which drugs returned no results

    for idx, (resp, (drug_namemain, drug_name)) in enumerate(zip(search_results, drug_metadata)):
        if resp and isinstance(resp, dict):
            drug_results_output = resp.get('drugResultsOutput', {})
            if isinstance(drug_results_output, dict):
                search_results_data = drug_results_output.get('SearchResults', {})

                if isinstance(search_results_data, dict):
                    # Try both 'DrugResult' and 'Drug' keys
                    drug_result_data = search_results_data.get('DrugResult') or search_results_data.get('Drug', [])
                    results = as_list(drug_result_data)

                    if results:
                        # Collect ALL drug IDs from search results
                        for result in results:
                            if isinstance(result, dict):
                                drug_id = result.get('@id')
                                if drug_id and drug_id not in all_drug_ids:
                                    all_drug_ids[drug_id] = (drug_namemain, drug_name)
                    else:
                        # Search succeeded but returned 0 results
                        no_results_drug_names.append(drug_namemain)
        else:
            # Search request failed (timeout/error)
            failed_drug_names.append(drug_namemain)

    # Track data loss
    data_loss_tracker['failed_drug_searches'].extend(failed_drug_names)
    data_loss_tracker['no_results_drug_names'].extend(no_results_drug_names)

    # Report on failures and no-results
    if failed_drug_names:
        print(f"  ⚠️  Failed drug searches (timeout/error): {len(failed_drug_names)}")
        for name in failed_drug_names[:10]:  # Show first 10
            print(f"      - {name}")
    if no_results_drug_names:
        print(f"  ℹ️  Drug searches with no results: {len(no_results_drug_names)}")
        if len(no_results_drug_names) <= 5:
            for name in no_results_drug_names:
                print(f"      - {name}")

    print(f"  Found {len(all_drug_ids)} unique drug IDs from search")

    # Step 2: Bulk query comprehensive drug records using /drugs?idList= endpoint (async)
    if all_drug_ids:
        drug_ids_list = list(all_drug_ids.keys())
        chunk_size = 30  # Bulk query in chunks of 30
        total_chunks = (len(drug_ids_list) + chunk_size - 1) // chunk_size
        print(f"  Fetching {len(drug_ids_list)} comprehensive drug records in {total_chunks} chunks (async)...")

        # Prepare all chunk URLs
        chunk_urls = []
        for i in range(0, len(drug_ids_list), chunk_size):
            chunk = drug_ids_list[i:i + chunk_size]
            drug_id_str = ",".join(chunk)

            drugs_url = requests.Request(
                'GET',
                f"{CORTELLIS_DRUGS_URL}/drugs",
                params={'idList': drug_id_str, 'fmt': 'json'}
            ).prepare().url
            chunk_urls.append(drugs_url)

        # Fetch all chunks concurrently with adaptive concurrency (25→20)
        results = fetch_with_adaptive_concurrency(chunk_urls, DRUGS_HEADER, initial_concurrent=25, fallback_concurrent=20)

        # Process results
        failed_chunks = 0
        for resp in results:
            if resp:
                drug_records_output = resp.get('drugRecordsOutput', {})
                drugs = as_list(drug_records_output.get('Drug', []))

                for drug_record in drugs:
                    if isinstance(drug_record, dict):
                        record_id = drug_record.get('@id')
                        if record_id:
                            # Add source drug name from Targets API
                            source_names = all_drug_ids.get(record_id, ('', ''))
                            drug_record['_source_drug_name'] = source_names[0]
                            drug_record['_source_display_name'] = source_names[1]
                            comprehensive_records[record_id] = drug_record
            else:
                failed_chunks += 1

        # Track data loss
        data_loss_tracker['failed_comprehensive_drug_chunks'] += failed_chunks

        if failed_chunks > 0:
            print(f"    Warning: {failed_chunks}/{total_chunks} comprehensive drug chunks failed (after retry)")

    print(f"  Retrieved {len(comprehensive_records)} comprehensive drug records")
    return comprehensive_records


def target_association(target_id: str) -> Dict:
    """Get condition gene associations (Target field renamed to association)."""
    # Get ConditionGeneAssociations
    target_disbiomarker_url = requests.Request(
        'GET',
        f"{CORTELLIS_TARGETS_URL}/target/conditionGeneAssociations",
        params={'idList': target_id, 'fmt': 'json'}
    ).prepare().url

    target_disbiomarker = requests.get(target_disbiomarker_url, headers=TARGET_HEADER, timeout=120).json()['TargetRecordsOutput']['Targets']

    # Get ConditionGeneVariantAssociations
    target_disgenevar_url = requests.Request(
        'GET',
        f"{CORTELLIS_TARGETS_URL}/target/conditionGeneVariantAssociations",
        params={'idList': target_id, 'fmt': 'json'}
    ).prepare().url

    target_disbiomarker['Target']['ConditionGeneVariantAssociations'] = requests.get(
        target_disgenevar_url, headers=TARGET_HEADER, timeout=120
    ).json()['TargetRecordsOutput']['Targets']['Target']['ConditionGeneVariantAssociations']

    return target_disbiomarker


def query_cortellis_biomarker_info(
    query_param: str,
    hits: int = 30,
    offset: int = 0,
    sort_direction: str = "ascending",
    output_format: str = "json",
    headers: Dict = BIOMARKER_HEADER
) -> Dict:
    """Query Cortellis Biomarker API."""
    params = {
        "query": query_param,
        "hits": hits,
        "offset": offset,
        "sortDirection": sort_direction,
        "fmt": output_format
    }

    req = requests.Request('GET', f"{CORTELLIS_BIOMARKER_URL}/biomarkerUse/search", params=params)
    query_url = req.prepare().url

    print(f"  Querying biomarker: {query_param}")

    response = requests.get(query_url, headers=headers, timeout=120)
    response.raise_for_status()

    if output_format == "json":
        return response.json()
    else:
        return response.text


def target_biomarker(target: str) -> Dict:
    """Get biomarker information."""
    bm_json = query_cortellis_biomarker_info(
        query_param=f"biomarkerSynonyms:{target}",
        hits=1000
    )

    sr = (bm_json or {}).get('biomarkerUseResultsOutput', {}).get('SearchResults', {})
    results = as_list(sr.get('BiomarkerUseResult'))

    if not results:
        return sr

    # Get BiomarkerUseRecords
    biomarkeruse_ids = [r.get('@id') for r in results if r.get('@id')]
    if biomarkeruse_ids:
        biomarkeruse_url = requests.Request(
            'GET',
            f"{CORTELLIS_BIOMARKER_URL}/biomarkerUses",
            params={'idList': join_ids(biomarkeruse_ids), 'fmt': 'json'}
        ).prepare().url

        use_payload = requests.get(biomarkeruse_url, headers=BIOMARKER_HEADER, timeout=120).json()
        sr['BiomarkerUse'] = (use_payload or {}).get('biomarkerUseRecordsOutput', {}).get('BiomarkerUse')

    # Get BiomarkerRecords
    biomarker_ids = []
    for r in results:
        b = r.get('Biomarker') or {}
        bid = b.get('@id')
        if bid:
            biomarker_ids.append(bid)

    if biomarker_ids:
        biomarkerecord_url = requests.Request(
            'GET',
            f"{CORTELLIS_BIOMARKER_URL}/biomarkers",
            params={'idList': join_ids(biomarker_ids), 'fmt': 'json'}
        ).prepare().url

        rec_payload = requests.get(biomarkerecord_url, headers=BIOMARKER_HEADER, timeout=120).json()
        sr['Biomarker'] = (rec_payload or {}).get('BiomarkerRecordsOutput', {}).get('Biomarker')

    return sr


def summarize_annotation(data: Dict) -> str:
    """Summarize annotation field in markdown format."""
    summary = []
    summary.append("## Annotation\n")

    # Handle both old and new API response formats
    target_name = (data.get('targetNameMain') or data.get('NameMain') or
                   data.get('@DisplayName') or data.get('@namemain', 'N/A'))

    symbol = data.get('Symbol', 'N/A')

    # Extract Gene ID
    gene_id = 'N/A'
    if 'EntrezgeneIdentifiers' in data:
        gene_id_data = data['EntrezgeneIdentifiers'].get('Identifier', '')
        gene_id = str(gene_id_data) if gene_id_data else 'N/A'
    elif 'GeneId' in data:
        gene_id = extract_value(data.get('GeneId'), 'N/A')

    # Extract UniProt ID
    uniprot_id = 'N/A'
    if 'ExternalIdentifiers' in data:
        identifiers = as_list(data['ExternalIdentifiers'].get('Identifier', []))
        for ident in identifiers:
            if isinstance(ident, dict) and ident.get('@type') == 'SwissProt':
                uniprot_id = extract_value(ident)
                break
    elif 'UniprotId' in data:
        uniprot_id = extract_value(data.get('UniprotId'), 'N/A')

    summary.append(f"- **Target Name:** {target_name}")
    summary.append(f"- **Symbol:** {symbol}")
    summary.append(f"- **Gene ID:** {gene_id}")
    summary.append(f"- **UniProt ID:** {uniprot_id}")

    # Target type
    target_types = as_list(data.get('TargetType') or data.get('Type', []))
    if target_types:
        types = [extract_value(t) for t in target_types]
        summary.append(f"- **Target Type:** {', '.join(types)}")

    # Description
    description = ''
    if 'Descriptions' in data:
        descriptions = as_list(data['Descriptions'].get('Description', []))
        if descriptions:
            description = extract_value(descriptions[0])
    elif 'Description' in data:
        description = extract_value(data.get('Description'))

    if description:
        summary.append(f"- **Description:** {description[:300]}...")

    # Organism
    organism = data.get('Organism', {})
    organism_name = extract_value(organism)
    if organism_name:
        summary.append(f"- **Organism:** {organism_name}")

    return "\n".join(summary)


def summarize_interaction(data: Dict) -> str:
    """Summarize interaction field in markdown format."""
    summary = []
    summary.append("## Interactions\n")

    target = data.get('TargetRecordsOutput', {}).get('Targets', {}).get('Target', {})
    interactions = target.get('Interactions', {}).get('Interaction', [])
    interactions = as_list(interactions)

    if not interactions:
        summary.append("*No interactions found*")
        return "\n".join(summary)

    summary.append(f"**Total Interactions:** {len(interactions)}\n")

    # Group by type
    interaction_types = {}
    for interaction in interactions:
        itype = interaction.get('Type', {}).get('$', 'Unknown')
        interaction_types[itype] = interaction_types.get(itype, 0) + 1

    summary.append("### By Type:\n")
    for itype, count in sorted(interaction_types.items(), key=lambda x: x[1], reverse=True):
        summary.append(f"- **{itype}:** {count}")

    return "\n".join(summary)


def summarize_drug(data: Dict) -> str:
    """Summarize drug field in markdown format showing both basic and comprehensive records."""
    summary = []
    summary.append("## Drugs\n")

    basic_drugs = data.get('Drug', [])
    comprehensive_records = data.get('DrugRecord', {})

    # Show basic drug count from Targets API
    if basic_drugs:
        summary.append(f"**Total Drugs (Targets API):** {len(basic_drugs)}")

    # Show comprehensive drug count from Investigational Drugs API
    if comprehensive_records:
        summary.append(f"**Comprehensive Drug Records (Investigational Drugs API):** {len(comprehensive_records)}\n")
        summary.append("### Top Comprehensive Drug Records:\n")

        for i, (drug_name_key, drug_record) in enumerate(list(comprehensive_records.items())[:5], 1):
            drug_name = drug_record.get('DrugName', drug_name_key)

            # Phase
            phase_highest = drug_record.get('PhaseHighest', {})
            phase = extract_value(phase_highest) if isinstance(phase_highest, dict) else phase_highest

            # Company
            originator = drug_record.get('CompanyOriginator', {})
            company = extract_value(originator) if originator else 'N/A'

            # Primary indication
            indications_primary = drug_record.get('IndicationsPrimary', {})
            indication = 'N/A'
            if isinstance(indications_primary, dict):
                inds = as_list(indications_primary.get('Indication', []))
                if inds:
                    ind_name = inds[0].get('IndicationName', {}) if isinstance(inds[0], dict) else inds[0]
                    indication = extract_value(ind_name) if isinstance(ind_name, dict) else ind_name

            summary.append(f"{i}. **{drug_name}**")
            summary.append(f"   - Phase: {phase}")
            summary.append(f"   - Originator: {company}")
            summary.append(f"   - Primary Indication: {indication}")

        if len(comprehensive_records) > 5:
            summary.append(f"\n*... and {len(comprehensive_records) - 5} more comprehensive records*")
    else:
        summary.append(f"**Comprehensive Drug Records (Investigational Drugs API):** 0")

    # Trials
    trials = data.get('Trial', [])
    if trials:
        summary.append(f"\n**Total Clinical Trials:** {len(trials)}")

    return "\n".join(summary)


def summarize_biomarker(data: Dict) -> str:
    """Summarize biomarker field in markdown format."""
    summary = []
    summary.append("## Biomarkers\n")

    biomarker_uses = as_list(data.get('BiomarkerUseResult', []))
    if not biomarker_uses:
        summary.append("*No biomarker uses found*")
        return "\n".join(summary)

    summary.append(f"**Total Biomarker Uses:** {len(biomarker_uses)}\n")

    # Group by biomarker application
    applications = {}
    for use in biomarker_uses:
        app = use.get('BiomarkerApplication', {}).get('$', 'Unknown')
        applications[app] = applications.get(app, 0) + 1

    summary.append("### By Application:\n")
    for app, count in sorted(applications.items(), key=lambda x: x[1], reverse=True):
        summary.append(f"- **{app}:** {count}")

    # Sample biomarkers
    biomarkers = as_list(data.get('Biomarker', []))
    if biomarkers:
        summary.append(f"\n### Sample Biomarkers:\n")
        for i, bm in enumerate(biomarkers[:3]):
            bm_name = bm.get('biomarkerName', {}).get('$', 'N/A')
            bm_type = bm.get('BiomarkerType', {}).get('$', 'N/A')
            summary.append(f"{i+1}. **{bm_name}** (Type: {bm_type})")

    return "\n".join(summary)


def summarize_association(data: Dict) -> str:
    """Summarize association field (formerly Target) in markdown format."""
    summary = []
    summary.append("## Associations\n")

    target = data.get('Target', {})

    # Condition Gene Associations
    cond_gene_assoc = target.get('ConditionGeneAssociations', {})
    conditions = as_list(cond_gene_assoc.get('Condition', []))

    if not conditions:
        summary.append("*No condition-gene associations found*")
    else:
        summary.append(f"**Total Condition-Gene Associations:** {len(conditions)}\n")

        # Sample conditions
        summary.append("### Top Associations:\n")
        for i, cond in enumerate(conditions[:5]):
            cond_name = cond.get('ConditionName', {}).get('$', 'N/A')
            association_type = cond.get('AssociationType', {}).get('$', 'N/A')
            summary.append(f"{i+1}. **{cond_name}** ({association_type})")

    # Condition Gene Variant Associations
    cond_var_assoc = target.get('ConditionGeneVariantAssociations', {})
    variants = as_list(cond_var_assoc.get('Condition', []))

    if variants:
        summary.append(f"\n**Total Condition-Gene-Variant Associations:** {len(variants)}")

    return "\n".join(summary)


def extract_value(obj, default=''):
    """Extract value from dict with '$' key or return the object itself."""
    if isinstance(obj, dict):
        return obj.get('$', default)
    return obj if obj is not None else default


def flatten_annotation_to_df(data: Dict) -> pd.DataFrame:
    """Flatten annotation data to DataFrame."""
    rows = []

    # Handle both old and new API response formats
    target_name = (data.get('targetNameMain') or data.get('NameMain') or
                   data.get('@DisplayName') or data.get('@namemain', ''))

    # Extract Gene ID from EntrezgeneIdentifiers
    gene_id = ''
    if 'EntrezgeneIdentifiers' in data:
        gene_id_data = data['EntrezgeneIdentifiers'].get('Identifier', '')
        gene_id = str(gene_id_data) if gene_id_data else ''
    elif 'GeneId' in data:
        gene_id = extract_value(data.get('GeneId'))

    # Extract UniProt ID from ExternalIdentifiers
    uniprot_id = ''
    if 'ExternalIdentifiers' in data:
        identifiers = as_list(data['ExternalIdentifiers'].get('Identifier', []))
        for ident in identifiers:
            if isinstance(ident, dict) and ident.get('@type') == 'SwissProt':
                uniprot_id = extract_value(ident)
                break
    elif 'UniprotId' in data:
        uniprot_id = extract_value(data.get('UniprotId'))

    # Get descriptions
    description = ''
    if 'Descriptions' in data:
        descriptions = as_list(data['Descriptions'].get('Description', []))
        if descriptions:
            description = extract_value(descriptions[0])
    elif 'Description' in data:
        description = extract_value(data.get('Description'))

    row = {
        'Target ID': data.get('@Id') or data.get('@id', ''),
        'Target Name': target_name,
        'Symbol': data.get('Symbol', ''),
        'Gene ID': gene_id,
        'UniProt ID': uniprot_id,
        'Description': description[:500] if description else '',
    }

    # Target types
    target_types = as_list(data.get('TargetType') or data.get('Type', []))
    types = [extract_value(t) for t in target_types]
    row['Target Types'] = ', '.join(types) if types else ''

    # Synonyms
    synonyms_data = data.get('targetSynonyms') or data.get('Synonyms', {})
    if isinstance(synonyms_data, dict):
        synonyms = as_list(synonyms_data.get('Synonym', []))
    else:
        synonyms = as_list(synonyms_data)
    syn_list = [extract_value(s) for s in synonyms]
    row['Synonyms'] = ', '.join(syn_list[:10]) if syn_list else ''

    # Organism
    organism = data.get('Organism', {})
    row['Organism'] = extract_value(organism)

    rows.append(row)
    return pd.DataFrame(rows)


def flatten_interaction_to_df(data: Dict) -> pd.DataFrame:
    """Flatten interaction data to DataFrame."""
    rows = []
    target = data.get('TargetRecordsOutput', {}).get('Targets', {}).get('Target', {})
    interactions = as_list(target.get('Interactions', {}).get('Interaction', []))

    for interaction in interactions:
        row = {
            'Interaction Type': extract_value(interaction.get('Type')),
            'Partner Name': extract_value(interaction.get('PartnerName')),
            'Partner Gene Symbol': extract_value(interaction.get('PartnerGeneSymbol')),
            'Partner UniProt ID': extract_value(interaction.get('PartnerUniprotId')),
            'Evidence': extract_value(interaction.get('Evidence')),
            'PubMed ID': extract_value(interaction.get('PubmedId')),
        }
        rows.append(row)

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=['Interaction Type', 'Partner Name', 'Partner Gene Symbol', 'Partner UniProt ID', 'Evidence', 'PubMed ID'])


def flatten_basic_drugs_to_df(data: Dict) -> pd.DataFrame:
    """Flatten basic drug records from Targets API to DataFrame."""
    rows = []

    # Use basic drug records from Targets API
    basic_drugs = data.get('Drug', [])

    if basic_drugs:
        for drug in basic_drugs:
            row = {
                'Drug ID': drug.get('@id', ''),
                'Drug Name': drug.get('@namemain', ''),
            }

            # Extract chemical/description name
            names_data = drug.get('NamesChemicalAndDescriptions', {})
            if isinstance(names_data, dict):
                name_entry = names_data.get('Name', '')
                if isinstance(name_entry, dict):
                    row['Chemical Name'] = extract_value(name_entry)
                elif isinstance(name_entry, str):
                    row['Chemical Name'] = name_entry
                else:
                    row['Chemical Name'] = ''
            else:
                row['Chemical Name'] = ''

            # Extract code name
            names_code = drug.get('NamesCode', {})
            if isinstance(names_code, dict):
                code_entry = names_code.get('Name', '')
                row['Code Name'] = extract_value(code_entry) if isinstance(code_entry, dict) else str(code_entry) if code_entry else ''
            else:
                row['Code Name'] = ''

            # Extract molecular mechanism
            mechanisms = drug.get('MechanismsMolecular', {})
            if isinstance(mechanisms, dict):
                mech_entry = mechanisms.get('Mechanism', {})
                row['Molecular Mechanism'] = extract_value(mech_entry) if mech_entry else ''
            else:
                row['Molecular Mechanism'] = ''

            rows.append(row)

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=['Drug ID', 'Drug Name', 'Chemical Name', 'Code Name', 'Molecular Mechanism'])


def flatten_drug_to_df(data: Dict) -> pd.DataFrame:
    """Flatten comprehensive drug records to DataFrame."""
    rows = []

    # Use comprehensive drug records from Investigational Drugs API
    comprehensive_records = data.get('DrugRecord', {})

    if comprehensive_records:
        for drug_id_key, drug_record in comprehensive_records.items():
            row = {
                'Drug ID': drug_id_key,
                'Drug Name': drug_record.get('DrugName', ''),
                'Source Drug Name (Targets API)': drug_record.get('_source_drug_name', ''),
                'Source Display Name': drug_record.get('_source_display_name', ''),
            }

            # Phase information
            phase_highest = drug_record.get('PhaseHighest', {})
            if isinstance(phase_highest, dict):
                row['Highest Phase'] = extract_value(phase_highest)
            else:
                row['Highest Phase'] = phase_highest if phase_highest else ''

            # Company Originator
            originator = drug_record.get('CompanyOriginator', {})
            row['Originator'] = extract_value(originator) if originator else ''

            # Primary Companies - extract from Company or direct value
            primary_companies = drug_record.get('CompaniesPrimary', {})
            company_list = []
            if isinstance(primary_companies, dict):
                companies = as_list(primary_companies.get('Company', []))
                for comp in companies:
                    if isinstance(comp, dict):
                        # Try extracting direct '$' value first
                        comp_value = extract_value(comp)
                        if comp_value:
                            company_list.append(comp_value)
                    elif isinstance(comp, str):
                        company_list.append(comp)
            row['Primary Companies'] = '; '.join([c for c in company_list if c]) if company_list else ''

            # Primary Indications - extract from list with '@id' and '$' keys
            indications_primary = drug_record.get('IndicationsPrimary', {})
            indication_list = []
            if isinstance(indications_primary, dict):
                indications = as_list(indications_primary.get('Indication', []))
                for ind in indications[:10]:  # Increase to first 10
                    if isinstance(ind, dict):
                        # Extract '$' value directly
                        ind_value = extract_value(ind)
                        if ind_value:
                            indication_list.append(ind_value)
                    elif isinstance(ind, str):
                        indication_list.append(ind)
            row['Primary Indications'] = '; '.join([i for i in indication_list if i]) if indication_list else ''

            # Primary Actions (mechanisms) - extract from list with '@id' and '$' keys
            actions_primary = drug_record.get('ActionsPrimary', {})
            action_list = []
            if isinstance(actions_primary, dict):
                actions = as_list(actions_primary.get('Action', []))
                for action in actions[:5]:  # Increase to first 5
                    if isinstance(action, dict):
                        # Extract '$' value directly
                        action_value = extract_value(action)
                        if action_value:
                            action_list.append(action_value)
                    elif isinstance(action, str):
                        action_list.append(action)
            row['Mechanisms of Action'] = '; '.join([a for a in action_list if a]) if action_list else ''

            # Therapy Areas - can be list of strings or dicts
            therapy_areas = drug_record.get('TherapyAreas', {})
            area_list = []
            if isinstance(therapy_areas, dict):
                areas = as_list(therapy_areas.get('TherapyArea', []))
                for area in areas:
                    if isinstance(area, dict):
                        area_value = extract_value(area)
                        if area_value:
                            area_list.append(area_value)
                    elif isinstance(area, str):
                        area_list.append(area)
            row['Therapy Areas'] = '; '.join([a for a in area_list if a]) if area_list else ''

            # Regulatory Designations (like Orphan Drug, Fast Track)
            reg_designations = drug_record.get('RegulatoryDesignations', {})
            designation_list = []
            if isinstance(reg_designations, dict):
                designations = as_list(reg_designations.get('RegulatoryDesignation', []))
                for desig in designations:
                    if isinstance(desig, dict):
                        desig_value = extract_value(desig)
                        if desig_value:
                            designation_list.append(desig_value)
                    elif isinstance(desig, str):
                        designation_list.append(desig)
            row['Regulatory Designations'] = '; '.join([d for d in designation_list if d]) if designation_list else ''

            # Technologies (drug modality)
            technologies = drug_record.get('Technologies', {})
            tech_list = []
            if isinstance(technologies, dict):
                techs = as_list(technologies.get('Technology', []))
                for tech in techs[:3]:  # First 3
                    if isinstance(tech, dict):
                        tech_value = extract_value(tech)
                        if tech_value:
                            tech_list.append(tech_value)
                    elif isinstance(tech, str):
                        tech_list.append(tech)
            row['Technologies'] = '; '.join([t for t in tech_list if t]) if tech_list else ''

            # Drug Synonyms - extract from Name list with Value keys
            synonyms = drug_record.get('DrugSynonyms', {})
            synonym_list = []
            if isinstance(synonyms, dict):
                names = as_list(synonyms.get('Name', []))
                for name_entry in names[:10]:  # First 10
                    if isinstance(name_entry, dict):
                        value = name_entry.get('Value', '')
                        if value:
                            synonym_list.append(str(value))
                    elif name_entry:
                        synonym_list.append(str(name_entry))
            row['Synonyms'] = '; '.join([s for s in synonym_list if s]) if synonym_list else ''

            # Structure SMILES (chemical structure)
            row['Structure (SMILES)'] = drug_record.get('StructureSmiles', '')

            rows.append(row)

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=[
        'Drug ID', 'Drug Name', 'Source Drug Name (Targets API)', 'Source Display Name',
        'Highest Phase', 'Originator', 'Primary Companies', 'Primary Indications',
        'Mechanisms of Action', 'Therapy Areas', 'Regulatory Designations',
        'Technologies', 'Synonyms', 'Structure (SMILES)'
    ])


def flatten_trials_to_df(data: Dict) -> pd.DataFrame:
    """Flatten clinical trials data to DataFrame."""
    rows = []
    trials = as_list(data.get('Trial', []))

    for trial in trials:
        row = {
            'Trial ID': trial.get('@id', ''),
            'Title (Display)': trial.get('TitleDisplay', ''),
            'Title (Official)': trial.get('TitleOfficial', ''),
        }

        # Indications - extract from nested structure
        indications_data = trial.get('Indications', {})
        indication_list = []
        if isinstance(indications_data, dict):
            indications = as_list(indications_data.get('Indication', []))
            for ind in indications:
                if isinstance(ind, dict):
                    ind_value = extract_value(ind)
                    if ind_value:
                        indication_list.append(ind_value)
                elif isinstance(ind, str):
                    indication_list.append(ind)
        row['Indications'] = '; '.join([i for i in indication_list if i]) if indication_list else ''

        # Note: ProtocolAndOutcomes appears to be empty in the data
        protocol = trial.get('ProtocolAndOutcomes', '')
        row['Protocol Info'] = protocol if protocol else ''

        rows.append(row)

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=['Trial ID', 'Title (Display)', 'Title (Official)', 'Indications', 'Protocol Info'])


def flatten_biomarker_to_df(data: Dict) -> pd.DataFrame:
    """Flatten biomarker data to DataFrame."""
    rows = []
    biomarker_uses = as_list(data.get('BiomarkerUseResult', []))

    for use_entry in biomarker_uses:
        # Get biomarker details
        biomarker = use_entry.get('Biomarker', {})
        biomarker_name = biomarker.get('@mainName', '') or extract_value(biomarker.get('biomarkerName'))
        biomarker_id = biomarker.get('@id', '')

        # Get biomarker use details
        use = use_entry.get('BiomarkerUse', {})

        # Extract indication
        indication_data = use.get('BiomarkerUseIndication', {})
        indication = extract_value(indication_data.get('Indication', {}))

        # Extract drug
        drugs_studied = use.get('BiomarkerUseDrugsStudied', {})
        drug = ''
        if isinstance(drugs_studied, dict):
            drug_list = as_list(drugs_studied.get('Drug', []))
            drug_names = [extract_value(d.get('DrugName', d)) if isinstance(d, dict) else str(d) for d in drug_list]
            drug = '; '.join([d for d in drug_names if d])

        # Extract role
        role_data = use.get('BiomarkerUseRole', {})
        role = extract_value(role_data)

        # Extract indication type
        indication_type_data = use.get('BiomarkerUseIndicationType', {})
        indication_type = extract_value(indication_type_data)

        # Extract validity
        validity_data = use.get('BiomarkerUseValidity', {})
        validity = extract_value(validity_data)

        row = {
            'Biomarker Use ID': use_entry.get('@id', ''),
            'Biomarker ID': biomarker_id,
            'Biomarker Name': biomarker_name,
            'Role': role,
            'Indication': indication,
            'Indication Type': indication_type,
            'Drugs Studied': drug,
            'Validity': validity,
        }

        rows.append(row)

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=['Biomarker Use ID', 'Biomarker ID', 'Biomarker Name', 'Role', 'Indication', 'Indication Type', 'Drugs Studied', 'Validity'])


def flatten_association_to_df(data: Dict) -> pd.DataFrame:
    """Flatten association (condition-gene) data to DataFrame."""
    rows = []
    target = data.get('Target', {})

    # Condition Gene Associations
    cond_gene_assoc = target.get('ConditionGeneAssociations', {})
    conditions = as_list(cond_gene_assoc.get('Condition', []))

    for cond in conditions:
        row = {
            'Condition Name': extract_value(cond.get('ConditionName')),
            'Condition ID': extract_value(cond.get('ConditionId')),
            'Association Type': extract_value(cond.get('AssociationType')),
            'Evidence': extract_value(cond.get('Evidence')),
            'Source': extract_value(cond.get('Source')),
            'PubMed ID': extract_value(cond.get('PubmedId')),
        }
        rows.append(row)

    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=['Condition Name', 'Condition ID', 'Association Type', 'Evidence', 'Source', 'PubMed ID'])


def generate_excel(gene: str, data: Dict, fields: List[str], output_file: Path):
    """Generate Excel file with one sheet per field."""
    print(f"Generating Excel file: {output_file}")

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for field in fields:
            if field not in data:
                continue

            df = None

            if field == 'annotation':
                df = flatten_annotation_to_df(data[field])
            elif field == 'interaction':
                df = flatten_interaction_to_df(data[field])
            elif field == 'drug':
                # Create sheet for all basic drugs from Drug list
                basic_drugs_df = flatten_basic_drugs_to_df(data[field])
                if not basic_drugs_df.empty:
                    basic_drugs_df.to_excel(writer, sheet_name='Drugs_All', index=False)

                # Create sheet for comprehensive drug records
                df = flatten_drug_to_df(data[field])
                if not df.empty:
                    df.to_excel(writer, sheet_name='Drugs_Comprehensive', index=False)
                    df = None  # Prevent duplicate export below

                # Also add trials sheet
                trials_df = flatten_trials_to_df(data[field])
                if not trials_df.empty:
                    trials_df.to_excel(writer, sheet_name='Trials', index=False)
            elif field == 'biomarker':
                df = flatten_biomarker_to_df(data[field])
            elif field == 'association':
                df = flatten_association_to_df(data[field])

            if df is not None and not df.empty:
                # Capitalize sheet name
                sheet_name = field.capitalize()
                df.to_excel(writer, sheet_name=sheet_name, index=False)

                # Auto-adjust column widths
                worksheet = writer.sheets[sheet_name]
                for idx, col in enumerate(df.columns):
                    max_length = max(
                        df[col].astype(str).apply(len).max(),
                        len(str(col))
                    )
                    worksheet.column_dimensions[chr(65 + idx)].width = min(max_length + 2, 50)


# Global data loss tracker
data_loss_tracker = {
    'failed_drug_searches': [],
    'no_results_drug_names': [],
    'failed_drug_chunks': 0,
    'failed_trial_chunks': 0,
    'failed_comprehensive_drug_chunks': 0,
}


def reset_data_loss_tracker():
    """Reset the global data loss tracker."""
    global data_loss_tracker
    data_loss_tracker = {
        'failed_drug_searches': [],
        'no_results_drug_names': [],
        'failed_drug_chunks': 0,
        'failed_trial_chunks': 0,
        'failed_comprehensive_drug_chunks': 0,
    }


def report_data_loss() -> bool:
    """
    Report data loss and return True if there was any data loss.

    Returns:
        bool: True if there was data loss, False otherwise
    """
    print("\n" + "=" * 60)
    print("DATA LOSS REPORT")
    print("=" * 60)

    has_loss = False

    # Failed drug searches
    if data_loss_tracker['failed_drug_searches']:
        has_loss = True
        print(f"\n⚠️  Failed Drug Searches (timeout/error): {len(data_loss_tracker['failed_drug_searches'])}")
        for name in data_loss_tracker['failed_drug_searches'][:10]:
            print(f"    - {name}")
        if len(data_loss_tracker['failed_drug_searches']) > 10:
            print(f"    ... and {len(data_loss_tracker['failed_drug_searches']) - 10} more")

    # Drug names with no results
    if data_loss_tracker['no_results_drug_names']:
        print(f"\nℹ️  Drug Names with No Results in Investigational Drugs API: {len(data_loss_tracker['no_results_drug_names'])}")
        print("   (These drugs exist in Targets API but not in Investigational Drugs database)")
        if len(data_loss_tracker['no_results_drug_names']) <= 10:
            for name in data_loss_tracker['no_results_drug_names']:
                print(f"    - {name}")
        else:
            print(f"    (Too many to display - {len(data_loss_tracker['no_results_drug_names'])} total)")

    # Failed chunks
    if data_loss_tracker['failed_drug_chunks'] > 0:
        has_loss = True
        print(f"\n⚠️  Failed Drug Record Chunks: {data_loss_tracker['failed_drug_chunks']}")

    if data_loss_tracker['failed_trial_chunks'] > 0:
        has_loss = True
        print(f"⚠️  Failed Trial Record Chunks: {data_loss_tracker['failed_trial_chunks']}")

    if data_loss_tracker['failed_comprehensive_drug_chunks'] > 0:
        has_loss = True
        print(f"⚠️  Failed Comprehensive Drug Chunks: {data_loss_tracker['failed_comprehensive_drug_chunks']}")

    if not has_loss:
        print("\n✅ No data loss detected - all requests succeeded!")

    print("=" * 60)
    return has_loss


def ask_retry(retry_count: int = 0, max_retries: int = 3) -> bool:
    """
    Automatically retry failed requests up to max_retries times.

    Args:
        retry_count: Current retry attempt number (0-indexed)
        max_retries: Maximum number of retry attempts (default: 3)

    Returns:
        bool: True if should retry, False if max retries exceeded
    """
    if not data_loss_tracker['failed_drug_searches']:
        return False

    if retry_count >= max_retries:
        print(f"\n❌ Max retries ({max_retries}) reached. Skipping further retry attempts.")
        return False

    print(f"\n🔄 Auto-retry {retry_count + 1}/{max_retries} for {len(data_loss_tracker['failed_drug_searches'])} failed drug searches...")
    return True


def retry_failed_drug_searches(headers: Dict) -> Dict:
    """
    Retry failed drug searches with more aggressive settings.

    Args:
        headers: HTTP headers for authentication

    Returns:
        Dict mapping drug_id -> (source_drug_name, source_display_name)
    """
    if not data_loss_tracker['failed_drug_searches']:
        return {}

    print(f"\nRetrying {len(data_loss_tracker['failed_drug_searches'])} failed drug searches...")
    print("Using sequential requests with extended timeout (240s)...")

    all_drug_ids = {}
    successful_retries = []
    still_failed = []

    for drug_namemain in data_loss_tracker['failed_drug_searches']:
        try:
            # Try with longer timeout and sequential (no concurrency stress)
            search_url = requests.Request(
                'GET',
                f"{CORTELLIS_DRUGS_URL}/drug/search",
                params={
                    'query': f'drugNamesAll:{drug_namemain}',
                    'hits': 2000,
                    'returnFilterCount': 5000,
                    'fmt': 'json'
                }
            ).prepare().url

            print(f"  Retrying: {drug_namemain}...", end='', flush=True)
            response = requests.get(search_url, headers=headers, timeout=240)

            if response.status_code == 200:
                resp_json = response.json()
                drug_results_output = resp_json.get('drugResultsOutput', {})
                search_results_data = drug_results_output.get('SearchResults', {})
                drug_result_data = search_results_data.get('DrugResult') or search_results_data.get('Drug', [])
                results = as_list(drug_result_data)

                if results:
                    for result in results:
                        if isinstance(result, dict):
                            drug_id = result.get('@id')
                            if drug_id:
                                all_drug_ids[drug_id] = (drug_namemain, drug_namemain)
                    successful_retries.append(drug_namemain)
                    print(f" ✅ Success ({len(results)} IDs)")
                else:
                    still_failed.append(drug_namemain)
                    print(" ⚠️  No results")
            else:
                still_failed.append(drug_namemain)
                print(f" ❌ Failed (HTTP {response.status_code})")
        except Exception as e:
            still_failed.append(drug_namemain)
            print(f" ❌ Error: {str(e)[:50]}")

    # Update tracker
    data_loss_tracker['failed_drug_searches'] = still_failed

    print(f"\nRetry Summary:")
    print(f"  ✅ Successful: {len(successful_retries)}")
    print(f"  ❌ Still Failed: {len(still_failed)}")

    return all_drug_ids


def query_gene(gene: str, fields: List[str]) -> Dict:
    """Query Cortellis for a single gene with specified fields."""
    print(f"\nQuerying gene: {gene}")
    print("=" * 60)

    # Reset data loss tracker for new query
    reset_data_loss_tracker()

    result = {}

    # Always get annotation first to get target_id
    print("Fetching annotation...")
    annotation, target_id = target_summaryrecord(gene)
    result['annotation'] = annotation

    # Query requested fields
    for field in fields:
        if field == 'annotation':
            continue  # Already fetched

        print(f"Fetching {field}...")

        if field == 'interaction':
            result['interaction'] = target_interaction(target_id)
        elif field == 'drug':
            # Pass annotation to include related targets
            result['drug'] = target_drug(target_id, annotation)
        elif field == 'biomarker':
            result['biomarker'] = target_biomarker(gene)
        elif field == 'association':
            result['association'] = target_association(target_id)

    return result


def generate_summary(gene: str, data: Dict, fields: List[str]) -> str:
    """Generate markdown summary of findings."""
    summary = []
    summary.append(f"# Cortellis Data Summary: {gene}\n")
    summary.append(f"*Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*\n")
    summary.append("---\n")

    for field in fields:
        if field in data:
            if field == 'annotation':
                summary.append(summarize_annotation(data[field]))
            elif field == 'interaction':
                summary.append(summarize_interaction(data[field]))
            elif field == 'drug':
                summary.append(summarize_drug(data[field]))
            elif field == 'biomarker':
                summary.append(summarize_biomarker(data[field]))
            elif field == 'association':
                summary.append(summarize_association(data[field]))
            summary.append("\n---\n")

    return "\n".join(summary)


def main():
    parser = argparse.ArgumentParser(
        description='Query Cortellis API for gene annotations and related data.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Query single gene with default field (drug)
  python cortellis_gene_query.py BRCA1

  # Query multiple genes with default field
  python cortellis_gene_query.py BRCA1 TP53 EGFR

  # Query with specific fields
  python cortellis_gene_query.py BRCA1 --fields drug interaction biomarker

  # Query with all fields and generate Excel
  python cortellis_gene_query.py BRCA1 --all --excel

  # Save to custom output directory
  python cortellis_gene_query.py BRCA1 --output-dir ./results --excel
        """
    )

    parser.add_argument('genes', nargs='+', help='Gene symbol(s) to query')
    parser.add_argument(
        '--fields',
        nargs='+',
        choices=['drug', 'biomarker', 'interaction', 'association'],
        help='Fields to query (annotation is always included). Default: drug'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Query all available fields'
    )
    parser.add_argument(
        '--output-dir',
        default=None,
        help='Output directory for JSON files (default: current working directory)'
    )
    parser.add_argument(
        '--no-summary',
        action='store_true',
        help='Skip printing summary to console'
    )
    parser.add_argument(
        '--excel',
        action='store_true',
        help='Generate Excel file with data tables (one sheet per field)'
    )

    args = parser.parse_args()

    # Determine which fields to query
    if args.all:
        fields = ['annotation', 'drug', 'biomarker', 'interaction', 'association']
    elif args.fields:
        fields = ['annotation'] + args.fields  # Always include annotation
    else:
        # Default: drug only
        fields = ['annotation', 'drug']

    # Determine output directory - use current working directory by default
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = Path.cwd()

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nQuerying {len(args.genes)} gene(s) with fields: {', '.join(fields)}")
    print(f"Output directory: {output_dir.absolute()}\n")

    # Process each gene
    for gene in args.genes:
        try:
            # Query API
            data = query_gene(gene, fields)

            # Save JSON
            output_file = output_dir / f"{gene}_cortellis_data.json"
            with open(output_file, 'w') as f:
                json.dump(data, f, indent=2)
            print(f"\nSaved: {output_file}")

            # Generate and print summary
            if not args.no_summary:
                summary = generate_summary(gene, data, fields)
                print(summary)

                # Save summary as markdown file
                summary_file = output_dir / f"{gene}_summary.md"
                with open(summary_file, 'w') as f:
                    f.write(summary)
                print(f"Summary saved: {summary_file}")

            # Generate Excel if requested
            if args.excel:
                excel_file = output_dir / f"{gene}_cortellis_data.xlsx"
                generate_excel(gene, data, fields, excel_file)
                print(f"Excel saved: {excel_file}")

            # Report data loss
            has_loss = report_data_loss()

            # Automatic retry loop (up to 3 attempts)
            retry_count = 0
            while has_loss and ask_retry(retry_count, max_retries=3):
                # Retry failed drug searches
                retry_drug_ids = retry_failed_drug_searches(DRUGS_HEADER)

                if retry_drug_ids:
                    print(f"\n  Fetching {len(retry_drug_ids)} recovered comprehensive drug records...")

                    # Bulk fetch comprehensive records for recovered drug IDs
                    drug_ids_list = list(retry_drug_ids.keys())
                    chunk_size = 30
                    chunk_urls = []
                    for i in range(0, len(drug_ids_list), chunk_size):
                        chunk = drug_ids_list[i:i + chunk_size]
                        drug_id_str = ",".join(chunk)
                        drugs_url = requests.Request(
                            'GET',
                            f"{CORTELLIS_DRUGS_URL}/drugs",
                            params={'idList': drug_id_str, 'fmt': 'json'}
                        ).prepare().url
                        chunk_urls.append(drugs_url)

                    # Fetch recovered drugs
                    results = fetch_with_adaptive_concurrency(chunk_urls, DRUGS_HEADER, initial_concurrent=25, fallback_concurrent=20)

                    recovered_records = {}
                    for resp in results:
                        if resp:
                            drug_records_output = resp.get('drugRecordsOutput', {})
                            drugs = as_list(drug_records_output.get('Drug', []))
                            for drug_record in drugs:
                                if isinstance(drug_record, dict):
                                    record_id = drug_record.get('@id')
                                    if record_id:
                                        source_names = retry_drug_ids.get(record_id, ('', ''))
                                        drug_record['_source_drug_name'] = source_names[0]
                                        drug_record['_source_display_name'] = source_names[1]
                                        recovered_records[record_id] = drug_record

                    print(f"  Recovered {len(recovered_records)} additional comprehensive drug records")

                    # Merge recovered records into existing data
                    if 'drug' in data and 'DrugRecord' in data['drug']:
                        data['drug']['DrugRecord'].update(recovered_records)
                        print(f"  Total comprehensive drug records: {len(data['drug']['DrugRecord'])}")

                        # Re-save JSON
                        print(f"\n  Updating saved files with recovered data...")
                        with open(output_file, 'w') as f:
                            json.dump(data, f, indent=2)
                        print(f"  Updated: {output_file}")

                        # Re-generate Excel if requested
                        if args.excel:
                            generate_excel(gene, data, fields, excel_file)
                            print(f"  Updated: {excel_file}")

                # Check if there are still failures
                has_loss = report_data_loss()
                retry_count += 1

                # If no more failures, break the loop
                if not has_loss or not data_loss_tracker['failed_drug_searches']:
                    break

        except Exception as e:
            print(f"\nERROR processing gene {gene}: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 60)
    print("Query completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
