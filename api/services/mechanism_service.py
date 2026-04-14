"""
mechanism_service.py — Local mechanism data lookup.

Loads curated mechanism records from api/data/mechanism_records.json
and filters them by SMILES match. Uses smiles_index.json for fast lookup
and falls back to canonical SMILES matching via RDKit if available.
"""

import json
from pathlib import Path
from api.models import MechanismRecord, EvidenceItem

_DATA_DIR = Path(__file__).resolve().parent.parent / "data"
_DATA_FILE = _DATA_DIR / "mechanism_records.json"
_INDEX_FILE = _DATA_DIR / "smiles_index.json"

# Module-level cache
_raw_data: list[dict] | None = None
_smiles_index: dict | None = None


def _load_data() -> list[dict]:
    """Read the mechanism records JSON file (cached)."""
    global _raw_data
    if _raw_data is None:
        with open(_DATA_FILE, "r") as f:
            _raw_data = json.load(f)
    return _raw_data


def _load_index() -> dict:
    """Read the SMILES index file (cached). Returns empty dict if not found."""
    global _smiles_index
    if _smiles_index is None:
        if _INDEX_FILE.exists():
            with open(_INDEX_FILE, "r") as f:
                _smiles_index = json.load(f)
        else:
            _smiles_index = {"by_smiles": {}, "by_drug_id": {}, "by_drug_name": {}}
    return _smiles_index


def _find_matching_entries(smiles: str, raw_data: list[dict]) -> list[dict]:
    """Find mechanism records matching a SMILES string."""

    # 1. Try exact SMILES match against records
    exact = [e for e in raw_data if e.get("smiles") == smiles]
    if exact:
        return exact

    # 2. Try SMILES index lookup
    index = _load_index()
    drug_id = index.get("by_smiles", {}).get(smiles)
    if drug_id:
        by_id = [e for e in raw_data if e.get("drug_id") == drug_id]
        if by_id:
            return by_id

    # 3. Try canonical SMILES matching via RDKit
    try:
        from rdkit import Chem

        query_mol = Chem.MolFromSmiles(smiles)
        if query_mol:
            query_canonical = Chem.MolToSmiles(query_mol)

            # Check against all records
            for entry in raw_data:
                entry_smiles = entry.get("smiles", "")
                if not entry_smiles:
                    continue
                entry_mol = Chem.MolFromSmiles(entry_smiles)
                if entry_mol and Chem.MolToSmiles(entry_mol) == query_canonical:
                    return [entry]

            # Check against index
            for stored_smiles, did in index.get("by_smiles", {}).items():
                stored_mol = Chem.MolFromSmiles(stored_smiles)
                if stored_mol and Chem.MolToSmiles(stored_mol) == query_canonical:
                    by_id = [e for e in raw_data if e.get("drug_id") == did]
                    if by_id:
                        return by_id
    except ImportError:
        pass  # RDKit not installed

    return []


async def lookup_mechanisms(
    smiles: str,
) -> tuple[list[MechanismRecord], list[EvidenceItem]]:
    """
    Look up curated mechanism records that match the given SMILES string.

    Tries in order:
      1. Exact SMILES match against records
      2. SMILES index lookup (by drug_id)
      3. Canonical SMILES match via RDKit (if installed)

    Returns empty lists if no match is found.
    """
    raw_data = _load_data()
    matched = _find_matching_entries(smiles, raw_data)

    mechanism_records: list[MechanismRecord] = []
    evidence_items: list[EvidenceItem] = []

    for entry in matched:
        # Pop evidence_items from a copy so we don't mutate the cache
        entry_copy = dict(entry)
        raw_evidence = entry_copy.pop("evidence_items", [])

        # Remove drug_id if present (not in the MechanismRecord model)
        entry_copy.pop("drug_id", None)

        mechanism_records.append(MechanismRecord(**entry_copy))
        evidence_items.extend(EvidenceItem(**ev) for ev in raw_evidence)

    return mechanism_records, evidence_items


async def lookup_by_name(
    drug_name: str,
) -> tuple[list[MechanismRecord], list[EvidenceItem]]:
    """Look up by drug name (case-insensitive)."""
    index = _load_index()
    drug_id = index.get("by_drug_name", {}).get(drug_name.lower())
    if not drug_id:
        return [], []

    drug_info = index.get("by_drug_id", {}).get(drug_id, {})
    smiles = drug_info.get("smiles", "")
    if smiles:
        return await lookup_mechanisms(smiles)
    return [], []


def get_all_drug_names() -> list[dict]:
    """Return list of all available drugs (useful for frontend dropdown)."""
    index = _load_index()
    drugs = []
    for drug_id, info in index.get("by_drug_id", {}).items():
        drugs.append(
            {
                "drug_id": drug_id,
                "drug_name": info["drug_name"],
                "smiles": info["smiles"],
            }
        )
    return sorted(drugs, key=lambda x: x["drug_name"])