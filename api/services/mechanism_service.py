"""
mechanism_service.py — Local mechanism data lookup.

Loads curated mechanism records from api/data/mechanism_records.json
and filters them by SMILES match.
"""

import json
from pathlib import Path
from api.models import MechanismRecord, EvidenceItem


_DATA_FILE = Path(__file__).resolve().parent.parent / "data" / "mechanism_records.json"


def _load_mechanism_data() -> list[dict]:
    """Read the JSON file and return raw records."""
    with open(_DATA_FILE, "r") as f:
        return json.load(f)


async def lookup_mechanisms(
    smiles: str,
) -> tuple[list[MechanismRecord], list[EvidenceItem]]:
    """
    Look up curated mechanism records that match the given SMILES string.

    Currently does exact SMILES match.

    Returns:
        A tuple of (mechanism_records, evidence_items) extracted from the data.
    """
    raw_data = _load_mechanism_data()

    mechanism_records: list[MechanismRecord] = []
    evidence_items: list[EvidenceItem] = []

    for entry in raw_data:
        if entry.get("smiles") == smiles:
            entry_copy = dict(entry)
            raw_evidence = entry_copy.pop("evidence_items", [])

            mechanism_records.append(MechanismRecord(**entry_copy))
            evidence_items.extend(EvidenceItem(**ev) for ev in raw_evidence)

    return mechanism_records, evidence_items