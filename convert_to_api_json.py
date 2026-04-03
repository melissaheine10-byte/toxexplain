"""
convert_to_api_json.py — Convert LLM extraction Excel into API-compatible JSON.

Matches the EXISTING format in api/data/mechanism_records.json, where
evidence_items are nested inside each drug record.

Reads:
  - llm_extraction_output.xlsx (LLM pipeline output)
  - Withdrawn_Drugs_Mechanistic_Curation_Template.xlsx (for withdrawal reasons)
  - drugbank.xml (for SMILES lookup)
  - api/data/mechanism_records.json (existing records to preserve)

Writes:
  - api/data/mechanism_records.json (merged: existing + new)
  - api/data/smiles_index.json (lookup table for SMILES matching)

USAGE:
  python convert_to_api_json.py \
    --llm llm_extraction_output.xlsx \
    --gt Withdrawn_Drugs_Mechanistic_Curation_Template.xlsx \
    --xml drugbank.xml \
    --existing api/data/mechanism_records.json \
    --outdir api/data
"""

import argparse
import json
import os
import pandas as pd
from lxml import etree

NS = {"db": "http://www.drugbank.ca"}


def extract_smiles_from_xml(xml_path, drug_ids):
    print(f"Extracting SMILES from {xml_path}...")
    smiles_map = {}
    drug_ids = set(drug_ids)
    context = etree.iterparse(xml_path, events=("end",), tag=f"{{{NS['db']}}}drug")
    for event, elem in context:
        db_id_elem = elem.find('db:drugbank-id[@primary="true"]', NS)
        if db_id_elem is None or db_id_elem.text not in drug_ids:
            elem.clear()
            continue
        db_id = db_id_elem.text
        if db_id in smiles_map:
            elem.clear()
            continue
        for prop in elem.findall("db:calculated-properties/db:property", NS):
            kind = prop.find("db:kind", NS)
            val = prop.find("db:value", NS)
            if kind is not None and kind.text == "SMILES" and val is not None:
                smiles_map[db_id] = val.text
                break
        if len(smiles_map) == len(drug_ids):
            elem.clear()
            break
        elem.clear()
    print(f"  Found SMILES for {len(smiles_map)}/{len(drug_ids)} drugs")
    return smiles_map


def build_records(llm_sheets, gt_sheets, smiles_map):
    records = []
    overview = llm_sheets["Overview"]
    gt_overview = gt_sheets["Overview"]
    withdrawal_reasons = dict(zip(gt_overview["drug_id"], gt_overview["withdrawal_reason"]))

    for _, row in overview.iterrows():
        drug_id = row["drug_id"]
        drug_name = row["drug_name"]

        targets_df = llm_sheets["Targets"]
        drug_targets = targets_df[targets_df["drug_id"] == drug_id]
        primary_targets = drug_targets[drug_targets["role"] == "primary"]["target_name"].tolist()
        off_targets = drug_targets[drug_targets["role"].isin(["off_target", "secondary"])]["target_name"].tolist()

        mech_df = llm_sheets["Mechanisms"]
        drug_mechs = mech_df[mech_df["drug_id"] == drug_id]
        pathways = drug_mechs["Mechanism_name"].dropna().tolist()

        organs = set()
        for sheet in ["Targets", "Mechanisms", "Adverse Events", "Risk Modifiers"]:
            df = llm_sheets[sheet]
            drug_rows = df[df["drug_id"] == drug_id]
            if "organ_system" in drug_rows.columns:
                organs.update(drug_rows["organ_system"].dropna().tolist())
        organs.discard("other")
        organs.discard("unknown")

        refs = set()
        for sheet in ["Targets", "Mechanisms", "Adverse Events", "Risk Modifiers"]:
            df = llm_sheets[sheet]
            drug_rows = df[df["drug_id"] == drug_id]
            if "source_id" in drug_rows.columns:
                refs.update(drug_rows["source_id"].dropna().tolist())

        evidence_items = []
        sheet_origins = {
            "Targets": "target annotation",
            "Mechanisms": "mechanism extraction",
            "Adverse Events": "adverse event report",
            "Risk Modifiers": "risk modifier analysis",
        }
        for sheet_name, origin in sheet_origins.items():
            df = llm_sheets[sheet_name]
            drug_rows = df[df["drug_id"] == drug_id]
            for _, r in drug_rows.iterrows():
                if sheet_name == "Targets":
                    label = f"{r.get('target_name', '')} — {r.get('action', '')} ({r.get('role', '')})"
                elif sheet_name == "Mechanisms":
                    label = str(r.get("Mechanism_name", ""))
                elif sheet_name == "Adverse Events":
                    label = f"{r.get('name', '')} — {r.get('organ_system', '')} ({r.get('severity', '')})"
                elif sheet_name == "Risk Modifiers":
                    label = f"{r.get('factor', '')} — {r.get('direction', '')} risk"
                if not label.strip() or label.strip() == "—  ()":
                    continue
                confidence = str(r.get("confidence", "medium")).lower()
                if confidence not in ("high", "medium", "low"):
                    confidence = "medium"
                evidence_items.append({
                    "label": label,
                    "source": str(r.get("source_id", "DrugBank")),
                    "origin": origin,
                    "group": "B",
                    "confidence": confidence,
                })

        wr = withdrawal_reasons.get(drug_id, "")
        if pd.isna(wr):
            wr = ""

        records.append({
            "drug_name": drug_name,
            "drug_id": drug_id,
            "smiles": smiles_map.get(drug_id, ""),
            "primary_targets": primary_targets,
            "off_targets": off_targets,
            "pathways": pathways,
            "organ_systems": sorted(organs),
            "withdrawal_reason": str(wr),
            "source": "Group B",
            "references": sorted(refs),
            "evidence_items": evidence_items,
        })

    return records


def build_smiles_index(all_records):
    index = {"by_smiles": {}, "by_drug_id": {}, "by_drug_name": {}}
    for rec in all_records:
        drug_id = rec.get("drug_id", rec["drug_name"].lower().replace(" ", "-"))
        smiles = rec.get("smiles", "")
        name = rec["drug_name"].lower()
        if smiles:
            index["by_smiles"][smiles] = drug_id
        index["by_drug_id"][drug_id] = {"drug_name": rec["drug_name"], "smiles": smiles}
        index["by_drug_name"][name] = drug_id
    return index


def main():
    parser = argparse.ArgumentParser(description="Convert LLM extraction to API JSON")
    parser.add_argument("--llm", required=True, help="Path to LLM extraction Excel")
    parser.add_argument("--gt", required=True, help="Path to ground truth Excel")
    parser.add_argument("--xml", required=True, help="Path to DrugBank XML")
    parser.add_argument("--existing", default=None, help="Path to existing mechanism_records.json to preserve")
    parser.add_argument("--outdir", default="api/data", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading Excel files...")
    llm_sheets = pd.read_excel(args.llm, sheet_name=None)
    gt_sheets = pd.read_excel(args.gt, sheet_name=None)

    drug_ids = llm_sheets["Overview"]["drug_id"].dropna().tolist()
    smiles_map = extract_smiles_from_xml(args.xml, drug_ids)

    print("Building mechanism records from LLM extraction...")
    new_records = build_records(llm_sheets, gt_sheets, smiles_map)
    print(f"  Built {len(new_records)} new records")

    existing_records = []
    if args.existing and os.path.exists(args.existing):
        print(f"Loading existing records from {args.existing}...")
        with open(args.existing) as f:
            existing_records = json.load(f)
        print(f"  Found {len(existing_records)} existing records")
        new_names = {r["drug_name"].lower() for r in new_records}
        kept = [r for r in existing_records if r["drug_name"].lower() not in new_names]
        replaced = len(existing_records) - len(kept)
        if replaced:
            print(f"  Replacing {replaced} existing record(s) with LLM extraction")
        existing_records = kept

    all_records = existing_records + new_records
    print(f"  Total: {len(all_records)} records ({len(existing_records)} existing + {len(new_records)} new)")

    print("Building SMILES index...")
    smiles_index = build_smiles_index(all_records)

    mech_path = os.path.join(args.outdir, "mechanism_records.json")
    with open(mech_path, "w") as f:
        json.dump(all_records, f, indent=2)
    print(f"  Wrote {len(all_records)} records -> {mech_path}")

    idx_path = os.path.join(args.outdir, "smiles_index.json")
    with open(idx_path, "w") as f:
        json.dump(smiles_index, f, indent=2)
    print(f"  Wrote index -> {idx_path}")

    print(f"\nDone. Summary:")
    for rec in all_records:
        n_ev = len(rec.get("evidence_items", []))
        n_pt = len(rec.get("primary_targets", []))
        n_ot = len(rec.get("off_targets", []))
        print(f"  {rec['drug_name']:20s} | {n_pt} primary, {n_ot} off-target | {n_ev} evidence items | source: {rec['source']}")


if __name__ == "__main__":
    main()
