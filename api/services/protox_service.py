from __future__ import annotations

import os
import re
from typing import Any

import httpx
from rdkit import Chem
from rdkit.Chem import rdDepictor


from api.models import ProToxPrediction, Endpoints


PROTOX_BASE_URL = "https://tox.charite.de/protox3"
PROTOX_TIMEOUT_SECONDS = float(os.getenv("PROTOX_TIMEOUT_SECONDS", "90"))
PROTOX_USE_STUB = os.getenv("PROTOX_USE_STUB", "false").lower() == "true"

_CLASS_LABELS = {
    1: "Fatal",
    2: "Fatal",
    3: "Toxic",
    4: "Harmful",
    5: "May Be Harmful",
    6: "Non-toxic",
}


def _build_stub(smiles: str) -> ProToxPrediction:
    return ProToxPrediction(
        smiles=smiles,
        ld50=200,
        toxicity_class=3,
        toxicity_class_label=_CLASS_LABELS[3],
        endpoints=Endpoints(
            hepatotoxicity=True,
            nephrotoxicity=False,
            immunotoxicity=False,
            mutagenicity=False,
            cytotoxicity=True,
        ),
    )


def _smiles_to_molblock(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError("Invalid SMILES string")
    mol = Chem.AddHs(mol)
    rdDepictor.Compute2DCoords(mol)
    return Chem.MolToMolBlock(mol)


def _extract_ld50(html: str) -> int:
    match = re.search(r"Predicted LD50:\s*([0-9]+)\s*mg/kg", html, re.IGNORECASE)
    if not match:
        raise RuntimeError("Could not extract LD50 from ProTox response")
    return int(match.group(1))


def _extract_toxicity_class(html: str) -> int:
    match = re.search(r"Predicted Toxicity Class:\s*([1-6])", html, re.IGNORECASE)
    if not match:
        raise RuntimeError("Could not extract toxicity class from ProTox response")
    return int(match.group(1))


def _extract_server_id(html: str) -> str | None:
    match = re.search(r"var\s+server_id\s*=\s*'([^']+)'", html)
    if match:
        return match.group(1)
    return None


def _parse_model_predictions(raw: dict[str, Any]) -> Endpoints:
    def is_active(key: str) -> bool:
        item = raw.get(key)
        if not item:
            return False
        pred = str(item.get("Prediction", "0")).strip()
        return pred in {"1", "1.0"}

    return Endpoints(
        hepatotoxicity=is_active("dili"),
        nephrotoxicity=is_active("nephro"),
        immunotoxicity=is_active("immuno"),
        mutagenicity=is_active("mutagen"),
        cytotoxicity=is_active("cyto"),
    )


async def _call_compound_search_similarity(
    client: httpx.AsyncClient,
    smiles: str,
    molblock: str,
) -> tuple[int, int, str | None]:
    response = await client.post(
        f"{PROTOX_BASE_URL}/index.php?site=compound_search_similarity",
        data={
            "smilesString": molblock,
            "defaultName": "User defined",
            "smiles": smiles,
            "pubchem_name": "",
        },
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
        },
    )
    response.raise_for_status()
    html = response.text

    ld50 = _extract_ld50(html)
    tox_class = _extract_toxicity_class(html)
    server_id = _extract_server_id(html)

    return ld50, tox_class, server_id


async def _call_run_models(
    client: httpx.AsyncClient,
    molblock: str,
    server_id: str,
) -> Endpoints:
    model_string = "dili nephro immuno mutagen cyto"

    response = await client.post(
        f"{PROTOX_BASE_URL}/src/run_models.php",
        data={
            "models": model_string,
            "sdfile": "empty",
            "mol": molblock,
            "id": server_id,
        },
        headers={
            "X-Requested-With": "XMLHttpRequest",
            "Content-Type": "application/x-www-form-urlencoded",
        },
    )
    response.raise_for_status()

    raw = response.json()
    return _parse_model_predictions(raw)


async def get_protox_prediction(smiles: str) -> ProToxPrediction:
    if PROTOX_USE_STUB:
        return _build_stub(smiles)

    molblock = _smiles_to_molblock(smiles)

    async with httpx.AsyncClient(
        timeout=PROTOX_TIMEOUT_SECONDS,
        follow_redirects=True,
    ) as client:
        try:
            ld50, tox_class, server_id = await _call_compound_search_similarity(
                client=client,
                smiles=smiles,
                molblock=molblock,
            )

            endpoints = Endpoints()

            if server_id:
                try:
                    endpoints = await _call_run_models(
                        client=client,
                        molblock=molblock,
                        server_id=server_id,
                    )
                except Exception:
                    # keep acute toxicity result even if optional endpoint models fail
                    endpoints = Endpoints()

            return ProToxPrediction(
                smiles=smiles,
                ld50=ld50,
                toxicity_class=tox_class,
                toxicity_class_label=_CLASS_LABELS.get(tox_class, "Unknown"),
                endpoints=endpoints,
            )

        except httpx.HTTPStatusError as e:
            raise RuntimeError(f"ProTox returned HTTP {e.response.status_code}") from e
        except httpx.RequestError as e:
            raise RuntimeError("Could not reach ProTox service") from e
        except Exception as e:
            raise RuntimeError(f"Failed to process ProTox response: {e}") from e