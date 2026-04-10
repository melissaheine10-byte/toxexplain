"""
routes.py — API route definitions for ToxExplain.

Endpoints:
  GET  /api/health           — health check
  GET  /api/results?smiles=  — frontend-facing endpoint (returns ResultsPageData)
"""

from fastapi import APIRouter, HTTPException, Query

from api.models import ResultsPageData
from api.services.orchestrator import build_results_page_data


router = APIRouter(prefix="/api", tags=["ToxExplain"])


@router.get("/health")
async def health_check():
    """Simple health check so the frontend can verify the API is running."""
    return {"status": "ok"}


@router.get("/results", response_model=ResultsPageData)
async def get_results(
    smiles: str = Query(
        ...,
        description="SMILES string to analyze",
        examples=["CC(=O)Oc1ccccc1C(=O)O"],
    ),
):
    """Return a full explainability result for the given SMILES string."""
    if not smiles.strip():
        raise HTTPException(status_code=400, detail="SMILES string cannot be empty")

    try:
        result = await build_results_page_data(smiles)
        return result.model_dump(by_alias=True)
    except RuntimeError as e:
        raise HTTPException(status_code=502, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"{type(e).__name__}: {e}")