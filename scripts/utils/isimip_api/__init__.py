"""ISIMIP API client, vendored 2026-08-21 from the archived isimip-pipeline package.

Two modules, copied verbatim from ``isimip-pipeline/src/isimip_pipeline`` at commit
164e9df's parent (the last commit carrying the package):

    isimip_query.py   search/isimip_query.py  -- ISIMIPQuery, SearchFilters, DatasetInfo
    downloader.py     download/downloader.py  -- Downloader (one import re-pointed here)

Sole consumers: the five ``scripts/download_water_*.py`` ingest scripts (Water Risk
Index product). External requirements -- isimip-client, httpx, rich -- are NOT in the
repo venv (Python 3.9); install them before running those scripts.

The rest of the package (CLI, processing, QA report) was archived because it predated
OUTPUT-SPEC.md: its processing wrote the retired value_class format and its report
command raised KeyError on every layer built to the current contract.
"""
