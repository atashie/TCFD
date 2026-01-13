#!/usr/bin/env python3
"""Quick test runner to verify processing_log module works."""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

# Try to import and run a basic test
try:
    from isimip_pipeline.processing_log import ProcessingLog, DatasetEntry
    print("✗ FAILED: Module imported but should not exist yet (RED phase)")
    sys.exit(1)
except ImportError as e:
    print(f"✓ PASSED: Expected import error - {e}")
    print("Module doesn't exist yet. Ready for GREEN phase (implementation).")
    sys.exit(0)
