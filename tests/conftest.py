"""Genuine data tests are explicitly skipped if the external kernel is absent."""
from pathlib import Path
import os
import pytest


@pytest.fixture(scope='session')
def de440_paths():
    bsp = Path(os.environ.get('DE440_BSP', 'data/de440.bsp'))
    gm = Path(os.environ.get('DE440_GM', 'data/gm_de440.tpc'))
    if not bsp.is_file() or not gm.is_file():
        pytest.skip('Genuine DE440 BSP and matching gm_de440.tpc absent; run scripts/download_data.py')
    return bsp, gm
