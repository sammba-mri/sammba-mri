"""
Helper functions to download NeuroImaging mouse datasets and atlases
"""

from .func import (fetch_zurich_test_retest, fetch_zurich_anesthesiant)
from .atlas import (fetch_atlas_dorr_2008, fetch_atlas_waxholm_rat_2014,
                    fetch_masks_dorr_2008, fetch_atlas_lemur_mircen_2017)

__all__ = ['fetch_zurich_test_retest', 'fetch_zurich_anesthesiant',
           'fetch_atlas_dorr_2008', 'fetch_atlas_waxholm_rat_2014',
           'fetch_masks_dorr_2008', 'fetch_atlas_lemur_mircen_2017']
