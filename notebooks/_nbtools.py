"""Helpers to write the VRE gallery notebooks with a common layout.

Layout of the Water Quality VRE (ocean-and-coastal-water-quality-vre.netlify.app):
Title / Requirements / Tutorial objectives / Environment / Resources / Imports / numbered
sections. Added here: "Data" (editable paths) and "Take-home messages" / "References".
"""
import nbformat as nbf

VERSION_BADGE = '![version](https://img.shields.io/badge/version-11.2026-blue)'

DATA_NOTE = (
    "> **Data paths.** On the CNES JupyterHub, point the variables below to the TREX datalake or to "
    "your project space. Outside CNES, download Sentinel-2 L1C with the *Access data - Copernicus "
    "API* tutorial. Each path can also be set through the environment variable named in the cell.")

COMMON_DATA = '''import os
from pathlib import Path

# ---- edit here (or set the environment variables) ------------------------------------
S2_L1C = os.environ.get('TRISHNA_TBX_S2_L1C',
    '/data/satellite/Sentinel-2/L1C/31TFJ/2025/05/26/'
    'S2B_MSIL1C_20250526T103629_N0511_R008_T31TFJ_20250526T123708.SAFE')
CAMS_FILE = os.environ.get('TRISHNA_TBX_CAMS', '/data/cams/world/cams_forecast_2025-05.nc')
# 30 km box (UTM 31N) over the Berre lagoon and the eastern Gulf of Fos
SUBSET = [(655000, 4800000, 685000, 4830000)]
L2A_DIR = os.environ.get('TRISHNA_TBX_L2A_DIR',
    '/home/harmel/Dropbox/Dropbox/work/projet/unesco/satellite/L2A')  # GRS L2A (beam profile)
L2A_TILE = 'T33PVQ'                      # Lake Chad (southern pool)
L2A_WINDOW = 'y=0:3600,x=0:960'          # 20 m pixel window with water
L2A_DATE = '2022-02-12'                  # date used by the single-image tutorials
# single-image example of the products tutorial: Loire estuary and plume, GRS v3.1.1 (Zarr pyramid,
# local S3 bucket grstbx-s3 of grstbx, read here through its folder)
L2A_EXAMPLE = os.environ.get('TRISHNA_TBX_L2A_EXAMPLE',
    '/data/satellite/s3/grstbx-s3/sentinel-2/l2a/30TWT/2017/04/12/'
    'S2A_MSIL2AGRS_20170412T110621_N0500_R137_T30TWT_20231111T110836_V3.1.1/'
    'S2A_MSIL2AGRS_20170412T110621_N0500_R137_T30TWT_20231111T110836_V3.1.1.zarr')
L2A_EXAMPLE_WINDOW = 'y=2700:5490,x=1500:4500'   # 20 m pixels: estuary, plume, Bourgneuf Bay
L2B_EXAMPLE = os.environ.get('TRISHNA_TBX_L2B_EXAMPLE',   # GRSl2bgen L2B of the same image
    '/data/satellite/s3/grstbx-s3/sentinel-2/l2b/30TWT/2017/04/12/'
    'S2A_MSIL2B_20170412T110621_N0500_R137_T30TWT_20231111T110836_V3.1.1/'
    'S2A_MSIL2B_20170412T110621_N0500_R137_T30TWT_20231111T110836_V3.1.1.zarr')
AERONET_HYP_DIR = os.environ.get('TRISHNA_TBX_AERONET_HYP', '/DATA/AERONET/OCv3_invRrs/data')
OUTPUT_DIR = Path(os.environ.get('TRISHNA_TBX_OUTPUT', './output')).resolve()
# ---------------------------------------------------------------------------------------
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
'''


def header(title, requirements, objectives, extra_env=(), resources='1 core, 4 GB', runtime='< 2 min'):
    env = ['GRSdriver', 'grs', 'GRSl2bgen', 'radcalnet_oc', 'invRrs', 'trishna_tbx'] + list(extra_env)
    md = [f'# {title}',
          '## Requirements\n' + '\n'.join(f'- {r}' for r in requirements),
          '## Tutorial objectives\n' + '\n'.join(f'- {o}' for o in objectives),
          '## Environment\n'
          'Use the `sandbox_env` kernel of the CNES JupyterHub. Otherwise install the packages below in '
          'your environment (the `trishna_tbx` toolbox: `pip install -e <path>/trishna_tbx`):\n\n'
          + ', '.join(f'`{e}`' for e in env) + ', `matplotlib`, `pandas`, `xarray`.',
          f'## Resources\nThe code can be executed with the CPU: {resources}. Run time: {runtime}.\n\n'
          + VERSION_BADGE]
    return [nbf.v4.new_markdown_cell(m) for m in md]


def md(text):
    return nbf.v4.new_markdown_cell(text.strip('\n'))


def code(text):
    return nbf.v4.new_code_cell(text.strip('\n'))


def data_cells():
    return [md('## Data\n' + DATA_NOTE), code(COMMON_DATA)]


IMPORTS = '''import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

import trishna_tbx
from trishna_tbx import TRISHNA, make_s2
print('trishna_tbx', trishna_tbx.__version__)
from pathlib import Path
THUMB_DIR = Path('../thumbnails').resolve() if Path('../thumbnails').exists() else Path('.').resolve()
# plot style: Times-like serif fonts with LaTeX-style math (STIX, no TeX installation needed),
# larger fonts, minor ticks on all axes
plt.rcParams.update({'figure.dpi': 90, 'font.size': 14, 'font.family': 'STIXGeneral',
                     'mathtext.fontset': 'stix', 'axes.titlesize': 15, 'axes.labelsize': 15,
                     'legend.fontsize': 12, 'xtick.labelsize': 13, 'ytick.labelsize': 13,
                     'xtick.minor.visible': True, 'ytick.minor.visible': True})
'''


def imports(extra=''):
    return [md('## Imports'), code(IMPORTS + extra)]


def write(path, cells, kernel='python3'):
    nb = nbf.v4.new_notebook()
    nb['cells'] = cells
    nb['metadata'] = {'kernelspec': {'name': kernel, 'display_name': 'Python 3', 'language': 'python'},
                      'language_info': {'name': 'python'}}
    nbf.write(nb, path)
