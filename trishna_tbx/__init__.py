"""trishna_tbx: what can TRISHNA VSWIR (5 VNIR + cirrus + 1610 nm, 60 m) monitor over water,
compared with Sentinel-2, before launch.

Building blocks reused: GRSdriver (L1C reading, SRF), grs (atmospheric correction),
GRSl2bgen (L2B algorithms, OWT references), radcalnet_oc (TOA simulation, gas LUT),
invRrs (semi-analytical inversion), hgrs (coupled water-vapour retrieval).
"""
__version__ = '0.1.1'

from .sensor import TRISHNA, TRISHNA_BANDS, S2_BANDS, make_s2, make_s2_trishna_like, make_trishna  # noqa
