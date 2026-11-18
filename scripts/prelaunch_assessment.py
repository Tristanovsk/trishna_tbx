"""Reproduce the pre-launch TRISHNA vs S2 numbers quoted in README.md (no image needed).

python scripts/prelaunch_assessment.py [--inversion]
"""
import sys

import numpy as np
import pandas as pd

from trishna_tbx import TRISHNA, algorithms, make_s2, owt
from trishna_tbx.emulate import SBAF
from trishna_tbx.water_vapor import compare_trishna_s2, noise_error

pd.set_option('display.width', 200)
s2 = make_s2()
S2W = s2.subset(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8A'], 'S2-water')
S2T = s2.subset(['B2', 'B3', 'B4', 'B8A'], 'S2-as-TRISHNA')
TRW = TRISHNA.subset(['Blue', 'Green', 'Red', 'NIR'], 'TRISHNA-water')
TR3 = TRISHNA.subset(['Blue', 'Green', 'Red'], 'TRISHNA-3vis')

print('== 1. algorithm coverage (tol 15 nm)')
print(algorithms.coverage([s2, TRISHNA]))

print('\n== 2. OWT (Bi et al. 2024, SAM) agreement with hyperspectral class, %')
df, _ = owt.benchmark({'S2-water': S2W, 'S2-as-TRISHNA': S2T, 'TRISHNA-water': TRW, 'TRISHNA-3vis': TR3}, n=300)
print('noise-free\n', (100 * df).round(0).T)
k = 1 / (np.pi * 0.75)  # TOA reflectance noise -> Rrs noise
nz = {'S2-water': S2W.noise_reflectance(30, n_aggregated=9).assign_coords(band=S2W.band_names) * k,
      'S2-as-TRISHNA': S2T.noise_reflectance(30, n_aggregated=9).assign_coords(band=S2T.band_names) * k,
      'TRISHNA-water': TRW.noise_reflectance(30).assign_coords(band=TRW.band_names) * k}
df, _ = owt.benchmark({k_: v for k_, v in {'S2-water': S2W, 'S2-as-TRISHNA': S2T, 'TRISHNA-water': TRW}.items()},
                      n=300, noise=nz)
print('with SNR-spec noise (S2 aggregated 20->60 m)\n', (100 * df).round(0).T)
mean, sd = owt.load_reference('Bi2024')
print('\nconfusable pairs on TRISHNA (spectral angle < 6 deg):')
sep = owt.separability(mean, TRW).to_pandas()
print([(i, j, round(sep.loc[i, j], 1)) for i in sep.index for j in sep.columns if i < j and sep.loc[i, j] < 6])

print('\n== 3. S2 -> TRISHNA band adjustment on Bi2024 members')
m = owt.sample_members(mean, sd, n=100).stack(s=('true_owt', 'member'))
m = m.drop_vars(['s', 'true_owt', 'member']).transpose('s', 'wl')
print(SBAF(m).table())

print('\n== 4. water vapour: TCWV error (kg m-2) for 1% ratio error, M=3')
c = compare_trishna_s2(3.)
print(c.dW_for_1pct_ratio_error.sel(W=[5, 20, 40, 60]).to_pandas().round(2))
print('noise-only TCWV error (kg m-2), flat surface reflectance rho, M=3')
for rho in (0.3, 0.05, 0.01):
    for s in ('TRISHNA', 'S2'):
        e = noise_error(s, rho, m=3., sza=48.19)
        print(f'  rho={rho:<5} {s:8s}', ' '.join(f'W={w}:{float(e.sel(W=w)):.2f}' for w in (5, 20, 40, 60)))

if '--inversion' in sys.argv:
    print('\n== 5. invRrs MultispectralSimple: linearised relative std (5% + 2e-4 sr-1 noise)')
    import invRrs.inverse as inv
    from invRrs import utils as u
    from trishna_tbx.inversion import information_content, trishna_forward_model
    aphy, bbphy = u.data().load_iop_phyto('Bricaud1995')
    models = {'TRISHNA': trishna_forward_model('MultispectralSimple', aphy_star=aphy, bbphy_star=bbphy, sza=30.),
              'S2': inv.MultispectralSimple(aphy, bbphy, 30., sat='s2a', band_idx=[0, 1, 2, 3, 4, 5, 6, 8])}
    for state in ([0.5, 0.05, 0.2, 0], [2, 0.5, 1.0, 0], [20, 1.0, 10., 0]):
        for name, mdl in models.items():
            r = np.asarray(mdl.forward_model(np.asarray(state, float), output='Rrs'))
            ic = information_content(mdl, state, 0.05 * r + 2e-4)
            rel = ic['std'].values[:3] / np.abs(state[:3])
            print(f'  {name:8s} chl={state[0]:<4} acdom={state[1]:<4} spm={state[2]:<4} -> '
                  f'chl {rel[0]:.2f}  acdom {rel[1]:.2f}  spm {rel[2]:.2f}  corr(chl,cdom) {float(ic.corr[0, 1]):+.2f}')
