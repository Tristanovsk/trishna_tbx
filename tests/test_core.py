import numpy as np
import pytest

from trishna_tbx import TRISHNA, make_s2, algorithms


def test_trishna_srf_fwhm():
    for b in TRISHNA.bands:
        r = TRISHNA.srf.sel(band=b.name)
        half = r.wl.where(r >= 0.5, drop=True)
        assert abs(float(half.max() - half.min()) - b.fwhm) <= 2.
        assert abs(float((r * r.wl).sum() / r.sum()) - b.center) < 0.5


def test_noise_positive_and_aggregation():
    n1 = TRISHNA.noise_reflectance(30.)
    n9 = TRISHNA.noise_reflectance(30., n_aggregated=9)
    assert (n1 > 0).all()
    np.testing.assert_allclose(n1 / n9, 3.)


def test_coverage_flags_red_edge():
    cov = algorithms.coverage([make_s2(), TRISHNA])
    assert not cov.loc[('Chla_index', 'NDCI'), 'TRISHNA']
    assert cov.loc[('SPM', 'Nechad10'), 'TRISHNA']


def test_strict_run_refuses_missing_band():
    import xarray as xr
    rrs = xr.DataArray([0.01, 0.01, 0.005, 0.001], dims='wl', coords={'wl': [485, 555, 670, 860]})
    ndci = next(a for a in algorithms.REGISTRY if a.name == 'NDCI')
    with pytest.raises(ValueError):
        ndci.run(rrs, TRISHNA)


@pytest.mark.parametrize('sensor', ['TRISHNA', 'S2'])
def test_water_vapor_closure(sensor):
    pytest.importorskip('radcalnet_oc')
    from trishna_tbx.water_vapor import RatioRetrieval
    r = RatioRetrieval(sensor)
    ang = np.degrees(np.arccos(1 / 1.5))  # M = 3
    w = r.retrieve(np.array([1.]), r.forward(25., 3.).values * np.ones(1), ang, ang)
    np.testing.assert_allclose(w, 25., atol=0.05)


def test_owt_benchmark_runs():
    pytest.importorskip('GRSl2bgen')
    from trishna_tbx import owt
    tri = TRISHNA.subset(['Blue', 'Green', 'Red', 'NIR'], 'T')
    df, conf = owt.benchmark({'T': tri}, n=20)
    assert 0.5 < df.loc['all', 'T'] <= 1.


def test_metrics_and_l1c_restriction():
    import xarray as xr
    from trishna_tbx.validation import metrics
    from trishna_tbx.emulate import l1c_dataset_as_trishna
    x = np.linspace(0.001, 0.01, 20)
    m = metrics(x, 1.1 * x)
    assert m['N'] == 20 and abs(m['slope'] - 1.1) < 1e-9 and abs(m['MdAPD (%)'] - 10) < 1e-6
    wl = [443, 490, 560, 665, 705, 740, 783, 842, 865, 945, 1375, 1610, 2190]
    ds = xr.Dataset({'bands': (('wl', 'y', 'x'), np.random.rand(13, 6, 6)),
                     'sza': (('y', 'x'), np.full((6, 6), 30.)),
                     'surfwater': (('y', 'x'), np.ones((6, 6), 'uint8'))},
                    coords={'wl': wl, 'x': np.arange(6) * 20., 'y': np.arange(6) * 20.})
    out = l1c_dataset_as_trishna(ds)
    assert list(out.wl.values) == [490, 560, 665, 865, 945, 1375, 1610, 2190]
    assert list(out.attrs['wl_to_process']) == [490, 560, 665, 865, 1610]
    assert out.sizes['x'] == 2 and int(out.surfwater.min()) == 1
