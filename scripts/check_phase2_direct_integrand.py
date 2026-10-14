"""Validate full-anchor direct JPLModel integration against term reconstruction."""
import json
from pathlib import Path
import time
import _common
import numpy as np
from tt_tdb.constants import T0_TDB, TDB0, A_SCALE
from tt_tdb.extended_ephemeris import ExtendedDE440Ephemeris
from tt_tdb.integration import SegmentedIntegrator
from tt_tdb.jpl_model import JPLModel
from tt_tdb.reference import difference_metrics


def main():
    started = time.perf_counter()
    table = np.genfromtxt('results/tables/phase2_time_series.csv', delimiter=',', names=True)
    times = table['tdb_j2000_seconds']
    cache = Path('results/cache/phase2_direct_gauss24.npz')
    with ExtendedDE440Ephemeris() as ephemeris:
        model = JPLModel(ephemeris)
        signature = model.fingerprint()
        parameters = dict(fingerprint=signature, anchor=T0_TDB.seconds_since_j2000(),
                          step_seconds=43200., method='gauss_legendre', order=24, channels=2)
        integrator = None
        if cache.exists():
            try:
                integrator = SegmentedIntegrator.load(cache, model, **parameters)
            except ValueError:
                pass
        if integrator is None:
            integrator = SegmentedIntegrator(model, **parameters)
        values = integrator.integral(times)
        integrator.save(cache)
        direct = TDB0 + A_SCALE * values.sum(axis=-1)
        direct_rate = A_SCALE * model(times).sum(axis=-1)
    result = {
        'command': '.venv/bin/python scripts/check_phase2_direct_integrand.py',
        'kind': 'direct_two_channel_final_JPLModel_from_original_1977_anchor',
        'model_signature': json.loads(signature), 'anchor_seconds': integrator.anchor,
        'sample_count': len(times), 'method': 'gauss_legendre', 'order': 24,
        'step_days': .5, 'channels': 2,
        'direct_minus_component_gauss': difference_metrics(times, direct-table['gauss_delta_s']),
        'direct_minus_component_romberg': difference_metrics(times, direct-table['extended_delta_s']),
        'direct_rate_minus_component_rate': difference_metrics(times, direct_rate-table['extended_rate'], value_unit='dimensionless'),
        'direct_raw_jpl_residual': difference_metrics(times, direct-table['jpl_delta_s']),
        'stats': integrator.stats, 'elapsed_wall_seconds': time.perf_counter()-started,
        'limitations': 'This measures numerical/source-path consistency, not completeness of unavailable small-body physics.'}
    _common.write_json('results/tables/phase2_direct_integrand_check.json', result)
    print(json.dumps({key: result[key] for key in ('direct_minus_component_gauss', 'direct_minus_component_romberg')}, indent=2))


if __name__ == '__main__':
    main()
