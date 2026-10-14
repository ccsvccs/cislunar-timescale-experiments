"""Independent analytic checks of each verified Park term and real SPK input."""
from dataclasses import replace
import json
from types import SimpleNamespace
import numpy as np
import pytest
import spiceypy as spice
from tt_tdb.constants import C_LIGHT, C_EXACT
from tt_tdb.extended_ephemeris import ExtendedDE440Ephemeris, ExtendedEphemerisState, type2_state_derivatives
from tt_tdb.jpl_model import JPLModel, JPLTermSwitches, solar_j2_potential, kuiper_ring_potential, vector_potential, pn_scalar_components, read_de440_constants
from tt_tdb.pn_rate import rate_c2,rate_c4


def simple_state():
    # Deliberately analytic unit-test input; never an ephemeris fallback.
    return ExtendedEphemerisState(np.array([0.,0.,0.]),np.array([0.,3.,0.]),
          {'sun':np.array([2.,0.,0.]),'moon':np.array([0.,4.,0.])},
          {'sun':10.,'moon':20.},
          {'sun':np.array([3.,4.,0.]),'moon':np.array([0.,2.,0.])},
          {'sun':np.array([-2.,0.,0.]),'moon':np.array([0.,-3.,0.])},30.)


def test_vector_potential_exact_velocity_sum():
    s=simple_state()
    np.testing.assert_array_equal(vector_potential(s.earth_position,s.body_positions,s.body_velocities,s.gm),[15.,30.,0.])
    with pytest.raises(ValueError,match='velocities'):
        vector_potential(s.earth_position,s.body_positions,{},s.gm)


def test_scalar_components_outer_exclusion_and_inner_earth_inclusion():
    s=simple_state()
    d=pn_scalar_components(s)
    # U_sun=U_moon=5. Distances sun-moon=sqrt(20). Inner potentials include Earth!
    expected={'velocity':5*(-2*25)+5*(-2*4),
              'body_potential':5*(30/2+20/np.sqrt(20))+5*(30/4+10/np.sqrt(20)),
              'radial_velocity':5*0.5*3**2+5*0.5*2**2,
              'acceleration':5*0.5*4+5*0.5*12}
    for name,value in expected.items(): assert d[name] == pytest.approx(value,rel=1e-15)
    # Removing only Earth GM changes precisely the sum of GM_i GM_E/r_iE².
    less=pn_scalar_components(replace(s,earth_gm=15))
    assert d['body_potential']-less['body_potential'] == pytest.approx(5*15/2+5*15/4)
    for name in ('velocity','radial_velocity','acceleration'):
        assert d[name] == less[name]
    with pytest.raises(ValueError,match='Earth GM'):
        pn_scalar_components(replace(s,earth_gm=0))


def test_solar_j2_equator_pole_and_radius_power():
    equator=solar_j2_potential(np.array([2.,0.,0.]),np.zeros(3),10.,0.2,3.,[0,0,1])
    pole=solar_j2_potential(np.array([0.,0.,2.]),np.zeros(3),10.,0.2,3.,[0,0,1])
    assert equator == 10*.2*9/(2*8)
    assert pole == -2*equator
    far=solar_j2_potential(np.array([4.,0.,0.]),np.zeros(3),10.,0.2,3.,[0,0,1])
    assert far == equator/8
    with pytest.raises(ValueError,match='unit vector'):
        solar_j2_potential(np.array([2.,0.,0.]),np.zeros(3),10.,0.2,3.,[0,0,2])


def test_header_parameters_are_original_de440():
    h=read_de440_constants('docs/sources/de440_header.440')
    assert h['ASUN'] == 696000
    assert h['J2SUN'] == pytest.approx(2.1961391516529825e-7,rel=1e-16)
    assert h['BETA'] == h['GAMMA'] == 1


def model_with_test_state():
    model=JPLModel(SimpleNamespace(states=lambda t:simple_state(),missing_potential_sources=('analytic test fixture',)))
    model.solar_radius=3.  # finite analytic fixture, not the actual solar radius
    return model


def test_every_switch_disables_only_its_term_and_no_constant_changes():
    m=model_with_test_state(); original=m.contribution_rates(0.)
    for name in original:
        m.switches=replace(JPLTermSwitches(),**{name:False})
        changed=m.contribution_rates(0.)
        assert changed[name] == 0
        for other in original:
            if other != name: assert changed[other] == original[other]
        # Rate difference is the disabled physical term, no C fit or zero point.
        m.switches=JPLTermSwitches(); full=m(0.)
        m.switches=replace(JPLTermSwitches(),**{name:False}); less=m(0.)
        channel=0 if name.endswith('_c2') else 1
        assert full[channel]-less[channel] == pytest.approx(original[name],rel=2e-13,abs=3e-24 if channel == 0 else 1e-44)


def test_park_rate_sign_and_scalar_terms():
    m=model_with_test_state();d=m.diagnostics(0.);s=simple_state()
    expected=rate_c2(9,10)-C_EXACT+rate_c4(9,10)
    expected+=solar_j2_potential(s.earth_position,s.body_positions['sun'],10,m.solar_j2,m.solar_radius,m.solar_axis)/C_LIGHT**2
    expected+=kuiper_ring_potential(s.earth_position,s.body_positions['sun'],m.ring_gm,radius=m.ring_radius)/C_LIGHT**2
    expected-=4*90/C_LIGHT**4
    expected-=sum(pn_scalar_components(s).values())/C_LIGHT**4
    assert m(0.).sum() == pytest.approx(expected,rel=1e-15)
    assert tuple(m.correction_rates(0)) == m.correction_names
    assert d['observer_position_term_seconds'] == 0
    for term in ('pn_velocity_c4','pn_body_potential_c4','pn_radial_velocity_c4','pn_acceleration_c4'):
        assert d['terms'][term] == -d['pn_scalar_components'][term[3:-3]]/C_LIGHT**4


class PolynomialSegment:
    data_type=2
    def load_array(self):
        # x-position=3 T0+2 T1+4 T2+5 T3, y,z=0; one two-day record.
        c=np.zeros((3,1,4));c[0,0]=[3,2,4,5]
        return 2451545.,2.,c


def test_type2_second_derivative_against_exact_polynomial():
    jd1=np.array([2451545.,2451546.,2451547.]); x=np.array([-1.,0.,1.])
    p,v,a=type2_state_derivatives(PolynomialSegment(),jd1)
    np.testing.assert_allclose(p[:,0],3+2*x+4*(2*x*x-1)+5*(4*x**3-3*x),rtol=0,atol=0)
    np.testing.assert_allclose(v[:,0],(2+16*x+5*(12*x*x-3))/86400,rtol=2e-16,atol=0)
    np.testing.assert_allclose(a[:,0],(16+120*x)/86400**2,rtol=2e-16,atol=0)
    with pytest.raises(ValueError,match='coverage'):
        type2_state_derivatives(PolynomialSegment(),2451548.)


@pytest.mark.de440
def test_extended_velocities_against_independent_spice_and_acceleration(de440_paths):
    spice.furnsh(str(de440_paths[0]))
    try:
        with ExtendedDE440Ephemeris(*de440_paths) as e:
            t=np.array([-725803167.8160655,0.,8e8]); s=e.states(t)
            for n,target in [('sun',10),('moon',301),('jupiter',5)]:
                for k,epoch in enumerate(t):
                    ref,_=spice.spkezr(str(target),float(epoch),'J2000','NONE','0')
                    np.testing.assert_allclose(s.body_positions[n][k],ref[:3]*1000,rtol=0,atol=.01)
                    np.testing.assert_allclose(s.body_velocities[n][k],ref[3:]*1000,rtol=0,atol=2e-10)
                    # Independent velocity central derivative; h/h2 tests
                    # reject a unit, sign or Chebyshev factor error in a.
                    h=10.
                    before,_=spice.spkezr(str(target),float(epoch-h),'J2000','NONE','0')
                    after,_=spice.spkezr(str(target),float(epoch+h),'J2000','NONE','0')
                    numeric=(after[3:]-before[3:])*1000/(2*h)
                    np.testing.assert_allclose(s.body_accelerations[n][k],numeric,rtol=2e-7,atol=1e-11)
            model=JPLModel(e)
            assert model(t).shape == (3,2)
            fp=json.loads(model.fingerprint());model.switches=replace(model.switches,vector_c4=False)
            assert json.loads(model.fingerprint())['switches'] != fp['switches']
            assert fp['missing']
            assert not any('jpl_reference' in n for n in fp['sources'])
    finally:
        spice.unload(str(de440_paths[0]))


def test_published_ring_centre_mass_radius_and_phase_invariance():
    gm=8.916783186240352e12;au=149597870700.;r=44*au
    centre=kuiper_ring_potential(np.zeros(3),np.zeros(3),gm,radius=r)
    assert centre == pytest.approx(gm/r,rel=2e-16)
    earth=np.array([au,0.,0.])
    values=[kuiper_ring_potential(earth,np.zeros(3),gm,radius=r,phase=p) for p in (0,.017,.3,1.9)]
    np.testing.assert_allclose(values,values[0],rtol=5e-16,atol=0)
    # For r_E/R=.024 the conservative phase bound is astronomically below ps.
    q=1.02/44
    rate_bound=4*(gm/r)*q**36/(1-q)/C_LIGHT**2
    assert rate_bound*1100*365.25*86400 < 1e-60
    with pytest.raises(ValueError,match='interior'):
        kuiper_ring_potential(np.array([r,0,0]),np.zeros(3),gm,radius=r)


@pytest.mark.de440
def test_anchored_converter_cache_and_term_ablation(de440_paths,tmp_path):
    from tt_tdb.integration import SegmentedIntegrator
    from tt_tdb.transforms import TTTDBConverter
    from tt_tdb.constants import T0_TDB,TDB0
    with ExtendedDE440Ephemeris(*de440_paths) as e:
        model=JPLModel(e);anchor=T0_TDB.seconds_since_j2000()
        integrator=SegmentedIntegrator(model,anchor=anchor,method='gauss_legendre',order=12,
                    fingerprint=model.fingerprint())
        converter=TTTDBConverter(integrator)
        assert converter.tdb_minus_tt(T0_TDB) == TDB0
        epoch=T0_TDB.add_seconds(2*86400)
        original=converter.tdb_minus_tt(epoch)
        tt=converter.tdb_to_tt(epoch)
        assert converter.tt_to_tdb(tt).difference_seconds(epoch) == pytest.approx(0.,abs=1e-14)
        path=tmp_path/'jpl_candidate.npz';integrator.save(path)
        cached=SegmentedIntegrator.load(path,model,fingerprint=model.fingerprint())
        assert TTTDBConverter(cached).tdb_minus_tt(epoch) == original
        model.switches=replace(model.switches,kuiper_ring_c2=False)
        with pytest.raises(ValueError,match='fingerprint'):
            SegmentedIntegrator.load(path,model,fingerprint=model.fingerprint())
        no_ring=TTTDBConverter(SegmentedIntegrator(model,anchor=anchor,
                  method='gauss_legendre',order=12,fingerprint=model.fingerprint()))
        difference=original-no_ring.tdb_minus_tt(epoch)
        assert 2.5e-12 < difference < 2.7e-12
        # Independent methods on an uncalibrated fixed-anchor short interval.
        romberg=TTTDBConverter(SegmentedIntegrator(model,anchor=anchor,method='romberg',
                 order=7,fingerprint=model.fingerprint()))
        assert abs(romberg.tdb_minus_tt(epoch)-no_ring.tdb_minus_tt(epoch)) < 1e-16
