import numpy as np
import pytest
from tt_tdb.potential import external_potential


def test_potential_units_contributions_and_translation_invariance():
    earth=np.array([0.,0.,0.])
    positions={'sun':np.array([1.5e11,0.,0.]),'moon':np.array([0.,4e8,0.])}
    gm={'sun':1.327e20,'moon':4.9e12}
    total,bodies=external_potential(earth,positions,gm,breakdown=True)
    assert bodies['sun']==gm['sun']/1.5e11
    assert bodies['moon']==12250.
    assert total==sum(bodies.values())
    shift=np.array([1e9,-2e9,3e9])
    assert external_potential(earth+shift,{k:v+shift for k,v in positions.items()},gm)==total
    with pytest.raises(ValueError):
        external_potential(earth,{'earth':positions['sun']},{'earth':1e15})
    with pytest.raises(ValueError):
        external_potential(earth,{'sun':earth},{'sun':1e20})
