"""External monopole gravitational potential, Eq. (21), SI inputs only."""
import numpy as np


def external_potential(earth_position, body_positions, gm, *, breakdown=False):
    """Return sum(GM/r) in m²/s²; never include Earth or EMB self-mass.

    Positions have shape (...,3) on the same TDB-compatible ICRF system.
    The finite body set is explicit; system barycenters use system GMs.
    """
    if not body_positions or body_positions.keys() != gm.keys():
        raise ValueError('Positions and GM must name the same nonempty body set')
    if {'earth','emb','399','3',399,3}.intersection(body_positions):
        raise ValueError('Earth/EMB cannot enter geocentric external potential')
    earth = np.asarray(earth_position, dtype=float)
    if earth.shape[-1:] != (3,) or not np.all(np.isfinite(earth)):
        raise ValueError('Earth position must be finite Cartesian SI components')
    contributions = {}
    for body, position in body_positions.items():
        mass = float(gm[body])
        if not np.isfinite(mass) or mass <= 0:
            raise ValueError('GM must be finite and positive in SI')
        position = np.asarray(position, dtype=float)
        if position.shape != earth.shape or not np.all(np.isfinite(position)):
            raise ValueError('Body and Earth position shapes must agree')
        radius = np.linalg.norm(position - earth, axis=-1)
        if np.any(radius <= 0):
            raise ValueError('External body coincides with Earth')
        contributions[body] = mass/radius
    total = np.sum(np.stack(list(contributions.values())), axis=0)
    return (total, contributions) if breakdown else total
