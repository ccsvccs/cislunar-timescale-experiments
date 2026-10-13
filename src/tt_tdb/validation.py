"""Explicitly synthetic analytic validation, NEVER a production ephemeris fallback."""
import numpy as np
from .constants import C_LIGHT


class KeplerReference:
    """Newtonian Sun-Earth ellipse with an exact analytic c^-2 residual integral.

    Derivation: docs/kepler_reference.md. Does not validate DE440 accuracy.
    """
    gm=1.32712440041279419e20
    a=149597870700.0
    e=0.0167086

    def __init__(self):
        self.n=np.sqrt(self.gm/self.a**3)
        self.period=2*np.pi/self.n
        self.mean=1.5*self.gm/(self.a*C_LIGHT**2)
        self.amplitude=2*self.gm*self.e/(self.a*C_LIGHT**2*self.n)

    def eccentric_anomaly(self,t):
        m=np.asarray(t)*self.n
        x=m.copy()
        for _ in range(12):
            step=(x-self.e*np.sin(x)-m)/(1-self.e*np.cos(x))
            x=x-step
            if np.max(np.abs(step))<1e-14:
                break
        return x

    def rate(self,t):
        e=self.eccentric_anomaly(t)
        radius=self.a*(1-self.e*np.cos(e))
        return (self.gm*(2/radius-1/(2*self.a))/C_LIGHT**2-self.mean)[...,None]

    def integral(self,a,b):
        return self.amplitude*(np.sin(self.eccentric_anomaly(b))-np.sin(self.eccentric_anomaly(a)))
