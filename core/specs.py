"""Conversão de especificações de desempenho em polos dominantes."""
from dataclasses import dataclass
import math
from typing import Optional

@dataclass
class DesignSpecs:
    zeta: Optional[float]=None
    omega_n: Optional[float]=None
    mp_max: Optional[float]=None
    ts_max: Optional[float]=None
    settling_band: float=0.02
    pole: Optional[complex]=None

    @classmethod
    def from_zeta_wn(cls,zeta,omega_n):
        return cls(zeta=float(zeta),omega_n=float(omega_n))

    @classmethod
    def from_mp_ts(cls,mp_max,ts_max,settling_band):
        band=float(settling_band)
        if band not in (0.02,0.05):
            raise ValueError("O critério de acomodação deve ser 2% ou 5%.")
        return cls(mp_max=float(mp_max),ts_max=float(ts_max),settling_band=band)

    @classmethod
    def from_pole(cls,real,imag):
        return cls(pole=complex(float(real),abs(float(imag))))

def zeta_from_overshoot(mp_percent):
    if not 0 < mp_percent < 100:
        raise ValueError("Mp deve estar entre 0 e 100%.")
    ln_mp=math.log(mp_percent/100.0)
    return -ln_mp/math.sqrt(math.pi**2+ln_mp**2)

def settling_sigma(ts,band):
    if ts<=0: raise ValueError("ts deve ser positivo.")
    if abs(band-0.05)<1e-12: return 3.0/ts
    if abs(band-0.02)<1e-12: return 4.0/ts
    raise ValueError("Use 2% ou 5%.")

def desired_pole(zeta,omega_n):
    if not 0<zeta<1: raise ValueError("Use 0 < ξ < 1.")
    if omega_n<=0: raise ValueError("ωn deve ser positivo.")
    sigma=zeta*omega_n
    wd=omega_n*math.sqrt(1-zeta**2)
    return complex(-sigma,wd)

def resolve_specs(specs):
    if specs.pole is not None:
        p=complex(specs.pole)
        if p.real>=0 or p.imag<=0: raise ValueError("Polo desejado inválido.")
        wn=abs(p); zeta=-p.real/wn
        return p,zeta,wn,-p.real
    if specs.zeta is not None and specs.omega_n is not None:
        p=desired_pole(specs.zeta,specs.omega_n)
        return p,specs.zeta,specs.omega_n,-p.real
    if specs.mp_max is not None and specs.ts_max is not None:
        zeta=zeta_from_overshoot(specs.mp_max)
        sigma=settling_sigma(specs.ts_max,specs.settling_band)
        wn=sigma/zeta
        return desired_pole(zeta,wn),zeta,wn,sigma
    raise ValueError("Especificações incompletas.")
