"""Geometria trigonométrica do projeto pelo LGR."""
import cmath
import math

import numpy as np
import sympy as sp

from .models import TrigContribution


def phase_deg(value):
    return math.degrees(cmath.phase(complex(value)))


def phase_deficiency(base_phase_deg):
    return (180.0-base_phase_deg)%360.0


def split_zero_phase(total_phase,nzeros):
    angle=total_phase/nzeros
    if not 0.0<angle<180.0:
        raise ValueError("Fase incompatível com zeros reais nesta configuração.")
    return angle


def real_zero_parameter(sd,phi_deg):
    sigma=-complex(sd).real
    wd=abs(complex(sd).imag)
    t=math.tan(math.radians(phi_deg))
    if abs(t)<1e-12:
        raise ValueError("Ângulo degenerado.")
    return sigma+wd/t


def _roots_with_multiplicity(expr,s):
    poly=sp.Poly(expr,s)
    if poly.degree()<=0:
        return []
    try:
        roots=sp.roots(poly.as_expr(),s)
        if sum(int(m) for m in roots.values())==poly.degree():
            out=[]
            for root,m in roots.items():
                out.extend([complex(sp.N(root,16))]*int(m))
            return out
    except Exception:
        pass
    coeffs=np.asarray([complex(sp.N(c,16)) for c in poly.all_coeffs()],dtype=complex)
    return [complex(x) for x in np.roots(coeffs)]


def singularity_contributions(num,den,sd):
    """Contribuições angulares da função de malha aberta já reduzida.

    Cancelamentos exatos entre G(s) e H(s) são removidos antes de listar polos
    e zeros. Assim a tabela visual coincide com a função usada no LGR.
    """
    s=sp.symbols("s")
    num=np.trim_zeros(np.asarray(num,dtype=float),"f")
    den=np.trim_zeros(np.asarray(den,dtype=float),"f")

    num_coeffs=[sp.nsimplify(float(c),rational=True,tolerance=1e-12) for c in num]
    den_coeffs=[sp.nsimplify(float(c),rational=True,tolerance=1e-12) for c in den]
    num_expr=sp.Poly.from_list(num_coeffs,gens=s).as_expr()
    den_expr=sp.Poly.from_list(den_coeffs,gens=s).as_expr()
    reduced=sp.cancel(num_expr/den_expr)
    rnum,rden=sp.fraction(reduced)

    zeros=_roots_with_multiplicity(rnum,s)
    poles=_roots_with_multiplicity(rden,s)

    out=[]
    for i,z in enumerate(zeros,1):
        d=sd-z
        out.append(TrigContribution(f"z{i}",complex(z),float(d.real),float(d.imag),phase_deg(d),"zero"))
    for i,p in enumerate(poles,1):
        d=sd-p
        out.append(TrigContribution(f"p{i}",complex(p),float(d.real),float(d.imag),phase_deg(d),"pole"))
    return out
