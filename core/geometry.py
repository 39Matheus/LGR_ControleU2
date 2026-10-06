"""Geometria trigonométrica do projeto pelo LGR."""
import cmath, math
import numpy as np
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
    sigma=-complex(sd).real; wd=abs(complex(sd).imag)
    t=math.tan(math.radians(phi_deg))
    if abs(t)<1e-12: raise ValueError("Ângulo degenerado.")
    return sigma+wd/t

def singularity_contributions(num,den,sd):
    num=np.trim_zeros(np.asarray(num,dtype=float),"f"); den=np.trim_zeros(np.asarray(den,dtype=float),"f")
    zeros=np.roots(num) if len(num)>1 else []
    poles=np.roots(den) if len(den)>1 else []
    out=[]
    for i,z in enumerate(zeros,1):
        d=sd-z; out.append(TrigContribution(f"z{i}",complex(z),float(d.real),float(d.imag),phase_deg(d),"zero"))
    for i,p in enumerate(poles,1):
        d=sd-p; out.append(TrigContribution(f"p{i}",complex(p),float(d.real),float(d.imag),phase_deg(d),"pole"))
    return out
