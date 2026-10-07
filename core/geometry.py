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


def magnitude_breakdown(design,num_g,den_g,num_h,den_h):
    """Decompõe a condição de módulo no formato usado na resolução manual.

    A_i são as distâncias do polo desejado aos polos da malha aberta
    compensada sem Kc; B_i são as distâncias aos zeros. K_T é a razão
    produto(A_i)/produto(B_i). O ganho do controlador é obtido de
    K_T = |K_G K_H| K_C.
    """
    sd=complex(design.desired_pole)
    contributions=design.plant_contributions+design.controller_contributions

    pole_terms=[]
    zero_terms=[]
    for c in contributions:
        distance=abs(sd-complex(c.singularity))
        item={
            "source_label":c.label,
            "singularity":complex(c.singularity),
            "dx":float(c.dx),
            "dy":float(c.dy),
            "distance":float(distance),
        }
        if c.kind=="pole":
            item["term"]=f"A{len(pole_terms)+1}"
            pole_terms.append(item)
        else:
            item["term"]=f"B{len(zero_terms)+1}"
            zero_terms.append(item)

    prod_a=float(np.prod([x["distance"] for x in pole_terms])) if pole_terms else 1.0
    prod_b=float(np.prod([x["distance"] for x in zero_terms])) if zero_terms else 1.0
    kt=prod_a/prod_b

    num_g=np.trim_zeros(np.asarray(num_g,dtype=float),"f")
    den_g=np.trim_zeros(np.asarray(den_g,dtype=float),"f")
    num_h=np.trim_zeros(np.asarray(num_h,dtype=float),"f")
    den_h=np.trim_zeros(np.asarray(den_h,dtype=float),"f")
    kg=float(num_g[0]/den_g[0])
    kh=float(num_h[0]/den_h[0])
    kgh=kg*kh
    if abs(kgh)<1e-15:
        raise ValueError("Ganho constante KG*KH nulo na decomposição de módulo.")

    kc=kt/abs(kgh)
    return {
        "poles":pole_terms,
        "zeros":zero_terms,
        "prod_a":prod_a,
        "prod_b":prod_b,
        "kt":float(kt),
        "kg":kg,
        "kh":kh,
        "kgh":kgh,
        "kc":float(kc),
    }


def angle_breakdown(design,num_g,den_g,num_h,den_h):
    """Dados da condição de ângulo no formato da resolução manual."""
    plant_zeros=[c for c in design.plant_contributions if c.kind=="zero"]
    plant_poles=[c for c in design.plant_contributions if c.kind=="pole"]
    controller_poles=[c for c in design.controller_contributions if c.kind=="pole"]
    controller_zeros=[c for c in design.controller_contributions if c.kind=="zero"]

    kg=float(np.trim_zeros(np.asarray(num_g,dtype=float),"f")[0]/np.trim_zeros(np.asarray(den_g,dtype=float),"f")[0])
    kh=float(np.trim_zeros(np.asarray(num_h,dtype=float),"f")[0]/np.trim_zeros(np.asarray(den_h,dtype=float),"f")[0])
    constant_phase=180.0 if kg*kh<0 else 0.0

    raw_base=(
        constant_phase
        +sum(c.angle_deg for c in plant_zeros)
        -sum(c.angle_deg for c in plant_poles)
        -sum(c.angle_deg for c in controller_poles)
    )
    controller_phase=sum(c.angle_deg for c in controller_zeros)
    final_raw=raw_base+controller_phase
    target=180.0+360.0*round((final_raw-180.0)/360.0)

    plant_ids={
        "zero":{c.label:i for i,c in enumerate(plant_zeros,1)},
        "pole":{c.label:i for i,c in enumerate(plant_poles,1)},
    }

    return {
        "plant_zeros":plant_zeros,
        "plant_poles":plant_poles,
        "controller_poles":controller_poles,
        "controller_zeros":controller_zeros,
        "constant_phase":constant_phase,
        "raw_base":raw_base,
        "controller_phase":controller_phase,
        "final_raw":final_raw,
        "target":target,
        "kg":kg,
        "kh":kh,
        "plant_ids":plant_ids,
    }
