"""Aproximação discreta de funções de transferência contínuas.

Implementa as substituições do capítulo 5 da apostila:
Forward/Euler, Backward e Tustin. A invariância ao degrau usa equivalência ZOH.
"""
from typing import Iterable, Tuple

import numpy as np
import sympy as sp

try:
    import control as ctrl
except ImportError:
    ctrl=None

from .models import DiscretizationResult


_METHOD_LABELS={
    "forward":"Euler / Forward",
    "euler":"Euler / Forward",
    "backward":"Backward",
    "tustin":"Tustin / Bilinear",
    "zoh":"Invariância ao degrau (SOZ)",
    "step":"Invariância ao degrau (SOZ)",
}


def normalize_method(method: str) -> str:
    key=str(method).strip().lower()
    aliases={
        "euler":"forward",
        "forward/euler":"forward",
        "euler / forward":"forward",
        "forward":"forward",
        "backward":"backward",
        "tustin":"tustin",
        "bilinear":"tustin",
        "tustin / bilinear":"tustin",
        "invariância ao degrau":"zoh",
        "invariancia ao degrau":"zoh",
        "zoh":"zoh",
        "step":"zoh",
    }
    if key not in aliases:
        raise ValueError("Método deve ser Forward/Euler, Backward, Tustin ou Invariância ao Degrau.")
    return aliases[key]


def _rational_coeffs(coeffs: Iterable[float]):
    vals=np.trim_zeros(np.asarray(list(coeffs),dtype=float),"f")
    if vals.size==0:
        raise ValueError("Polinômio nulo não é válido.")
    return [sp.nsimplify(float(v),rational=True,tolerance=1e-12) for v in vals]


def transfer_expr(num,den,var):
    n=sp.Poly.from_list(_rational_coeffs(num),gens=var).as_expr()
    d=sp.Poly.from_list(_rational_coeffs(den),gens=var).as_expr()
    return sp.cancel(n/d)


def _numeric_coeffs(expr,var) -> Tuple[list,list]:
    num,den=sp.fraction(sp.cancel(expr))
    pnum=sp.Poly(sp.expand(num),var)
    pden=sp.Poly(sp.expand(den),var)
    n=np.asarray([float(sp.N(v,16)) for v in pnum.all_coeffs()],dtype=float)
    d=np.asarray([float(sp.N(v,16)) for v in pden.all_coeffs()],dtype=float)
    # remove resíduos minúsculos e normaliza o denominador
    n[np.abs(n)<1e-12]=0.0
    d[np.abs(d)<1e-12]=0.0
    if abs(d[0])<1e-15:
        raise ValueError("Denominador discreto degenerado.")
    n=n/d[0]
    d=d/d[0]
    return n.tolist(),d.tolist()


def _format_coeff(value):
    value=float(value)
    if abs(value)<1e-12:
        return "0"
    if abs(value-round(value))<1e-10:
        return str(int(round(value)))
    return f"{value:.8g}"


def _poly_latex(coeffs,var):
    """Polinômio legível sem reconstruir racionais gigantes a partir de floats."""
    vals=np.asarray(coeffs,dtype=float)
    degree=len(vals)-1
    terms=[]
    for i,value in enumerate(vals):
        if abs(value)<1e-12:
            continue
        power=degree-i
        sign="-" if value<0 else "+"
        mag=abs(float(value))
        coeff="" if power>0 and abs(mag-1.0)<1e-12 else _format_coeff(mag)
        if power==0:
            body=coeff or "1"
        elif power==1:
            body=(coeff + r"\," if coeff else "") + sp.latex(var)
        else:
            body=(coeff + r"\," if coeff else "") + sp.latex(var) + rf"^{{{power}}}"
        terms.append((sign,body))
    if not terms:
        return "0"
    out=""
    for i,(sign,body) in enumerate(terms):
        if i==0:
            out+=("-" if sign=="-" else "")+body
        else:
            out+=f" {sign} {body}"
    return out


def transfer_latex(num,den,var):
    return rf"\frac{{{_poly_latex(num,var)}}}{{{_poly_latex(den,var)}}}"


def difference_equation_latex(num,den):
    """Forma causal em potências de z^-1, assumindo denominador normalizado."""
    den=np.asarray(den,dtype=float)
    num=np.asarray(num,dtype=float)
    n=max(len(den),len(num))-1
    den=np.pad(den,(n+1-len(den),0))
    num=np.pad(num,(n+1-len(num),0))
    if abs(den[0])<1e-15:
        # Os coeficientes chegam em ordem decrescente. Se o numerador/denominador
        # tiver grau menor que n, o padding deve ser à esquerda, mas o denominador
        # principal precisa continuar em 1. Recalcula sem elevar artificialmente n.
        n=len(den)-1
    # Melhor: use o grau real do denominador como referência.
    den=np.asarray(den,dtype=float)
    nz=np.flatnonzero(np.abs(den)>1e-12)
    if len(nz)==0:
        return ""
    den=den[nz[0]:]
    order=len(den)-1
    num=np.asarray(num,dtype=float)
    if len(num)<order+1:
        num=np.pad(num,(order+1-len(num),0))
    elif len(num)>order+1:
        # função imprópria: ainda apresenta a relação algébrica, mas sem tentar
        # impor forma causal de ordem menor.
        order=max(order,len(num)-1)
        den=np.pad(den,(order+1-len(den),0))
    if abs(den[0]-1)>1e-10:
        num=num/den[0]
        den=den/den[0]

    lhs=[r"y[k]"]
    rhs=[]
    for i,a in enumerate(den[1:],1):
        coef=-a
        if abs(coef)>1e-12:
            rhs.append((coef,rf"y[k-{i}]"))
    # num em ordem decrescente de z. Dividir por z^order gera b_i z^-i.
    if len(num)<order+1:
        num=np.pad(num,(order+1-len(num),0))
    for i,b in enumerate(num,0):
        if abs(b)>1e-12:
            rhs.append((b,rf"u[k-{i}]"))

    def fmt_term(coef,symbol,first=False):
        sign="-" if coef<0 else ("+" if not first else "")
        val=abs(coef)
        if abs(val-1)<1e-12:
            body=symbol
        else:
            body=rf"{val:.8g}\,{symbol}"
        return f"{sign}{body}"

    if not rhs:
        return r"y[k]=0"
    rendered=[]
    for i,(coef,sym) in enumerate(rhs):
        rendered.append(fmt_term(coef,sym,first=(i==0)))
    return r"y[k]="+"".join(rendered)


def _symbolic_discretize(num,den,T,method):
    s,z=sp.symbols("s z")
    Tq=sp.nsimplify(float(T),rational=True,tolerance=1e-12)
    g=transfer_expr(num,den,s)

    if method=="forward":
        sub=(z-1)/Tq
        sub_latex=rf"s\approx\frac{{z-1}}{{T}}=\frac{{z-1}}{{{sp.latex(Tq)}}}"
    elif method=="backward":
        sub=(z-1)/(Tq*z)
        sub_latex=rf"s\approx\frac{{z-1}}{{Tz}}=\frac{{z-1}}{{{sp.latex(Tq)}z}}"
    elif method=="tustin":
        sub=2*(z-1)/(Tq*(z+1))
        sub_latex=rf"s\approx\frac{{2}}{{T}}\frac{{z-1}}{{z+1}}=\frac{{2}}{{{sp.latex(Tq)}}}\frac{{z-1}}{{z+1}}"
    else:
        raise ValueError("Método simbólico inválido.")

    substituted=sp.cancel(sp.together(g.subs(s,sub)))
    num_d,den_d=_numeric_coeffs(substituted,z)
    return g,sub,sub_latex,substituted,num_d,den_d


def _zoh_discretize(num,den,T):
    if ctrl is None:
        raise RuntimeError("A invariância ao degrau requer o pacote python-control.")
    sys=ctrl.TransferFunction(num,den)
    ds=ctrl.sample_system(sys,float(T),method="zoh")
    num_d=np.asarray(ds.num[0][0],dtype=float)
    den_d=np.asarray(ds.den[0][0],dtype=float)
    if abs(den_d[0])>1e-15:
        num_d=num_d/den_d[0]
        den_d=den_d/den_d[0]
    num_d[np.abs(num_d)<1e-12]=0.0
    den_d[np.abs(den_d)<1e-12]=0.0
    return num_d.tolist(),den_d.tolist()


def discretize_transfer(num,den,sample_time,method,target_label="G(s)") -> DiscretizationResult:
    T=float(sample_time)
    if T<=0:
        raise ValueError("O período de amostragem T deve ser positivo.")
    key=normalize_method(method)
    s,z=sp.symbols("s z")
    continuous=transfer_expr(num,den,s)

    notes=[]
    if key=="zoh":
        num_d,den_d=_zoh_discretize(num,den,T)
        sub_latex=r"H(z)=\frac{z-1}{z}\,\mathcal{Z}\left\{\frac{G(s)}{s}\right\}"
        substituted_latex=(
            r"\text{Equivalência por segurador de ordem zero (invariância ao degrau).}"
        )
        notes.append(
            "A invariância ao degrau foi calculada pela equivalência ZOH, correspondente à Eq. (5.8)."
        )
    else:
        _,sub,sub_latex,substituted,num_d,den_d=_symbolic_discretize(num,den,T,key)
        # A forma simbólica exata pode carregar enormes denominadores comuns
        # quando os parâmetros do controlador são decimais. Para apresentação,
        # mostramos a mesma FT já normalizada pelos coeficientes numéricos.
        substituted_latex=transfer_latex(num_d,den_d,z)
    zeros=[complex(v) for v in np.roots(num_d)] if len(num_d)>1 else []
    poles=[complex(v) for v in np.roots(den_d)] if len(den_d)>1 else []

    return DiscretizationResult(
        method=_METHOD_LABELS[key],
        sample_time=T,
        target_label=str(target_label),
        continuous_num=[float(v) for v in np.asarray(num,dtype=float)],
        continuous_den=[float(v) for v in np.asarray(den,dtype=float)],
        discrete_num=[float(v) for v in num_d],
        discrete_den=[float(v) for v in den_d],
        continuous_latex=transfer_latex(num,den,s),
        substitution_latex=sub_latex,
        substituted_latex=substituted_latex,
        discrete_latex=transfer_latex(num_d,den_d,z),
        difference_equation_latex=difference_equation_latex(num_d,den_d),
        poles=poles,
        zeros=zeros,
        notes=notes+["Coeficientes apresentados com o denominador normalizado para coeficiente líder igual a 1."],
    )


def multiply_transfers(*pairs):
    num=np.asarray([1.0])
    den=np.asarray([1.0])
    for n,d in pairs:
        num=np.polymul(num,np.asarray(n,dtype=float))
        den=np.polymul(den,np.asarray(d,dtype=float))
    return np.trim_zeros(num,"f").tolist(),np.trim_zeros(den,"f").tolist()


def project_discretization_target(design,data):
    target=str(data.get("discretization_target","Controlador Gc(s)"))
    if target=="Controlador Gc(s)":
        return design.controller_num,design.controller_den,"G_c(s)"
    if target=="G(s)":
        return data["num_g"],data["den_g"],"G(s)"
    if target=="G(s)H(s)":
        n,d=multiply_transfers(
            (data["num_g"],data["den_g"]),
            (data["num_h"],data["den_h"]),
        )
        return n,d,"G(s)H(s)"
    if target=="G(s)H(s)Gc(s)":
        n,d=multiply_transfers(
            (data["num_g"],data["den_g"]),
            (data["num_h"],data["den_h"]),
            (design.controller_num,design.controller_den),
        )
        return n,d,"G(s)H(s)G_c(s)"
    raise ValueError("Alvo de discretização desconhecido.")


def discretize_project(design,data):
    num,den,label=project_discretization_target(design,data)
    return discretize_transfer(
        num,den,
        data.get("sample_time",1.0),
        data.get("discretization_method","Tustin"),
        target_label=label,
    )
