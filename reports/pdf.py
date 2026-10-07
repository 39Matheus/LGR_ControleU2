"""Relatórios compactos para reprodução manual."""
from io import BytesIO
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle,getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Paragraph,SimpleDocTemplate,Spacer,Table,TableStyle

from core.geometry import magnitude_breakdown

def _styles():
    s=getSampleStyleSheet()
    return (
        ParagraphStyle("T",parent=s["Title"],fontSize=16,alignment=TA_CENTER,spaceAfter=8),
        ParagraphStyle("H",parent=s["Heading2"],fontSize=10.5,spaceBefore=6,spaceAfter=3),
        ParagraphStyle("B",parent=s["BodyText"],fontSize=8.5,leading=10.5,spaceAfter=3)
    )

def _fmtc(z,d=5):
    z=complex(z)
    if abs(z.imag)<10**(-d): return f"{z.real:.{d}f}"
    return f"{z.real:.{d}f} {'+' if z.imag>=0 else '-'} {abs(z.imag):.{d}f}j"

def build_controller_pdf(design,num_g,den_g,num_h,den_h):
    buf=BytesIO(); doc=SimpleDocTemplate(buf,pagesize=A4,leftMargin=1.35*cm,rightMargin=1.35*cm,topMargin=1.2*cm,bottomMargin=1.2*cm)
    title,h,b=_styles()
    st=[Paragraph(f"Projeto {design.controller_type} pelo LGR — Modo Prova",title)]
    s=design.specification_summary
    source=s.get("source")
    initial_pole=complex(s.get("initial_pole",design.desired_pole))
    zeta0=float(s.get("initial_zeta",s.get("zeta",0.0)))
    wn0=float(s.get("initial_omega_n",s.get("omega_n",0.0)))
    sigma0=float(s.get("initial_sigma",-initial_pole.real))
    wd0=float(s.get("initial_omega_d",abs(initial_pole.imag)))

    st.append(Paragraph("1. Obtenção do polo desejado",h))
    if source=="mp_ts":
        mp=float(s["mp_max"]); ts=float(s["ts_max"]); band=float(s["settling_band"])
        c=3.0 if abs(band-0.05)<1e-12 else 4.0
        st += [
            Paragraph(
                f"Mp={mp:.6g}% -> Mp={mp/100.0:.6g}; "
                f"xi=-ln(Mp)/sqrt(pi^2+ln(Mp)^2)={zeta0:.6f}.",
                b,
            ),
            Paragraph(
                f"ts({100*band:.0f}%)≈{c:g}/sigma -> "
                f"sigma={c:g}/{ts:.6g}={sigma0:.6f}.",
                b,
            ),
            Paragraph(
                f"wn=sigma/xi={wn0:.6f}; "
                f"wd=wn*sqrt(1-xi^2)={wd0:.6f}; "
                f"s_d0={_fmtc(initial_pole)}.",
                b,
            ),
        ]
        if s.get("adjusted"):
            st.append(Paragraph(
                f"Após verificação da planta completa, mantendo xi e aumentando sigma: "
                f"s_d={_fmtc(design.desired_pole)}; "
                f"sigma_f={s['sigma']:.6f}; wn_f={s['omega_n']:.6f}.",
                b,
            ))
    elif source=="zeta_wn":
        st += [
            Paragraph(
                f"Dados xi={zeta0:.6f} e wn={wn0:.6f}: "
                f"sigma=xi*wn={sigma0:.6f}.",
                b,
            ),
            Paragraph(
                f"wd=wn*sqrt(1-xi^2)={wd0:.6f}; "
                f"s_d=-sigma±jwd={_fmtc(initial_pole)} e conjugado.",
                b,
            ),
        ]
    else:
        st += [
            Paragraph(
                f"Polo fornecido: s_d={_fmtc(initial_pole)} e conjugado; "
                f"wn=|s_d|={wn0:.6f}; xi=-Re(s_d)/wn={zeta0:.6f}.",
                b,
            )
        ]

    plant_zeros=[c for c in design.plant_contributions if c.kind=="zero"]
    plant_poles=[c for c in design.plant_contributions if c.kind=="pole"]
    controller_poles=[c for c in design.controller_contributions if c.kind=="pole"]
    controller_zeros=[c for c in design.controller_contributions if c.kind=="zero"]
    kg=float(num_g[0]/den_g[0]); kh=float(num_h[0]/den_h[0])
    gain_phase=180.0 if kg*kh<0 else 0.0
    raw_base=(
        gain_phase
        +sum(c.angle_deg for c in plant_zeros)
        -sum(c.angle_deg for c in plant_poles)
        -sum(c.angle_deg for c in controller_poles)
    )
    ctrl_phase=sum(c.angle_deg for c in controller_zeros)
    final_phase=raw_base+ctrl_phase
    target=180.0+360.0*round((final_phase-180.0)/360.0)

    st += [
        Paragraph("2. Condição de ângulo — trigonometria",h),
        Paragraph(
            "Convenção: soma(phi_zeros) - soma(theta_polos) = (2q+1)180 graus.",
            b,
        ),
    ]
    rows=[["Elemento","tipo","dRe","dIm","ângulo (graus)"]]+[
        [c.label,c.kind,f"{c.dx:.5f}",f"{c.dy:.5f}",f"{c.angle_deg:.5f}"]
        for c in design.plant_contributions+design.controller_contributions
    ]
    t=Table(rows,hAlign="LEFT")
    t.setStyle(TableStyle([
        ("GRID",(0,0),(-1,-1),0.3,colors.grey),
        ("BACKGROUND",(0,0),(-1,0),colors.lightgrey),
        ("FONTSIZE",(0,0),(-1,-1),7.5),
    ]))
    st += [t,Spacer(1,4)]
    st.append(Paragraph(
        f"Fase sem os zeros do controlador = {raw_base:.5f} graus. "
        f"Contribuição total necessária dos zeros do controlador = "
        f"{target-raw_base:.5f} graus.",
        b,
    ))
    if len(controller_zeros)==1:
        st.append(Paragraph(
            f"{raw_base:.5f} + phi_c = {target:.0f} -> "
            f"phi_c={design.zero_angle_deg:.5f} graus.",
            b,
        ))
    else:
        st.append(Paragraph(
            f"{raw_base:.5f} + {len(controller_zeros)} phi_c = {target:.0f} -> "
            f"phi_c={design.zero_angle_deg:.5f} graus.",
            b,
        ))
    st.append(Paragraph(
        f"tan(phi_c)=wd/(zc-sigma) -> "
        f"zc=sigma+wd/tan(phi_c)={design.zero_parameter:.8g}.",
        b,
    ))
    st.append(Paragraph(
        f"Conferência da fase: {final_phase:.5f} graus, equivalente a "
        "180 graus módulo 360.",
        b,
    ))

    mag=magnitude_breakdown(design,num_g,den_g,num_h,den_h)
    st.append(Paragraph("3. Condição de módulo",h))
    st.append(Paragraph("|Gc(sd)G(sd)H(sd)|=1.",b))
    mag_rows=[["Dist.","Elemento","dRe","dIm","módulo"]]
    for item in mag["poles"]+mag["zeros"]:
        mag_rows.append([
            item["term"],item["source_label"],
            f"{item['dx']:.5f}",f"{item['dy']:.5f}",f"{item['distance']:.6f}",
        ])
    mt=Table(mag_rows,hAlign="LEFT")
    mt.setStyle(TableStyle([
        ("GRID",(0,0),(-1,-1),0.3,colors.grey),
        ("BACKGROUND",(0,0),(-1,0),colors.lightgrey),
        ("FONTSIZE",(0,0),(-1,-1),7.5),
    ]))
    st += [mt,Spacer(1,4)]
    for item in mag["poles"]+mag["zeros"]:
        st.append(Paragraph(
            f"{item['term']}=sqrt(({item['dx']:.6f})^2+"
            f"({item['dy']:.6f})^2)={item['distance']:.6f}.",
            b,
        ))
    a_names=".".join(x["term"] for x in mag["poles"]) or "1"
    b_names=".".join(x["term"] for x in mag["zeros"]) or "1"
    st += [
        Paragraph(
            f"KT=prod(Ai)/prod(Bi)=({a_names})/({b_names})="
            f"{mag['prod_a']:.8g}/{mag['prod_b']:.8g}={mag['kt']:.8g}.",
            b,
        ),
        Paragraph(
            f"KG={mag['kg']:.8g}; KH={mag['kh']:.8g}; "
            "KT=|KG KH| KC.",
            b,
        ),
        Paragraph(
            f"KC=KT/|KG KH|={mag['kt']:.8g}/"
            f"|({mag['kg']:.8g})({mag['kh']:.8g})|={mag['kc']:.8g}.",
            b,
        ),
    ]
    z=design.zero_parameter
    if design.controller_type=="PD": ctrl=f"Gc(s)={design.kc:.8g}(s+{z:.8g}); Kp={design.kp:.8g}; Kd={design.kd:.8g}"
    elif design.controller_type=="PI": ctrl=f"Gc(s)={design.kc:.8g}(s+{z:.8g})/s; Kp={design.kp:.8g}; Ki={design.ki:.8g}"
    else: ctrl=f"Gc(s)={design.kc:.8g}(s+{z:.8g})²/s; Kp={design.kp:.8g}; Ki={design.ki:.8g}; Kd={design.kd:.8g}"
    st += [Paragraph("4. Controlador",h),Paragraph(ctrl,b),Paragraph("5. Verificação",h),Paragraph("Polos MF: "+"; ".join(_fmtc(p) for p in design.closed_loop_poles),b)]
    if design.metrics:
        m=design.metrics; tx=f"Estável={m.stable}"
        if m.overshoot_percent is not None: tx+=f"; Mp={m.overshoot_percent:.4f}%"
        if m.settling_time is not None: tx+=f"; ts({100*m.settling_band:.0f}%)={m.settling_time:.5f}s"
        st.append(Paragraph(tx,b))
    doc.build(st); return buf.getvalue()

def build_lgr_pdf(a,p1,p4,p7,p8,p9,p10,ptest,point):
    buf=BytesIO(); doc=SimpleDocTemplate(buf,pagesize=A4,leftMargin=1.35*cm,rightMargin=1.35*cm,topMargin=1.2*cm,bottomMargin=1.2*cm)
    title,h,b=_styles()
    st=[Paragraph("LGR — Modo Prova",title),Paragraph("1–2. Polinômio, polos e zeros",h),Paragraph(f"Φ(s,K)={p1['char_expr']}=0",b),
        Paragraph("Polos: "+", ".join(_fmtc(p) for p in a.polos),b),Paragraph("Zeros: "+(", ".join(_fmtc(z) for z in a.zeros) or "nenhum"),b),
        Paragraph("4. Segmentos reais",h),Paragraph("; ".join(str(x) for x in p4["segmentos"]) or "nenhum",b),
        Paragraph("7. Assíntotas",h),Paragraph(f"centro={p7.get('centro')}; ângulos={p7.get('angulos',[])}",b),
        Paragraph("9. Cruzamento imaginário",h),Paragraph("; ".join(f"s=±j{c['w']:.5f}, K={c['K']:.5f}" for c in p9["cruzamentos"]) or "nenhum",b),
        Paragraph("11–12. Ponto teste",h),Paragraph(f"s_t={_fmtc(point)}; fase={ptest['fase_mod']:.5f}°; pertence={ptest['pertence']}; K={ptest['K']}",b)]
    doc.build(st); return buf.getvalue()
