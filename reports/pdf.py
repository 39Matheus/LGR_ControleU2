"""Relatórios compactos para reprodução manual."""
from io import BytesIO
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle,getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Paragraph,SimpleDocTemplate,Spacer,Table,TableStyle

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
    st += [Paragraph("1. Polo desejado",h),Paragraph(f"s_d={_fmtc(design.desired_pole)}; ξ={s.get('zeta',0):.6f}; ωn={s.get('omega_n',0):.6f}.",b)]
    st += [Paragraph("2. Condição de ângulo — trigonometria",h),
           Paragraph(f"Fase base={design.base_phase_deg:.6f}°. Fase necessária={design.required_phase_deg:.6f}°. Cada zero: φ={design.zero_angle_deg:.6f}°.",b),
           Paragraph(f"tan(φ)=ωd/(z−σ) ⇒ z={design.zero_parameter:.8g}.",b)]
    rows=[["Elemento","ΔRe","ΔIm","ângulo (°)"]]+[[c.label,f"{c.dx:.5f}",f"{c.dy:.5f}",f"{c.angle_deg:.5f}"] for c in design.plant_contributions+design.controller_contributions]
    t=Table(rows,hAlign="LEFT"); t.setStyle(TableStyle([("GRID",(0,0),(-1,-1),0.3,colors.grey),("BACKGROUND",(0,0),(-1,0),colors.lightgrey),("FONTSIZE",(0,0),(-1,-1),7.5)]))
    st += [t,Spacer(1,4),Paragraph("3. Condição de módulo",h),Paragraph(f"|GcG H|=1 ⇒ Kc={design.kc:.8g}.",b)]
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
