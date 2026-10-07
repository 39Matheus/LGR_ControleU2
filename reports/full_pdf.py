"""Relatório completo em PDF com a mesma sequência da interface."""
from io import BytesIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    Image, PageBreak, Paragraph, SimpleDocTemplate,
    Spacer, Table, TableStyle,
)

from core.geometry import angle_breakdown, magnitude_breakdown
from core.visualization import geometry_figure, root_locus_figure, step_response_figure


def _fmtc(z,d=6):
    z=complex(z)
    if abs(z.imag)<10**(-d):
        return f"{z.real:.{d}f}"
    return f"{z.real:.{d}f} {'+' if z.imag>=0 else '-'} {abs(z.imag):.{d}f}j"


def _styles():
    styles=getSampleStyleSheet()
    return {
        "title":ParagraphStyle(
            "FullTitle",parent=styles["Title"],fontSize=18,leading=21,
            alignment=TA_CENTER,spaceAfter=10,textColor=colors.HexColor("#111111")
        ),
        "h2":ParagraphStyle(
            "FullH2",parent=styles["Heading2"],fontSize=12,leading=15,
            spaceBefore=10,spaceAfter=6,textColor=colors.HexColor("#111111")
        ),
        "h3":ParagraphStyle(
            "FullH3",parent=styles["Heading3"],fontSize=10,leading=12,
            spaceBefore=6,spaceAfter=4,textColor=colors.HexColor("#222222")
        ),
        "body":ParagraphStyle(
            "FullBody",parent=styles["BodyText"],fontSize=8.5,leading=11,
            spaceAfter=4,textColor=colors.HexColor("#222222")
        ),
    }


def _eq(tex,max_width=17.2*cm,font_size=13):
    """Renderiza uma equação MathText como imagem para o ReportLab."""
    buf=BytesIO()
    fig=plt.figure(figsize=(10,0.55))
    fig.patch.set_facecolor("white")
    fig.text(0.5,0.5,"$"+tex+"$",ha="center",va="center",fontsize=font_size,color="black")
    plt.axis("off")
    fig.savefig(buf,format="png",dpi=180,bbox_inches="tight",pad_inches=0.08,facecolor="white")
    plt.close(fig)
    buf.seek(0)
    img=Image(buf)
    scale=min(1.0,max_width/img.imageWidth)
    img.drawWidth=img.imageWidth*scale
    img.drawHeight=img.imageHeight*scale
    img.hAlign="CENTER"
    return img


def _table(headers,rows,widths=None,font_size=7.2):
    t=Table([headers]+rows,colWidths=widths,repeatRows=1,hAlign="LEFT")
    t.setStyle(TableStyle([
        ("GRID",(0,0),(-1,-1),0.35,colors.HexColor("#999999")),
        ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#eeeeee")),
        ("TEXTCOLOR",(0,0),(-1,-1),colors.HexColor("#111111")),
        ("FONTSIZE",(0,0),(-1,-1),font_size),
        ("VALIGN",(0,0),(-1,-1),"MIDDLE"),
        ("LEFTPADDING",(0,0),(-1,-1),4),
        ("RIGHTPADDING",(0,0),(-1,-1),4),
        ("TOPPADDING",(0,0),(-1,-1),3),
        ("BOTTOMPADDING",(0,0),(-1,-1),3),
    ]))
    return t


def _mpl_marker(symbol):
    return {
        "x":"x","circle-open":"o","circle":"o","star":"*",
        "diamond":"D","square":"s","triangle-up":"^",
    }.get(str(symbol),"o")


def _mpl_dash(dash):
    return {
        "dot":":","dash":"--","dashdot":"-.","longdash":"--","solid":"-",
    }.get(str(dash or "solid"),"-")


def _plotly_png(fig,width=7.7,height=4.6,dpi=150):
    """Renderiza os traces Plotly em Matplotlib para um PDF portátil."""
    out=BytesIO()
    mf,ax=plt.subplots(figsize=(width,height))
    mf.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for tr in fig.data:
        x=getattr(tr,"x",None)
        y=getattr(tr,"y",None)
        if x is None or y is None:
            continue
        try:
            xf=np.asarray(x,dtype=float)
            yf=np.asarray(y,dtype=float)
        except Exception:
            continue
        if xf.size==0 or yf.size==0:
            continue

        mode=str(getattr(tr,"mode","lines") or "lines")
        name=str(getattr(tr,"name","") or "")
        line=getattr(tr,"line",None)
        marker=getattr(tr,"marker",None)
        color=None
        if line is not None and getattr(line,"color",None):
            color=line.color
        elif marker is not None and isinstance(getattr(marker,"color",None),str):
            color=marker.color

        if "lines" in mode:
            ax.plot(
                xf,yf,label=name or None,
                linewidth=float(getattr(line,"width",1.4) or 1.4),
                linestyle=_mpl_dash(getattr(line,"dash",None)),
                color=color,
            )
        if "markers" in mode:
            symbol=getattr(marker,"symbol","circle") if marker is not None else "circle"
            size=float(getattr(marker,"size",7) or 7)
            face="none" if "open" in str(symbol) else color
            ax.scatter(
                xf,yf,label=(name or None) if "lines" not in mode else None,
                s=max(16,size**2),marker=_mpl_marker(symbol),
                facecolors=face,edgecolors=color if color else None,linewidths=1.1,
            )
        if "text" in mode:
            texts=getattr(tr,"text",None)
            if texts is not None:
                for xx,yy,tt in zip(xf,yf,texts):
                    if tt:
                        ax.annotate(str(tt),(xx,yy),xytext=(3,4),textcoords="offset points",fontsize=6)

    shapes=getattr(fig.layout,"shapes",None) or []
    for sh in shapes:
        try:
            if sh.type=="line":
                if sh.x0==sh.x1:
                    ax.axvline(float(sh.x0),linewidth=float(sh.line.width or 1),linestyle=_mpl_dash(sh.line.dash),color=sh.line.color or "#777777",alpha=float(sh.opacity or 1))
                elif sh.y0==sh.y1:
                    ax.axhline(float(sh.y0),linewidth=float(sh.line.width or 1),linestyle=_mpl_dash(sh.line.dash),color=sh.line.color or "#777777",alpha=float(sh.opacity or 1))
        except Exception:
            pass

    title=getattr(getattr(fig.layout,"title",None),"text",None)
    if title:
        ax.set_title(str(title),fontsize=10)
    xt=getattr(getattr(getattr(fig.layout,"xaxis",None),"title",None),"text",None)
    yt=getattr(getattr(getattr(fig.layout,"yaxis",None),"title",None),"text",None)
    if xt: ax.set_xlabel(str(xt),fontsize=8)
    if yt: ax.set_ylabel(str(yt),fontsize=8)

    xr=getattr(getattr(fig.layout,"xaxis",None),"range",None)
    yr=getattr(getattr(fig.layout,"yaxis",None),"range",None)
    if xr and len(xr)==2:
        ax.set_xlim(float(xr[0]),float(xr[1]))
    if yr and len(yr)==2:
        ax.set_ylim(float(yr[0]),float(yr[1]))

    ax.grid(True,alpha=.18)
    ax.tick_params(labelsize=7)
    handles,labels=ax.get_legend_handles_labels()
    if labels:
        unique={}
        for h,l in zip(handles,labels):
            if l and l not in unique:
                unique[l]=h
        ax.legend(unique.values(),unique.keys(),fontsize=5.7,loc="best",frameon=False)

    mf.tight_layout()
    mf.savefig(out,format="png",dpi=dpi,bbox_inches="tight",facecolor="white")
    plt.close(mf)
    out.seek(0)
    return out


def _chart_image(fig,max_width,max_height=None,width=7.7,height=4.6):
    img=Image(_plotly_png(fig,width=width,height=height))
    scale=max_width/img.imageWidth
    if max_height is not None:
        scale=min(scale,max_height/img.imageHeight)
    img.drawWidth=img.imageWidth*scale
    img.drawHeight=img.imageHeight*scale
    img.hAlign="CENTER"
    return img


def _spec_story(design,styles):
    s=design.specification_summary
    source=s.get("source")
    initial=complex(s.get("initial_pole",design.desired_pole))
    zeta=float(s.get("initial_zeta",s.get("zeta",0.0)))
    wn=float(s.get("initial_omega_n",s.get("omega_n",0.0)))
    sigma=float(s.get("initial_sigma",-initial.real))
    wd=float(s.get("initial_omega_d",abs(initial.imag)))
    story=[Paragraph("1 - Obtenção do polo desejado",styles["h2"])]

    if source=="mp_ts":
        mp=float(s["mp_max"]); ts=float(s["ts_max"]); band=float(s["settling_band"])
        c=3.0 if abs(band-0.05)<1e-12 else 4.0
        story += [
            Paragraph("a) Amortecimento a partir do máximo sobresinal",styles["h3"]),
            _eq(rf"M_p={mp:.6g}\%={mp/100.0:.6g},\quad \xi=\frac{{-\ln(M_p)}}{{\sqrt{{\pi^2+\ln^2(M_p)}}}}={zeta:.6f}"),
            Paragraph("b) Parte real a partir do tempo de acomodação",styles["h3"]),
            _eq(rf"t_s({100*band:.0f}\%)\approx\frac{{{c:g}}}{{\xi\omega_n}}=\frac{{{c:g}}}{{\sigma}},\quad \sigma=\frac{{{c:g}}}{{{ts:.6g}}}={sigma:.6f}"),
            Paragraph("c) Frequências e polo de fronteira",styles["h3"]),
            _eq(rf"\omega_n=\frac{{\sigma}}{{\xi}}={wn:.6f},\quad \omega_d=\omega_n\sqrt{{1-\xi^2}}={wd:.6f}"),
            _eq(rf"s_{{d,0}}=-\sigma\pm j\omega_d={initial.real:.6f}\pm j{abs(initial.imag):.6f}"),
        ]
        if s.get("adjusted"):
            final=complex(design.desired_pole)
            story.append(Paragraph(
                "Ajuste fino habilitado: após a verificação da planta completa, o polo foi deslocado mantendo xi.",
                styles["body"],
            ))
            story.append(_eq(rf"s_d={final.real:.6f}\pm j{abs(final.imag):.6f}"))
        else:
            story.append(Paragraph(
                "Ajuste fino desabilitado: o projeto usa diretamente o polo de fronteira.",
                styles["body"],
            ))
    elif source=="zeta_wn":
        story += [
            Paragraph("Como xi e omega_n são fornecidos, usa-se diretamente o modelo dominante de 2ª ordem.",styles["body"]),
            _eq(rf"\sigma=\xi\omega_n={zeta:.6f}\cdot {wn:.6f}={sigma:.6f}"),
            _eq(rf"\omega_d=\omega_n\sqrt{{1-\xi^2}}={wd:.6f}"),
            _eq(rf"s_d=-\xi\omega_n\pm j\omega_n\sqrt{{1-\xi^2}}={initial.real:.6f}\pm j{abs(initial.imag):.6f}"),
        ]
    else:
        story += [
            Paragraph("Os polos desejados foram fornecidos diretamente.",styles["body"]),
            _eq(rf"s_d={initial.real:.6f}\pm j{abs(initial.imag):.6f},\quad \omega_n=|s_d|={wn:.6f},\quad \xi=\frac{{-\mathrm{{Re}}(s_d)}}{{\omega_n}}={zeta:.6f}"),
        ]
    return story


def build_controller_full_pdf(design,data):
    buf=BytesIO()
    styles=_styles()
    doc=SimpleDocTemplate(
        buf,pagesize=A4,leftMargin=1.25*cm,rightMargin=1.25*cm,
        topMargin=1.15*cm,bottomMargin=1.15*cm,
        title=f"Resolução completa {design.controller_type} - LGR",
    )
    story=[
        Paragraph(f"Projeto {design.controller_type} pelo LGR - Resolução completa",styles["title"]),
        Paragraph(
            "Relatório detalhado com o mesmo encadeamento da interface. "
            "Os gráficos são incluídos como imagens estáticas no PDF.",
            styles["body"],
        ),
    ]
    story += _spec_story(design,styles)

    info=angle_breakdown(design,data["num_g"],data["den_g"],data["num_h"],data["den_h"])
    allc=info["plant_zeros"]+info["plant_poles"]+info["controller_poles"]+info["controller_zeros"]
    story.append(Paragraph("2 - Condição de ângulo",styles["h2"]))
    story.append(Paragraph(
        "Convenção: soma dos ângulos dos zeros menos soma dos ângulos dos polos.",
        styles["body"],
    ))
    story.append(_eq(r"\sum\phi_{\mathrm{zeros}}-\sum\theta_{\mathrm{polos}}=(2q+1)180^\circ"))
    story.append(_table(
        ["Elemento","Tipo","Singularidade","dRe","dIm","Ângulo"],
        [[c.label,c.kind,_fmtc(c.singularity,5),f"{c.dx:.5f}",f"{c.dy:.5f}",f"{c.angle_deg:.5f}°"] for c in allc],
        widths=[2.0*cm,1.8*cm,3.4*cm,2.5*cm,2.5*cm,2.8*cm],
    ))
    story.append(Spacer(1,5))

    zero_sum=" + ".join(f"{c.angle_deg:.5f}" for c in info["plant_zeros"]) or "0"
    pole_sum=" + ".join(f"{c.angle_deg:.5f}" for c in info["plant_poles"]+info["controller_poles"]) or "0"
    if info["constant_phase"]:
        story.append(_eq(rf"{info['constant_phase']:.5f}+({zero_sum})-({pole_sum})={info['raw_base']:.5f}^\circ"))
    else:
        story.append(_eq(rf"({zero_sum})-({pole_sum})={info['raw_base']:.5f}^\circ"))
    ncz=len(info["controller_zeros"])
    if ncz==1:
        story += [
            _eq(rf"{info['raw_base']:.5f}^\circ+\phi_c={info['target']:.0f}^\circ"),
            _eq(rf"\phi_c={design.zero_angle_deg:.5f}^\circ"),
        ]
    else:
        story += [
            _eq(rf"{info['raw_base']:.5f}^\circ+{ncz}\phi_c={info['target']:.0f}^\circ"),
            _eq(rf"\phi_c={design.zero_angle_deg:.5f}^\circ"),
        ]
    sigma=-design.desired_pole.real
    wd=abs(design.desired_pole.imag)
    story += [
        _eq(r"\tan(\phi_c)=\frac{\omega_d}{z_c-\sigma}"),
        _eq(rf"z_c=\sigma+\frac{{\omega_d}}{{\tan(\phi_c)}}={sigma:.6f}+\frac{{{wd:.6f}}}{{\tan({design.zero_angle_deg:.5f}^\circ)}}={design.zero_parameter:.6f}"),
        Paragraph("Zero(s) do controlador: "+", ".join(f"s={z:.6f}" for z in design.zero_locations),styles["body"]),
    ]

    mag=magnitude_breakdown(design,data["num_g"],data["den_g"],data["num_h"],data["den_h"])
    story.append(Paragraph("3 - Condição de módulo",styles["h2"]))
    story.append(_eq(r"|G_c(s_d)G(s_d)H(s_d)|=1"))
    story.append(Paragraph(
        "A_i são as distâncias do polo desejado aos polos; B_i são as distâncias aos zeros.",
        styles["body"],
    ))
    story.append(_table(
        ["Dist.","Elemento","Singularidade","dRe","dIm","Módulo"],
        [[x["term"],x["source_label"],_fmtc(x["singularity"],5),f"{x['dx']:.5f}",f"{x['dy']:.5f}",f"{x['distance']:.6f}"] for x in mag["poles"]+mag["zeros"]],
        widths=[1.6*cm,2.0*cm,3.4*cm,2.5*cm,2.5*cm,3.0*cm],
    ))
    story.append(Spacer(1,5))
    for x in mag["poles"]+mag["zeros"]:
        story.append(_eq(rf"{x['term']}=\sqrt{{({x['dx']:.6f})^2+({x['dy']:.6f})^2}}={x['distance']:.6f}",font_size=11.5))
    a_prod=r"\,\cdot\,".join(x["term"] for x in mag["poles"]) or "1"
    b_prod=r"\,\cdot\,".join(x["term"] for x in mag["zeros"]) or "1"
    story += [
        _eq(rf"K_T=\frac{{{a_prod}}}{{{b_prod}}}=\frac{{{mag['prod_a']:.8g}}}{{{mag['prod_b']:.8g}}}={mag['kt']:.8g}"),
        _eq(rf"K_G={mag['kg']:.8g},\quad K_H={mag['kh']:.8g},\quad K_T=|K_GK_H|K_C"),
        _eq(rf"K_C=\frac{{K_T}}{{|K_GK_H|}}=\frac{{{mag['kt']:.8g}}}{{|({mag['kg']:.8g})({mag['kh']:.8g})|}}={mag['kc']:.8g}"),
    ]

    story.append(Paragraph("4 - Controlador",styles["h2"]))
    z=design.zero_parameter
    if design.controller_type=="PD":
        story += [_eq(rf"G_c(s)={design.kc:.8g}(s+{z:.8g})"),_eq(rf"K_p={design.kp:.8g},\quad K_d={design.kd:.8g}")]
    elif design.controller_type=="PI":
        story += [_eq(rf"G_c(s)={design.kc:.8g}\frac{{s+{z:.8g}}}{{s}}"),_eq(rf"K_p={design.kp:.8g},\quad K_i={design.ki:.8g}")]
    else:
        story += [_eq(rf"G_c(s)={design.kc:.8g}\frac{{(s+{z:.8g})^2}}{{s}}"),_eq(rf"K_p={design.kp:.8g},\quad K_i={design.ki:.8g},\quad K_d={design.kd:.8g}")]

    story.append(PageBreak())
    story.append(Paragraph("5 - Representação visual do projeto",styles["h2"]))
    geometry=geometry_figure(design)
    before=root_locus_figure(
        data["num_g"],data["den_g"],data["num_h"],data["den_h"],
        desired_pole=design.desired_pole,title="LGR antes do controlador",
    )
    after=root_locus_figure(
        data["num_g"],data["den_g"],data["num_h"],data["den_h"],
        desired_pole=design.desired_pole,
        controller_num=design.controller_num,controller_den=design.controller_den,
        title="LGR com o controlador projetado",
    )
    story.append(_chart_image(geometry,17.2*cm,9.2*cm,width=8.0,height=4.3))
    pair=Table(
        [[_chart_image(before,8.25*cm,7.2*cm,width=4.3,height=4.0),
          _chart_image(after,8.25*cm,7.2*cm,width=4.3,height=4.0)]],
        colWidths=[8.4*cm,8.4*cm],
    )
    pair.setStyle(TableStyle([
        ("VALIGN",(0,0),(-1,-1),"TOP"),
        ("LEFTPADDING",(0,0),(-1,-1),0),
        ("RIGHTPADDING",(0,0),(-1,-1),0),
    ]))
    story.append(pair)
    step=step_response_figure(data["num_g"],data["den_g"],data["num_h"],data["den_h"],design)
    if step is not None:
        story.append(_chart_image(step,17.2*cm,8.3*cm,width=8.0,height=4.0))

    story.append(Paragraph("6 - Verificação numérica",styles["h2"]))
    story.append(Paragraph("Polos de malha fechada",styles["h3"]))
    story.append(_table(
        ["Polo"],
        [[_fmtc(p,6)] for p in sorted(design.closed_loop_poles,key=lambda x:(x.real,x.imag))],
        widths=[16.7*cm],
    ))
    if design.metrics:
        m=design.metrics
        story.append(Spacer(1,5))
        story.append(_table(
            ["Métrica","Valor"],
            [["Estável","Sim" if m.stable else "Não"],
             ["Mp","-" if m.overshoot_percent is None else f"{m.overshoot_percent:.4f}%"],
             [f"ts ({100*m.settling_band:.0f}%)","-" if m.settling_time is None else f"{m.settling_time:.6f} s"]],
            widths=[8.35*cm,8.35*cm],
        ))
    if len(design.refinement_history)>1:
        story.append(Paragraph(
            f"Ajuste fino aplicado em {len(design.refinement_history)-1} iteração(ões), mantendo xi.",
            styles["body"],
        ))

    doc.build(story)
    return buf.getvalue()
