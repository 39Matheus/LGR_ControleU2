"""Relatório completo em HTML, preservando equações e gráficos Plotly."""
from html import escape
from io import StringIO

import plotly.io as pio

from core.geometry import angle_breakdown,magnitude_breakdown
from core.visualization import focus_description,geometry_figure,root_locus_figure,step_response_figure


def _fmtc(z,d=6):
    z=complex(z)
    if abs(z.imag)<10**(-d):
        return f"{z.real:.{d}f}"
    return f"{z.real:.{d}f} {'+' if z.imag>=0 else '-'} {abs(z.imag):.{d}f}j"


def _table(headers,rows):
    out=["<div class='table-wrap'><table><thead><tr>"]
    out += [f"<th>{escape(str(h))}</th>" for h in headers]
    out.append("</tr></thead><tbody>")
    for row in rows:
        out.append("<tr>")
        out += [f"<td>{escape(str(v))}</td>" for v in row]
        out.append("</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def _math(expr):
    return f"<div class='math'>\\[{expr}\\]</div>"


def _spec_section(design):
    s=design.specification_summary
    source=s.get("source")
    initial=complex(s.get("initial_pole",design.desired_pole))
    zeta=float(s.get("initial_zeta",s.get("zeta",0.0)))
    wn=float(s.get("initial_omega_n",s.get("omega_n",0.0)))
    sigma=float(s.get("initial_sigma",-initial.real))
    wd=float(s.get("initial_omega_d",abs(initial.imag)))

    parts=["<section><h2>1 — Obtenção do polo desejado</h2>"]
    if source=="mp_ts":
        mp=float(s["mp_max"])
        ts=float(s["ts_max"])
        band=float(s["settling_band"])
        c=3.0 if abs(band-0.05)<1e-12 else 4.0
        parts.append("<h3>a) Amortecimento a partir do máximo sobresinal</h3>")
        parts.append(_math(
            rf"M_p={mp:.6g}\%={mp/100.0:.6g},\qquad "
            rf"\xi=\frac{{-\ln(M_p)}}{{\sqrt{{\pi^2+\ln^2(M_p)}}}}={zeta:.6f}"
        ))
        parts.append("<h3>b) Parte real a partir do tempo de acomodação</h3>")
        parts.append(_math(
            rf"t_s({100*band:.0f}\%)\approx\frac{{{c:g}}}{{\xi\omega_n}}"
            rf"=\frac{{{c:g}}}{{\sigma}},\qquad "
            rf"\sigma=\frac{{{c:g}}}{{{ts:.6g}}}={sigma:.6f}"
        ))
        parts.append("<h3>c) Frequências e polo de fronteira</h3>")
        parts.append(_math(
            rf"\omega_n=\frac{{\sigma}}{{\xi}}={wn:.6f},\qquad "
            rf"\omega_d=\omega_n\sqrt{{1-\xi^2}}={wd:.6f}"
        ))
        parts.append(_math(
            rf"s_{{d,0}}=-\sigma\pm j\omega_d={initial.real:.6f}\pm j{abs(initial.imag):.6f}"
        ))
        if s.get("adjusted"):
            final=complex(design.desired_pole)
            parts.append(
                "<div class='note'>Ajuste fino habilitado: após verificar a planta completa, "
                "o polo foi deslocado mantendo ξ.</div>"
            )
            parts.append(_math(
                rf"s_d={final.real:.6f}\pm j{abs(final.imag):.6f}"
            ))
        else:
            parts.append(
                "<div class='ok'>Ajuste fino desabilitado: o projeto usa diretamente o polo de fronteira.</div>"
            )
    elif source=="zeta_wn":
        parts.append(_math(
            rf"\sigma=\xi\omega_n={zeta:.6f}\cdot {wn:.6f}={sigma:.6f}"
        ))
        parts.append(_math(
            rf"\omega_d=\omega_n\sqrt{{1-\xi^2}}={wd:.6f}"
        ))
        parts.append(_math(
            rf"s_d=-\xi\omega_n\pm j\omega_n\sqrt{{1-\xi^2}}"
            rf"={initial.real:.6f}\pm j{abs(initial.imag):.6f}"
        ))
    else:
        parts.append("<p>Os polos desejados foram fornecidos diretamente.</p>")
        parts.append(_math(
            rf"s_d={initial.real:.6f}\pm j{abs(initial.imag):.6f},\qquad "
            rf"\omega_n=|s_d|={wn:.6f},\qquad "
            rf"\xi=\frac{{-\operatorname{{Re}}(s_d)}}{{\omega_n}}={zeta:.6f}"
        ))
    parts.append("</section>")
    return "".join(parts)


def _angle_section(design,data):
    info=angle_breakdown(
        design,data["num_g"],data["den_g"],data["num_h"],data["den_h"]
    )
    allc=(
        info["plant_zeros"]+info["plant_poles"]
        +info["controller_poles"]+info["controller_zeros"]
    )
    rows=[]
    for c in allc:
        rows.append([
            c.label,c.kind,_fmtc(c.singularity,5),
            f"{c.dx:.6f}",f"{c.dy:.6f}",f"{c.angle_deg:.6f}°",
        ])

    zero_sum=" + ".join(f"{c.angle_deg:.5f}" for c in info["plant_zeros"]) or "0"
    pole_sum=" + ".join(
        f"{c.angle_deg:.5f}" for c in info["plant_poles"]+info["controller_poles"]
    ) or "0"
    ncz=len(info["controller_zeros"])

    parts=[
        "<section><h2>2 — Condição de ângulo</h2>",
        "<p>Convenção usada: soma dos ângulos dos zeros menos soma dos ângulos dos polos.</p>",
        _math(r"\sum\phi_{\mathrm{zeros}}-\sum\theta_{\mathrm{polos}}=(2q+1)180^\circ"),
        _table(
            ["Elemento","Tipo","Singularidade","ΔRe","ΔIm","Ângulo"],
            rows,
        ),
    ]

    if info["constant_phase"]:
        parts.append(_math(
            rf"{info['constant_phase']:.5f}+({zero_sum})-({pole_sum})"
            rf"={info['raw_base']:.5f}^\circ"
        ))
    else:
        parts.append(_math(
            rf"({zero_sum})-({pole_sum})={info['raw_base']:.5f}^\circ"
        ))

    if ncz==1:
        parts.append(_math(
            rf"{info['raw_base']:.5f}^\circ+\phi_c={info['target']:.0f}^\circ"
        ))
        parts.append(_math(
            rf"\phi_c={design.zero_angle_deg:.5f}^\circ"
        ))
    else:
        parts.append(_math(
            rf"{info['raw_base']:.5f}^\circ+{ncz}\phi_c={info['target']:.0f}^\circ"
        ))
        parts.append(_math(
            rf"\phi_c={design.zero_angle_deg:.5f}^\circ"
        ))

    sigma=-design.desired_pole.real
    wd=abs(design.desired_pole.imag)
    parts.append(_math(
        rf"\tan(\phi_c)=\frac{{\omega_d}}{{z_c-\sigma}}"
    ))
    parts.append(_math(
        rf"z_c=\sigma+\frac{{\omega_d}}{{\tan(\phi_c)}}"
        rf"={sigma:.6f}+\frac{{{wd:.6f}}}{{\tan({design.zero_angle_deg:.5f}^\circ)}}"
        rf"={design.zero_parameter:.6f}"
    ))
    parts.append(
        f"<div class='ok'>Zero(s) do controlador: "
        +", ".join(escape(f"s={z:.6f}") for z in design.zero_locations)
        +"</div>"
    )
    parts.append("</section>")
    return "".join(parts)


def _magnitude_section(design,data):
    m=magnitude_breakdown(
        design,data["num_g"],data["den_g"],data["num_h"],data["den_h"]
    )
    rows=[]
    for item in m["poles"]+m["zeros"]:
        rows.append([
            item["term"],item["source_label"],_fmtc(item["singularity"],5),
            f"{item['dx']:.6f}",f"{item['dy']:.6f}",f"{item['distance']:.6f}",
        ])

    parts=[
        "<section><h2>3 — Condição de módulo</h2>",
        _math(r"|G_c(s_d)G(s_d)H(s_d)|=1"),
        "<p>Chamando de <b>A</b> as distâncias aos polos e de <b>B</b> as distâncias aos zeros:</p>",
        _table(
            ["Distância","Elemento","Singularidade","ΔRe","ΔIm","Módulo"],
            rows,
        ),
    ]
    for item in m["poles"]+m["zeros"]:
        parts.append(_math(
            rf"{item['term']}=\sqrt{{({item['dx']:.6f})^2+({item['dy']:.6f})^2}}"
            rf"={item['distance']:.6f}"
        ))

    a_prod=r"\,\cdot\,".join(x["term"] for x in m["poles"]) or "1"
    b_prod=r"\,\cdot\,".join(x["term"] for x in m["zeros"]) or "1"
    parts.append(_math(
        rf"K_T=\frac{{\prod A_i}}{{\prod B_i}}"
        rf"=\frac{{{a_prod}}}{{{b_prod}}}"
        rf"=\frac{{{m['prod_a']:.8g}}}{{{m['prod_b']:.8g}}}"
        rf"={m['kt']:.8g}"
    ))
    parts.append(_math(
        rf"K_G={m['kg']:.8g},\qquad K_H={m['kh']:.8g},\qquad "
        rf"K_T=|K_GK_H|K_C"
    ))
    parts.append(_math(
        rf"K_C=\frac{{K_T}}{{|K_GK_H|}}"
        rf"=\frac{{{m['kt']:.8g}}}{{|({m['kg']:.8g})({m['kh']:.8g})|}}"
        rf"={m['kc']:.8g}"
    ))
    parts.append(f"<div class='ok'>Kc = {design.kc:.8g}</div>")
    parts.append("</section>")
    return "".join(parts)


def _controller_section(design):
    z=design.zero_parameter
    parts=["<section><h2>4 — Controlador</h2>"]
    if design.controller_type=="PD":
        parts.append(_math(rf"G_c(s)={design.kc:.8g}(s+{z:.8g})"))
        parts.append(_math(rf"K_p={design.kp:.8g},\qquad K_d={design.kd:.8g}"))
    elif design.controller_type=="PI":
        parts.append(_math(rf"G_c(s)={design.kc:.8g}\frac{{s+{z:.8g}}}{{s}}"))
        parts.append(_math(rf"K_p={design.kp:.8g},\qquad K_i={design.ki:.8g}"))
    else:
        parts.append(_math(rf"G_c(s)={design.kc:.8g}\frac{{(s+{z:.8g})^2}}{{s}}"))
        parts.append(_math(
            rf"K_p={design.kp:.8g},\qquad K_i={design.ki:.8g},\qquad K_d={design.kd:.8g}"
        ))
    parts.append("</section>")
    return "".join(parts)


def _verification_section(design):
    rows=[[_fmtc(p,6)] for p in sorted(design.closed_loop_poles,key=lambda x:(x.real,x.imag))]
    parts=[
        "<section><h2>6 — Verificação numérica</h2>",
        "<h3>Polos de malha fechada</h3>",
        _table(["Polo"],rows),
    ]
    m=design.metrics
    if m is not None:
        metric_rows=[
            ["Estável","Sim" if m.stable else "Não"],
            ["Mp","—" if m.overshoot_percent is None else f"{m.overshoot_percent:.4f}%"],
            [f"ts ({100*m.settling_band:.0f}%)","—" if m.settling_time is None else f"{m.settling_time:.6f} s"],
        ]
        parts.append(_table(["Métrica","Valor"],metric_rows))
    if len(design.refinement_history)>1:
        parts.append(
            f"<div class='note'>Ajuste fino aplicado em "
            f"{len(design.refinement_history)-1} iteração(ões), mantendo ξ.</div>"
        )
    parts.append("</section>")
    return "".join(parts)


def build_controller_html(design,data):
    focus=data.get("plot_focus")
    figures=[
        geometry_figure(design,focus=focus),
        root_locus_figure(
            data["num_g"],data["den_g"],data["num_h"],data["den_h"],
            desired_pole=design.desired_pole,
            title="LGR antes do controlador",
            focus=focus,
        ),
        root_locus_figure(
            data["num_g"],data["den_g"],data["num_h"],data["den_h"],
            desired_pole=design.desired_pole,
            controller_num=design.controller_num,
            controller_den=design.controller_den,
            title="LGR com o controlador projetado",
            padding=padding,
        ),
    ]
    step=step_response_figure(
        data["num_g"],data["den_g"],data["num_h"],data["den_h"],design
    )
    if step is not None:
        figures.append(step)

    plot_html=[]
    for i,fig in enumerate(figures):
        fig.update_layout(
            template="plotly_white",
            paper_bgcolor="white",
            plot_bgcolor="white",
            font=dict(color="#111111"),
        )
        plot_html.append(
            pio.to_html(
                fig,
                full_html=False,
                include_plotlyjs="inline" if i==0 else False,
                config={"responsive":True,"displaylogo":False},
            )
        )

    css="""
    :root { color-scheme: light; }
    html, body {
        background: #ffffff !important;
        color: #111111 !important;
    }
    body {
        font-family: Arial, Helvetica, sans-serif;
        max-width: 1500px; margin: 0 auto; padding: 24px;
        line-height: 1.45;
    }
    h1 { margin-bottom: 4px; }
    h2 { margin-top: 34px; border-bottom: 1px solid #aaa; padding-bottom: 6px; }
    h3 { margin-top: 20px; }
    .meta { opacity: .72; margin-bottom: 24px; }
    .math { text-align: center; overflow-x: auto; margin: 13px 0; }
    .table-wrap { overflow-x: auto; margin: 12px 0 18px; }
    table { width: 100%; border-collapse: collapse; font-size: 14px; }
    th, td { border: 1px solid #aaa; padding: 7px 9px; text-align: left; }
    th { background: #f2f2f2; color: #111111; }
    td { background: #ffffff; color: #111111; }
    .ok { background: rgba(44,160,44,.12); padding: 10px 12px; border-radius: 6px; margin: 10px 0; }
    .note { background: rgba(31,119,180,.12); padding: 10px 12px; border-radius: 6px; margin: 10px 0; }
    .charts { display: grid; grid-template-columns: 1fr; gap: 18px; }
    @media (min-width: 1100px) {
        .charts.two { grid-template-columns: 1fr 1fr; }
    }
    @media print {
        body { max-width: none; padding: 0; color: #000; background: #fff; }
        .plotly-graph-div { break-inside: avoid; }
        section { break-inside: auto; }
    }
    """

    chart_block=[
        "<section><h2>5 — Representação visual do projeto</h2>",
        "<div class='charts'>",plot_html[0],"</div>",
        "<div class='charts two'>",plot_html[1],plot_html[2],"</div>",
    ]
    if len(plot_html)>3:
        chart_block += ["<div class='charts'>",plot_html[3],"</div>"]
    chart_block.append("</section>")

    html_parts=[
        "<!doctype html><html lang='pt-BR'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        f"<title>Resolução completa - {escape(design.controller_type)}</title>",
        f"<style>{css}</style>",
        "<script>window.MathJax={tex:{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]},svg:{fontCache:'global'}};</script>",
        "<script async src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js'></script>",
        "</head><body>",
        f"<h1>Projeto {escape(design.controller_type)} pelo LGR</h1>",
        f"<div class='meta'>Relatório completo gerado pelo aplicativo LGR Controle U2. Enquadramento: {escape(focus_description(focus))}.</div>",
        _spec_section(design),
        _angle_section(design,data),
        _magnitude_section(design,data),
        _controller_section(design),
        "".join(chart_block),
        _verification_section(design),
        "</body></html>",
    ]
    return "".join(html_parts).encode("utf-8")
