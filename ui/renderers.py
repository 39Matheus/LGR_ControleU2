import math
import streamlit as st

from core.geometry import angle_breakdown,magnitude_breakdown
from core.visualization import focus_description,geometry_figure,root_locus_figure,step_response_figure

def fmt_complex(z,digits=5):
    z=complex(z)
    if abs(z.imag)<10**(-digits): return f"{z.real:.{digits}f}"
    return f"{z.real:.{digits}f} {'+' if z.imag>=0 else '-'} {abs(z.imag):.{digits}f}j"

def render_trig_table(contribs):
    rows=[]
    for c in contribs:
        base=90.0 if abs(c.dx)<1e-12 else math.degrees(math.atan(abs(c.dy/c.dx)))
        rows.append({"Elemento":c.label,"Singularidade":fmt_complex(c.singularity),
                     "ΔRe":f"{c.dx:.5f}","ΔIm":f"{c.dy:.5f}",
                     "atan(|ΔIm/ΔRe|)":f"{base:.5f}°","ângulo":f"{c.angle_deg:.5f}°"})
    st.dataframe(rows,hide_index=True,use_container_width=True)


def _angle_symbol(c,plant_ids,controller_zero_count):
    if c.kind=="zero":
        if c.label.startswith("z_c"):
            if controller_zero_count==1:
                return r"\phi_c"
            idx=int(c.label.replace("z_c",""))
            return rf"\phi_{{c{idx}}}"
        idx=plant_ids["zero"].get(c.label,1)
        return rf"\phi_{{{idx}}}"
    if c.label=="p_c=0":
        return r"\theta_c"
    idx=plant_ids["pole"].get(c.label,1)
    return rf"\theta_{{{idx}}}"


def _render_angle_condition(d,data):
    info=angle_breakdown(d,data["num_g"],data["den_g"],data["num_h"],data["den_h"])
    all_contrib=(
        info["plant_zeros"]
        +info["plant_poles"]
        +info["controller_poles"]
        +info["controller_zeros"]
    )
    ncz=len(info["controller_zeros"])

    rows=[]
    for c in all_contrib:
        symbol=_angle_symbol(c,info["plant_ids"],ncz)
        base=90.0 if abs(c.dx)<1e-12 else math.degrees(math.atan(abs(c.dy/c.dx)))
        rows.append({
            "Termo":symbol.replace("\\",""),
            "Elemento":c.label,
            "Singularidade":fmt_complex(c.singularity),
            "ΔRe":f"{c.dx:.5f}",
            "ΔIm":f"{c.dy:.5f}",
            "atan(|ΔIm/ΔRe|)":f"{base:.5f}°",
            "ângulo usado":f"{c.angle_deg:.5f}°",
        })
    st.dataframe(rows,hide_index=True,use_container_width=True)

    st.markdown("**Montagem do critério do ângulo**")
    st.latex(
        r"\sum\phi_{\mathrm{zeros}}-\sum\theta_{\mathrm{polos}}"
        r"=(2q+1)180^\circ"
    )

    zero_terms=[f"{c.angle_deg:.5f}" for c in info["plant_zeros"]]
    pole_terms=[f"{c.angle_deg:.5f}" for c in info["plant_poles"]+info["controller_poles"]]
    zero_sum=" + ".join(zero_terms) if zero_terms else "0"
    pole_sum=" + ".join(pole_terms) if pole_terms else "0"
    gain_phase=info["constant_phase"]

    if gain_phase:
        st.latex(
            rf"{gain_phase:.5f}+({zero_sum})-({pole_sum})"
            rf"={info['raw_base']:.5f}^\circ"
        )
        st.caption("O ganho constante negativo acrescenta 180° à fase.")
    else:
        st.latex(
            rf"({zero_sum})-({pole_sum})"
            rf"={info['raw_base']:.5f}^\circ"
        )

    if ncz==1:
        st.latex(
            rf"{info['raw_base']:.5f}^\circ+\phi_c"
            rf"={info['target']:.0f}^\circ"
        )
        st.latex(
            rf"\phi_c={info['target']:.0f}^\circ-({info['raw_base']:.5f}^\circ)"
            rf"={d.zero_angle_deg:.5f}^\circ"
        )
    else:
        st.latex(
            rf"{info['raw_base']:.5f}^\circ+{ncz}\phi_c"
            rf"={info['target']:.0f}^\circ"
        )
        st.latex(
            rf"\phi_c=\frac{{{info['target']:.0f}^\circ-({info['raw_base']:.5f}^\circ)}}"
            rf"{{{ncz}}}={d.zero_angle_deg:.5f}^\circ"
        )

    st.caption(
        f"Conferência: a soma final vale {info['final_raw']:.5f}°, "
        "equivalente a 180° módulo 360°."
    )

    sigma=-d.desired_pole.real
    wd=abs(d.desired_pole.imag)
    st.markdown("**Posição do zero do controlador pela geometria**")
    st.latex(
        rf"\tan(\phi_c)=\frac{{\omega_d}}{{z_c-\sigma}}"
        rf"=\frac{{{wd:.6f}}}{{z_c-{sigma:.6f}}}"
    )
    st.latex(
        rf"z_c=\sigma+\frac{{\omega_d}}{{\tan(\phi_c)}}"
        rf"={sigma:.6f}+\frac{{{wd:.6f}}}{{\tan({d.zero_angle_deg:.5f}^\circ)}}"
        rf"={d.zero_parameter:.6f}"
    )
    st.success(
        "Zero(s) do controlador: "
        +", ".join(f"s={z:.6f}" for z in d.zero_locations)
    )


def _render_magnitude_condition(d,data):
    m=magnitude_breakdown(
        d,data["num_g"],data["den_g"],data["num_h"],data["den_h"]
    )

    st.latex(r"|G_c(s_d)G(s_d)H(s_d)|=1")
    st.markdown(
        "Chamando de **A** as distâncias do polo desejado aos polos e de **B** "
        "as distâncias aos zeros:"
    )

    rows=[]
    for item in m["poles"]+m["zeros"]:
        rows.append({
            "Distância":item["term"],
            "Elemento":item["source_label"],
            "Singularidade":fmt_complex(item["singularity"]),
            "ΔRe":f"{item['dx']:.6f}",
            "ΔIm":f"{item['dy']:.6f}",
            "Módulo":f"{item['distance']:.6f}",
        })
    st.dataframe(rows,hide_index=True,use_container_width=True)

    for item in m["poles"]+m["zeros"]:
        term=item["term"]
        st.latex(
            rf"{term}=|s_d-s_{{{item['source_label']}}}|"
            rf"=\sqrt{{({item['dx']:.6f})^2+({item['dy']:.6f})^2}}"
            rf"={item['distance']:.6f}"
        )

    a_prod=r"\,\cdot\,".join(item["term"] for item in m["poles"]) or "1"
    b_prod=r"\,\cdot\,".join(item["term"] for item in m["zeros"]) or "1"
    st.markdown("**Ganho total exigido pelo critério de módulo**")
    st.latex(
        rf"K_T=\frac{{\prod A_i}}{{\prod B_i}}"
        rf"=\frac{{{a_prod}}}{{{b_prod}}}"
        rf"=\frac{{{m['prod_a']:.8g}}}{{{m['prod_b']:.8g}}}"
        rf"={m['kt']:.8g}"
    )

    st.markdown("**Separação dos ganhos da planta, realimentação e controlador**")
    st.latex(
        rf"K_G=\frac{{a_{{0,G}}}}{{b_{{0,G}}}}"
        rf"={m['kg']:.8g},\qquad "
        rf"K_H=\frac{{a_{{0,H}}}}{{b_{{0,H}}}}"
        rf"={m['kh']:.8g}"
    )
    st.latex(r"K_T=|K_GK_H|K_C")
    st.latex(
        rf"K_C=\frac{{K_T}}{{|K_GK_H|}}"
        rf"=\frac{{{m['kt']:.8g}}}{{|({m['kg']:.8g})({m['kh']:.8g})|}}"
        rf"={m['kc']:.8g}"
    )
    st.success(f"Kc = {d.kc:.8g}")

def _render_desired_pole_derivation(d):
    s=d.specification_summary
    source=s.get("source")
    initial_pole=complex(s.get("initial_pole",d.desired_pole))
    zeta0=float(s.get("initial_zeta",s.get("zeta",float("nan"))))
    wn0=float(s.get("initial_omega_n",s.get("omega_n",float("nan"))))
    sigma0=float(s.get("initial_sigma",-initial_pole.real))
    wd0=float(s.get("initial_omega_d",abs(initial_pole.imag)))

    st.subheader("1 — Obtenção do polo desejado")

    if source=="mp_ts":
        mp=float(s["mp_max"])
        ts=float(s["ts_max"])
        band=float(s["settling_band"])
        c=3.0 if abs(band-0.05)<1e-12 else 4.0
        st.markdown("**a) Amortecimento a partir do máximo sobresinal**")
        st.latex(
            rf"M_p={mp:.6g}\%={mp/100.0:.6g}"
        )
        st.latex(
            r"\xi="
            r"\frac{-\ln(M_p)}{\sqrt{\pi^2+\ln^2(M_p)}}"
            rf"={zeta0:.6f}"
        )

        st.markdown("**b) Parte real a partir do tempo de acomodação**")
        st.latex(
            rf"t_s({100*band:.0f}\%)\approx\frac{{{c:g}}}{{\xi\omega_n}}"
            rf"=\frac{{{c:g}}}{{\sigma}}"
        )
        st.latex(
            rf"\sigma=\frac{{{c:g}}}{{t_s}}"
            rf"=\frac{{{c:g}}}{{{ts:.6g}}}={sigma0:.6f}"
        )

        st.markdown("**c) Frequência natural, frequência amortecida e polos**")
        st.latex(
            rf"\omega_n=\frac{{\sigma}}{{\xi}}"
            rf"=\frac{{{sigma0:.6f}}}{{{zeta0:.6f}}}={wn0:.6f}"
        )
        st.latex(
            rf"\omega_d=\omega_n\sqrt{{1-\xi^2}}={wd0:.6f}"
        )
        st.latex(
            rf"s_{{d,0}}=-\sigma\pm j\omega_d"
            rf"={initial_pole.real:.6f}\pm j{abs(initial_pole.imag):.6f}"
        )

        if s.get("adjusted"):
            final_pole=complex(d.desired_pole)
            st.info(
                "Esse é o polo de fronteira obtido diretamente das especificações. "
                "Como Mp e ts são limites e a planta completa não é exatamente de 2ª ordem, "
                "o aplicativo verificou a resposta e deslocou o polo mantendo ξ até satisfazer "
                "as especificações."
            )
            st.latex(
                rf"s_d={final_pole.real:.6f}\pm j{abs(final_pole.imag):.6f}"
            )
            st.latex(
                rf"\sigma_f={s['sigma']:.6f},\quad "
                rf"\omega_{{n,f}}={s['omega_n']:.6f},\quad "
                rf"\xi_f={s['zeta']:.6f}"
            )
        else:
            st.success("O polo de fronteira já foi aceito como polo de projeto.")

    elif source=="zeta_wn":
        st.markdown(r"Como \(\xi\) e \(\omega_n\) são dados, usa-se diretamente o modelo dominante de 2ª ordem.")
        st.latex(
            rf"\sigma=\xi\omega_n={zeta0:.6f}\cdot{wn0:.6f}={sigma0:.6f}"
        )
        st.latex(
            rf"\omega_d=\omega_n\sqrt{{1-\xi^2}}"
            rf"={wn0:.6f}\sqrt{{1-{zeta0:.6f}^2}}={wd0:.6f}"
        )
        st.latex(
            rf"s_d=-\xi\omega_n\pm j\omega_n\sqrt{{1-\xi^2}}"
            rf"={initial_pole.real:.6f}\pm j{abs(initial_pole.imag):.6f}"
        )

    else:
        st.markdown("Os polos desejados foram fornecidos diretamente no enunciado.")
        st.latex(
            rf"s_d={initial_pole.real:.6f}\pm j{abs(initial_pole.imag):.6f}"
        )
        st.latex(
            rf"\omega_n=|s_d|={wn0:.6f},\qquad "
            rf"\xi=\frac{{-\operatorname{{Re}}(s_d)}}{{\omega_n}}={zeta0:.6f}"
        )


def render_controller_design(d,data):
    _render_desired_pole_derivation(d)
    focus=data.get("plot_focus")
    st.plotly_chart(geometry_figure(d,focus=focus),use_container_width=True,key="geometry_controller")

    st.subheader("2 — Condição de ângulo")
    st.caption(
        "Convenção usada em aula: soma dos ângulos dos zeros menos soma dos "
        "ângulos dos polos. O atan2 é usado internamente apenas para determinar o quadrante."
    )
    _render_angle_condition(d,data)

    st.subheader("3 — Condição de módulo")
    _render_magnitude_condition(d,data)

    st.subheader("4 — Controlador")
    z=d.zero_parameter
    if d.controller_type=="PD":
        st.latex(rf"G_c(s)={d.kc:.8g}(s+{z:.8g})")
        st.latex(rf"K_p={d.kp:.8g},\quad K_d={d.kd:.8g}")
    elif d.controller_type=="PI":
        st.latex(rf"G_c(s)={d.kc:.8g}\frac{{s+{z:.8g}}}{{s}}")
        st.latex(rf"K_p={d.kp:.8g},\quad K_i={d.ki:.8g}")
    else:
        st.latex(rf"G_c(s)={d.kc:.8g}\frac{{(s+{z:.8g})^2}}{{s}}")
        st.latex(rf"K_p={d.kp:.8g},\quad K_i={d.ki:.8g},\quad K_d={d.kd:.8g}")

    st.subheader("5 — Representação visual do projeto")
    st.caption(
        "Enquadramento inicial: "
        +focus_description(focus)
        +". O zoom e o autoscale do Plotly continuam disponíveis."
    )
    left,right=st.columns(2)
    with left:
        st.plotly_chart(
            root_locus_figure(
                data["num_g"],data["den_g"],data["num_h"],data["den_h"],
                desired_pole=d.desired_pole,
                title="LGR antes do controlador",
                focus=focus,
            ),
            use_container_width=True,
            key="lgr_before_controller",
        )
    with right:
        st.plotly_chart(
            root_locus_figure(
                data["num_g"],data["den_g"],data["num_h"],data["den_h"],
                desired_pole=d.desired_pole,
                controller_num=d.controller_num,
                controller_den=d.controller_den,
                title="LGR com o controlador projetado",
                focus=focus,
            ),
            use_container_width=True,
            key="lgr_after_controller",
        )

    step_fig=step_response_figure(data["num_g"],data["den_g"],data["num_h"],data["den_h"],d)
    if step_fig is not None:
        st.plotly_chart(step_fig,use_container_width=True,key="step_controller")

    st.subheader("6 — Verificação numérica")
    st.write("Polos de malha fechada:")
    st.code("\n".join(fmt_complex(p,6) for p in sorted(d.closed_loop_poles,key=lambda x:(x.real,x.imag))))
    if d.metrics:
        c1,c2,c3=st.columns(3)
        c1.metric("Estável","Sim" if d.metrics.stable else "Não")
        c2.metric("Mp","—" if d.metrics.overshoot_percent is None else f"{d.metrics.overshoot_percent:.3f}%")
        c3.metric(f"ts ({100*d.metrics.settling_band:.0f}%)","—" if d.metrics.settling_time is None else f"{d.metrics.settling_time:.4f}s")
    if len(d.refinement_history)>1:
        st.info(f"Ajuste fino automático aplicado em {len(d.refinement_history)-1} iteração(ões), aumentando σ e mantendo ξ.")
