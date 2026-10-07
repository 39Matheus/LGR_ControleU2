import math
import streamlit as st

from core.visualization import geometry_figure,root_locus_figure,step_response_figure

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
    st.plotly_chart(geometry_figure(d),use_container_width=True,key="geometry_controller")

    st.subheader("2 — Condição de ângulo")
    st.caption("Apresentação trigonométrica; atan2 é usado internamente apenas para o quadrante correto.")
    render_trig_table(d.plant_contributions)
    st.latex(rf"\angle L_{{base}}={d.base_phase_deg:.6f}^\circ")
    st.latex(rf"\phi_{{nec}}={d.required_phase_deg:.6f}^\circ")
    if d.controller_type=="PID":
        st.latex(rf"\phi_1=\phi_2={d.zero_angle_deg:.6f}^\circ")
    else:
        st.latex(rf"\phi={d.zero_angle_deg:.6f}^\circ")
    st.latex(rf"\tan\phi=\frac{{\omega_d}}{{z-\sigma}}\Rightarrow z={d.zero_parameter:.6f}")
    st.success("Zero(s): "+", ".join(f"s={z:.6f}" for z in d.zero_locations))

    st.subheader("3 — Condição de módulo")
    st.latex(r"|G_c(s_d)G(s_d)H(s_d)|=1")
    st.latex(rf"K_c={d.kc:.8g}")

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
        "Os LGRs abrem inicialmente em uma região de interesse com margem de 10 unidades "
        "além dos polos, zeros e polo desejado. O zoom e o autoscale do Plotly continuam disponíveis."
    )
    left,right=st.columns(2)
    with left:
        st.plotly_chart(
            root_locus_figure(
                data["num_g"],data["den_g"],data["num_h"],data["den_h"],
                desired_pole=d.desired_pole,
                title="LGR antes do controlador",
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
