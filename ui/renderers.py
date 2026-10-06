import math
import streamlit as st

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

def render_controller_design(d):
    s=d.specification_summary
    st.subheader("1 — Especificações e polo desejado")
    st.latex(rf"s_d={d.desired_pole.real:.6f}+j{d.desired_pole.imag:.6f}")
    st.latex(rf"\xi={s.get('zeta',float('nan')):.6f},\quad \omega_n={s.get('omega_n',float('nan')):.6f}")
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
    st.subheader("5 — Verificação")
    st.code("\n".join(fmt_complex(p,6) for p in sorted(d.closed_loop_poles,key=lambda x:(x.real,x.imag))))
    if d.metrics:
        c1,c2,c3=st.columns(3)
        c1.metric("Estável","Sim" if d.metrics.stable else "Não")
        c2.metric("Mp","—" if d.metrics.overshoot_percent is None else f"{d.metrics.overshoot_percent:.3f}%")
        c3.metric(f"ts ({100*d.metrics.settling_band:.0f}%)","—" if d.metrics.settling_time is None else f"{d.metrics.settling_time:.4f}s")
    if len(d.refinement_history)>1:
        st.info(f"Ajuste fino automático aplicado em {len(d.refinement_history)-1} iteração(ões), aumentando σ e mantendo ξ.")
