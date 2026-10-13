import streamlit as st

from core.discretization import discretize_transfer


_METHOD_OPTIONS=[
    "Euler / Forward",
    "Backward",
    "Tustin / Bilinear",
    "Invariância ao degrau",
]


def render_discretization_result(result,heading="7 — Discretização"):
    st.subheader(heading)
    st.write(f"**Alvo:** {result.target_label}")
    c1,c2=st.columns(2)
    c1.metric("Método",result.method)
    c2.metric("Período de amostragem T",f"{result.sample_time:g} s")

    st.markdown("**Função contínua**")
    st.latex(rf"{result.target_label}={result.continuous_latex}")

    st.markdown("**Relação de aproximação**")
    st.latex(result.substitution_latex)

    st.markdown("**Substituição e simplificação**")
    if result.substituted_latex.startswith(r"\text"):
        st.latex(result.substituted_latex)
    else:
        st.latex(rf"H(z)={result.substituted_latex}")

    st.markdown("**Função discretizada**")
    st.latex(rf"H(z)={result.discrete_latex}")

    rows=[]
    n=max(len(result.discrete_num),len(result.discrete_den))
    for i in range(n):
        rows.append({
            "índice":i,
            "numerador":(
                f"{result.discrete_num[i]:.10g}"
                if i<len(result.discrete_num) else ""
            ),
            "denominador":(
                f"{result.discrete_den[i]:.10g}"
                if i<len(result.discrete_den) else ""
            ),
        })
    st.dataframe(rows,hide_index=True,use_container_width=True)

    if result.difference_equation_latex:
        st.markdown("**Equação de diferenças correspondente**")
        st.latex(result.difference_equation_latex)

    if result.poles or result.zeros:
        c1,c2=st.columns(2)
        with c1:
            st.markdown("**Polos em z**")
            if result.poles:
                st.code("\n".join(
                    f"{p.real:.6f} {'+' if p.imag>=0 else '-'} {abs(p.imag):.6f}j"
                    for p in result.poles
                ))
            else:
                st.write("Nenhum.")
        with c2:
            st.markdown("**Zeros em z**")
            if result.zeros:
                st.code("\n".join(
                    f"{z.real:.6f} {'+' if z.imag>=0 else '-'} {abs(z.imag):.6f}j"
                    for z in result.zeros
                ))
            else:
                st.write("Nenhum.")

    for note in result.notes:
        st.info(note)


def render_standalone_discretization(data):
    st.subheader("Discretização de função de transferência")
    st.caption(
        "Métodos do capítulo 5: Forward/Euler, Backward, Tustin e invariância ao degrau."
    )

    target=st.selectbox(
        "Função a discretizar",
        ["G(s)","H(s)","G(s)H(s)"],
        key="input_discrete_standalone_target",
    )
    method=st.selectbox(
        "Método",
        _METHOD_OPTIONS,
        key="input_discrete_method",
    )
    T=st.number_input(
        "Período de amostragem T (s)",
        min_value=1e-9,
        value=float(st.session_state.get("input_discrete_T",1.0)),
        key="input_discrete_T",
        format="%.6g",
    )

    if st.button("Discretizar",type="primary",use_container_width=True):
        from core.discretization import multiply_transfers
        try:
            if target=="G(s)":
                num,den=data["num_g"],data["den_g"]
            elif target=="H(s)":
                num,den=data["num_h"],data["den_h"]
            else:
                num,den=multiply_transfers(
                    (data["num_g"],data["den_g"]),
                    (data["num_h"],data["den_h"]),
                )
            result=discretize_transfer(num,den,T,method,target_label=target)
            st.session_state["last_discretization_result"]=result
        except Exception as exc:
            st.error(f"Falha na discretização: {exc}")

    if "last_discretization_result" in st.session_state:
        st.markdown("---")
        render_discretization_result(
            st.session_state["last_discretization_result"],
            heading="Resultado da discretização",
        )
