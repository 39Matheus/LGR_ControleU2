import sympy as sp
import streamlit as st

from core.presets import EXERCISES
from core.specs import DesignSpecs
from .controller_view import render_controller_result,solve_controller
from .lgr_view import render_lgr_results


_DEFAULTS={
    "input_mode":"Projeto de controlador",
    "input_num_g":"1",
    "input_den_g":"1;3;2",
    "input_num_h":"1",
    "input_den_h":"1",
    "input_controller":"PD",
    "input_spec_kind":"Mp + ts",
    "input_mp":10.0,
    "input_ts":4.0,
    "input_band":5,
    "input_zeta":0.7,
    "input_wn":0.5,
    "input_pole_real":-4.0,
    "input_pole_imag":4.0,
    "input_test_real":-2.0,
    "input_test_imag":2.0,
}


def _hydrate_inputs():
    cache=st.session_state.setdefault("_input_cache",{})
    for key,default in _DEFAULTS.items():
        if key in st.session_state:
            cache[key]=st.session_state[key]
        else:
            st.session_state[key]=cache.get(key,default)


def _cache_inputs():
    cache=st.session_state.setdefault("_input_cache",{})
    for key in _DEFAULTS:
        if key in st.session_state:
            cache[key]=st.session_state[key]


def _coeff_text(values):
    return ";".join(f"{float(v):g}" for v in values)


def _parse(text):
    vals=[float(x.strip().replace(",",".")) for x in text.split(";") if x.strip()]
    if not vals:
        raise ValueError("Informe ao menos um coeficiente.")
    if abs(vals[0])<1e-15:
        raise ValueError("O primeiro coeficiente não pode ser zero. Remova zeros à esquerda.")
    return vals


def _poly_latex(coeffs):
    s=sp.symbols("s")
    degree=len(coeffs)-1
    expr=sum(sp.nsimplify(c,rational=True)*s**(degree-i) for i,c in enumerate(coeffs))
    return sp.latex(sp.expand(expr))


def _show_transfer_functions(data):
    c1,c2=st.columns(2)
    with c1:
        st.caption("Planta interpretada")
        st.latex(
            rf"G(s)=\frac{{{_poly_latex(data['num_g'])}}}"
            rf"{{{_poly_latex(data['den_g'])}}}"
        )
    with c2:
        st.caption("Realimentação interpretada")
        st.latex(
            rf"H(s)=\frac{{{_poly_latex(data['num_h'])}}}"
            rf"{{{_poly_latex(data['den_h'])}}}"
        )


def _load_preset(name):
    data=EXERCISES[name]
    specs=data["specs"]()

    values={
        "input_mode":"Projeto de controlador",
        "input_num_g":_coeff_text(data["num_g"]),
        "input_den_g":_coeff_text(data["den_g"]),
        "input_num_h":_coeff_text(data["num_h"]),
        "input_den_h":_coeff_text(data["den_h"]),
        "input_controller":data["controller"],
    }

    if specs.pole is not None:
        values.update({
            "input_spec_kind":"Polos desejados",
            "input_pole_real":float(specs.pole.real),
            "input_pole_imag":float(abs(specs.pole.imag)),
        })
    elif specs.zeta is not None and specs.omega_n is not None:
        values.update({
            "input_spec_kind":"ξ + ωn",
            "input_zeta":float(specs.zeta),
            "input_wn":float(specs.omega_n),
        })
    else:
        values.update({
            "input_spec_kind":"Mp + ts",
            "input_mp":float(specs.mp_max),
            "input_ts":float(specs.ts_max),
            "input_band":int(round(100*specs.settling_band)),
        })

    st.session_state.update(values)
    st.session_state.setdefault("_input_cache",{}).update(values)


def _common_data():
    return {
        "num_g":_parse(st.session_state["input_num_g"]),
        "den_g":_parse(st.session_state["input_den_g"]),
        "num_h":_parse(st.session_state["input_num_h"]),
        "den_h":_parse(st.session_state["input_den_h"]),
    }


def _controller_specs():
    kind=st.session_state["input_spec_kind"]
    if kind=="Mp + ts":
        return DesignSpecs.from_mp_ts(
            st.session_state["input_mp"],
            st.session_state["input_ts"],
            st.session_state["input_band"]/100.0,
        )
    if kind=="ξ + ωn":
        return DesignSpecs.from_zeta_wn(
            st.session_state["input_zeta"],
            st.session_state["input_wn"],
        )
    return DesignSpecs.from_pole(
        st.session_state["input_pole_real"],
        st.session_state["input_pole_imag"],
    )


def _controller_error(exc,data=None):
    message=str(exc)
    st.error(f"Falha no projeto: {message}")
    if data is not None:
        _show_transfer_functions(data)
    if "Fase incompatível" in message or "zero fora" in message:
        st.warning(
            "A geometria dessa função de transferência não permite o projeto na "
            "configuração escolhida para o polo especificado. Confira também zeros "
            "de coeficientes ausentes: s³+4s²+4s deve ser 1;4;4;0."
        )


def render_main_view():
    _hydrate_inputs()

    st.header("Entrada do problema")

    p1,p2=st.columns([3,1])
    with p1:
        preset=st.selectbox(
            "Preset opcional",
            ["— Nenhum —"]+list(EXERCISES.keys()),
            key="preset_choice",
            help="O preset apenas preenche os campos abaixo; depois você pode editar qualquer valor.",
        )
    with p2:
        st.write("")
        st.write("")
        if st.button("Preencher campos",use_container_width=True,disabled=preset=="— Nenhum —"):
            _load_preset(preset)
            st.rerun()

    st.radio(
        "Operação",
        ["Projeto de controlador","Análise do LGR (12 passos)"],
        horizontal=True,
        key="input_mode",
        help="A planta permanece nos mesmos campos ao trocar de operação.",
    )

    st.caption(
        "Coeficientes em ordem decrescente, separados por ponto e vírgula. "
        "Inclua 0 para potências ausentes; por exemplo, s³+4s²+4s → 1;4;4;0."
    )
    c1,c2=st.columns(2)
    c1.text_input("Numerador G(s)",key="input_num_g")
    c1.text_input("Denominador G(s)",key="input_den_g")
    c2.text_input("Numerador H(s)",key="input_num_h")
    c2.text_input("Denominador H(s)",key="input_den_h")

    parsed=None
    try:
        parsed=_common_data()
        with st.expander("Conferir função de transferência interpretada",expanded=False):
            _show_transfer_functions(parsed)
    except Exception as exc:
        st.warning(f"Entrada ainda inválida: {exc}")

    if st.session_state["input_mode"]=="Projeto de controlador":
        st.subheader("Parâmetros do projeto")
        st.selectbox("Controlador",["PD","PI","PID"],key="input_controller")
        st.radio(
            "Especificação",
            ["Mp + ts","ξ + ωn","Polos desejados"],
            horizontal=True,
            key="input_spec_kind",
        )

        if st.session_state["input_spec_kind"]=="Mp + ts":
            a,b,c=st.columns(3)
            a.number_input("Mp máx (%)",0.01,99.0,key="input_mp")
            b.number_input("ts máx (s)",0.001,key="input_ts")
            c.selectbox("Critério ts",[5,2],key="input_band")
        elif st.session_state["input_spec_kind"]=="ξ + ωn":
            a,b=st.columns(2)
            a.number_input("ξ",0.001,0.999,key="input_zeta")
            b.number_input("ωn (rad/s)",0.001,key="input_wn")
        else:
            a,b=st.columns(2)
            a.number_input("Re(s_d)",key="input_pole_real")
            b.number_input("Im(s_d)",0.001,key="input_pole_imag")

        if st.button("Projetar controlador",type="primary",use_container_width=True):
            try:
                data=_common_data()
                data["controller"]=st.session_state["input_controller"]
                data["specs"]=_controller_specs()
                design=solve_controller(data)
                st.session_state["last_controller_result"]=(design,data)
            except Exception as exc:
                _controller_error(exc,parsed)

        if "last_controller_result" in st.session_state:
            st.markdown("---")
            st.caption("Resultado da última execução do modo Projeto de controlador.")
            design,data=st.session_state["last_controller_result"]
            render_controller_result(design,data)

    else:
        st.subheader("Parâmetros da análise LGR")
        a,b=st.columns(2)
        a.number_input("Re(s_t)",key="input_test_real")
        b.number_input("Im(s_t)",key="input_test_imag")

        if st.button("Analisar LGR",type="primary",use_container_width=True):
            try:
                data=_common_data()
                point=complex(
                    st.session_state["input_test_real"],
                    st.session_state["input_test_imag"],
                )
                st.session_state["last_lgr_result"]=(data,point)
            except Exception as exc:
                st.error(f"Falha na entrada: {exc}")

        if "last_lgr_result" in st.session_state:
            st.markdown("---")
            st.caption("Resultado da última execução do modo Análise do LGR.")
            data,point=st.session_state["last_lgr_result"]
            try:
                render_lgr_results(data,point)
            except Exception as exc:
                st.error("Falha durante a análise do LGR.")
                st.exception(exc)

    _cache_inputs()
