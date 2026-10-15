import sympy as sp
import streamlit as st

from core.presets import EXERCISES,get_exercise
from core.specs import DesignSpecs
from .controller_view import render_controller_result,solve_controller
from .discretization_view import render_standalone_discretization
from .lgr_view import render_lgr_results


_DEFAULTS={
    "input_mode":"Projeto de controlador",
    "input_num_g":"1",
    "input_den_g":"1;3;2",
    "input_num_h":"1",
    "input_den_h":"1",

    "input_controller":"PD",
    "input_pid_zero_mode":"Iguais",
    "input_pid_fixed_zero":1.0,
    "input_generic_fixed_zeros":"",
    "input_generic_fixed_poles":"",
    "input_generic_unknown_kind":"Polo real",
    "input_generic_unknown_count":1,
    "input_generic_require_lhp":True,

    "input_spec_kind":"Mp + ts",
    "input_mp":10.0,
    "input_ts":4.0,
    "input_band":5,
    "input_auto_refine":False,
    "input_zeta":0.7,
    "input_wn":0.5,
    "input_pole_real":-4.0,
    "input_pole_imag":4.0,

    "input_discretize_after_design":False,
    "input_discretization_target":"Controlador Gc(s)",
    "input_discretization_method":"Tustin / Bilinear",
    "input_sample_time":1.0,
    "input_discrete_standalone_target":"G(s)",
    "input_discrete_method":"Tustin / Bilinear",
    "input_discrete_T":1.0,

    "input_plot_focus_mode":"Automático proporcional",
    "input_plot_x_percent":20.0,
    "input_plot_y_percent":10.0,
    "input_plot_x_padding":7.0,
    "input_plot_y_padding":7.0,

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


def _parse_optional(text):
    return [float(x.strip().replace(",",".")) for x in str(text).split(";") if x.strip()]


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
    data=get_exercise(name)
    specs=data["specs"]()
    cfg=data.get("controller_config",{})

    values={
        "input_mode":"Projeto de controlador",
        "input_num_g":_coeff_text(data["num_g"]),
        "input_den_g":_coeff_text(data["den_g"]),
        "input_num_h":_coeff_text(data["num_h"]),
        "input_den_h":_coeff_text(data["den_h"]),
        "input_controller":data["controller"],
        "input_auto_refine":False,
        "input_discretize_after_design":bool(data.get("discretize_after_design",False)),
        "input_discretization_target":data.get("discretization_target","Controlador Gc(s)"),
        "input_discretization_method":data.get("discretization_method","Tustin / Bilinear"),
        "input_sample_time":float(data.get("sample_time",1.0)),
    }

    if data["controller"]=="PID":
        mode=str(cfg.get("pid_zero_mode","equal")).lower()
        values["input_pid_zero_mode"]="Distintos" if mode=="distinct" else "Iguais"
        if "fixed_zero" in cfg:
            values["input_pid_fixed_zero"]=float(cfg["fixed_zero"])
    elif data["controller"]=="Genérico":
        values["input_generic_fixed_zeros"]=_coeff_text(cfg.get("fixed_zeros",[])) if cfg.get("fixed_zeros") else ""
        values["input_generic_fixed_poles"]=_coeff_text(cfg.get("fixed_poles",[])) if cfg.get("fixed_poles") else ""
        unknown=cfg.get("unknown_kind")
        values["input_generic_unknown_kind"]={
            "zero":"Zero real","pole":"Polo real",None:"Nenhuma — somente ganho"
        }.get(unknown,"Polo real")
        values["input_generic_unknown_count"]=int(cfg.get("unknown_count",1))
        values["input_generic_require_lhp"]=bool(cfg.get("require_lhp",False))

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


def _plot_focus_from_state():
    mode_map={
        "Automático proporcional":"auto_proportional",
        "Proporcional personalizado":"proportional",
        "Linear personalizado":"linear",
    }
    return {
        "mode":mode_map[st.session_state["input_plot_focus_mode"]],
        "x_percent":float(st.session_state["input_plot_x_percent"]),
        "y_percent":float(st.session_state["input_plot_y_percent"]),
        "x_padding":float(st.session_state["input_plot_x_padding"]),
        "y_padding":float(st.session_state["input_plot_y_padding"]),
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


def _controller_config_from_state():
    kind=st.session_state["input_controller"]
    if kind=="PID":
        if st.session_state["input_pid_zero_mode"]=="Iguais":
            return {"pid_zero_mode":"equal"}
        return {
            "pid_zero_mode":"distinct",
            "fixed_zero":float(st.session_state["input_pid_fixed_zero"]),
        }
    if kind=="Genérico":
        unknown={
            "Zero real":"zero",
            "Polo real":"pole",
            "Nenhuma — somente ganho":None,
        }[st.session_state["input_generic_unknown_kind"]]
        return {
            "fixed_zeros":_parse_optional(st.session_state["input_generic_fixed_zeros"]),
            "fixed_poles":_parse_optional(st.session_state["input_generic_fixed_poles"]),
            "unknown_kind":unknown,
            "unknown_count":(
                int(st.session_state["input_generic_unknown_count"])
                if unknown is not None else 0
            ),
            "require_lhp":bool(st.session_state["input_generic_require_lhp"]),
        }
    return {}


def _controller_error(exc,data=None):
    message=str(exc)
    st.error(f"Falha no projeto: {message}")
    if data is not None:
        _show_transfer_functions(data)
    if "fase" in message.lower() or "zero fora" in message.lower() or "polo do controlador" in message.lower():
        st.warning(
            "A geometria dessa função de transferência não permite o projeto na "
            "configuração escolhida para o polo especificado. Confira também zeros "
            "de coeficientes ausentes: s³+4s²+4s deve ser 1;4;4;0."
        )


def _render_plot_options():
    st.selectbox(
        "Enquadramento inicial dos gráficos",
        [
            "Automático proporcional",
            "Proporcional personalizado",
            "Linear personalizado",
        ],
        key="input_plot_focus_mode",
        help=(
            "O automático usa +20% do maior |Re| em x e +10% do maior |Im| em y. "
            "Você também pode definir percentuais ou margens lineares manualmente."
        ),
    )
    if st.session_state["input_plot_focus_mode"]=="Automático proporcional":
        st.caption("Padrão: +20% em x e +10% em y.")
    elif st.session_state["input_plot_focus_mode"]=="Proporcional personalizado":
        px,py=st.columns(2)
        px.number_input("Margem proporcional em x (%)",min_value=0.0,max_value=500.0,step=1.0,key="input_plot_x_percent")
        py.number_input("Margem proporcional em y (%)",min_value=0.0,max_value=500.0,step=1.0,key="input_plot_y_percent")
    else:
        px,py=st.columns(2)
        px.number_input("Margem linear em x",min_value=0.0,max_value=1000.0,step=1.0,key="input_plot_x_padding")
        py.number_input("Margem linear em y",min_value=0.0,max_value=1000.0,step=1.0,key="input_plot_y_padding")
    st.caption("O limite mínimo dos eixos continua sendo ±10.")


def _render_controller_structure():
    st.selectbox(
        "Controlador",
        ["P","PD","PI","PID","Genérico"],
        key="input_controller",
        help=(
            "O modo Genérico usa Gc(s)=Kc·Π(s+zi)/Π(s+pj), permitindo zeros/polos "
            "fixos e uma singularidade real a ser calculada pelo critério de ângulo."
        ),
    )
    kind=st.session_state["input_controller"]
    if kind=="PID":
        st.radio(
            "Zeros do PID",
            ["Iguais","Distintos"],
            horizontal=True,
            key="input_pid_zero_mode",
        )
        if st.session_state["input_pid_zero_mode"]=="Distintos":
            st.number_input(
                "Primeiro zero: z1 em s=-z1",
                value=float(st.session_state["input_pid_fixed_zero"]),
                key="input_pid_fixed_zero",
                help="z1 é informado; z2 é calculado pela condição de ângulo.",
            )
    elif kind=="Genérico":
        st.caption(
            r"Estrutura: \(G_c(s)=K_c\prod(s+z_i)/\prod(s+p_j)\). "
            "Os parâmetros z_i e p_j são as posições positivas das singularidades em -z_i e -p_j."
        )
        c1,c2=st.columns(2)
        c1.text_input(
            "Zeros fixos zi (separados por ;)",
            key="input_generic_fixed_zeros",
            placeholder="Ex.: 1;4",
        )
        c2.text_input(
            "Polos fixos pj (separados por ;)",
            key="input_generic_fixed_poles",
            placeholder="Ex.: 0;5",
        )
        c3,c4=st.columns(2)
        c3.selectbox(
            "Singularidade a calcular",
            ["Zero real","Polo real","Nenhuma — somente ganho"],
            key="input_generic_unknown_kind",
        )
        if st.session_state["input_generic_unknown_kind"]!="Nenhuma — somente ganho":
            c4.number_input(
                "Quantidade de singularidades incógnitas iguais",
                min_value=1,max_value=4,step=1,
                key="input_generic_unknown_count",
            )
        st.checkbox(
            "Exigir singularidade no semiplano esquerdo",
            key="input_generic_require_lhp",
        )


def _render_project_discretization_options():
    st.checkbox(
        "Discretizar após o projeto",
        key="input_discretize_after_design",
        help="Útil para as questões da Lista U2.2.",
    )
    if st.session_state["input_discretize_after_design"]:
        c1,c2,c3=st.columns(3)
        c1.selectbox(
            "Alvo da discretização",
            ["Controlador Gc(s)","G(s)","G(s)H(s)","G(s)H(s)Gc(s)"],
            key="input_discretization_target",
        )
        c2.selectbox(
            "Método",
            ["Euler / Forward","Backward","Tustin / Bilinear","Invariância ao degrau"],
            key="input_discretization_method",
        )
        c3.number_input(
            "Período T (s)",
            min_value=1e-9,
            key="input_sample_time",
            format="%.6g",
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

    selected=get_exercise(preset) if preset!="— Nenhum —" else None
    if selected:
        st.info(selected["description"])

    st.radio(
        "Operação",
        ["Projeto de controlador","Análise do LGR (12 passos)","Discretização"],
        horizontal=True,
        key="input_mode",
        help="G(s) e H(s) permanecem preenchidos ao trocar de operação.",
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

    mode=st.session_state["input_mode"]

    if mode=="Projeto de controlador":
        st.subheader("Parâmetros do projeto")
        _render_controller_structure()
        _render_plot_options()

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
            st.checkbox(
                "Aplicar ajuste fino após o projeto de fronteira",
                key="input_auto_refine",
                help=(
                    "Desmarcado: usa exatamente o polo de fronteira obtido de Mp e ts. "
                    "Marcado: simula a planta completa e desloca o polo mantendo ξ, se necessário."
                ),
            )
            if not st.session_state["input_auto_refine"]:
                st.caption("Modo prova: o controlador será calculado diretamente com o polo de fronteira.")
        elif st.session_state["input_spec_kind"]=="ξ + ωn":
            a,b=st.columns(2)
            a.number_input("ξ",0.001,0.999,key="input_zeta")
            b.number_input("ωn (rad/s)",0.001,key="input_wn")
        else:
            a,b=st.columns(2)
            a.number_input("Re(s_d)",key="input_pole_real")
            b.number_input("Im(s_d)",0.001,key="input_pole_imag")

        st.markdown("#### Etapa discreta opcional")
        _render_project_discretization_options()

        if st.button("Projetar controlador",type="primary",use_container_width=True):
            try:
                data=_common_data()
                data["controller"]=st.session_state["input_controller"]
                data["controller_config"]=_controller_config_from_state()
                data["specs"]=_controller_specs()
                data["auto_refine"]=(
                    st.session_state["input_auto_refine"]
                    if st.session_state["input_spec_kind"]=="Mp + ts"
                    else False
                )
                data["plot_focus"]=_plot_focus_from_state()
                data["discretize_after_design"]=bool(st.session_state["input_discretize_after_design"])
                data["discretization_target"]=st.session_state["input_discretization_target"]
                data["discretization_method"]=st.session_state["input_discretization_method"]
                data["sample_time"]=float(st.session_state["input_sample_time"])
                design=solve_controller(data)
                st.session_state["last_controller_result"]=(design,data)
            except Exception as exc:
                _controller_error(exc,parsed)

        if "last_controller_result" in st.session_state:
            st.markdown("---")
            st.caption("Resultado da última execução do modo Projeto de controlador.")
            design,data=st.session_state["last_controller_result"]
            data["plot_focus"]=_plot_focus_from_state()
            render_controller_result(design,data)

    elif mode=="Análise do LGR (12 passos)":
        st.subheader("Parâmetros da análise LGR")
        a,b=st.columns(2)
        a.number_input("Re(s_t)",key="input_test_real")
        b.number_input("Im(s_t)",key="input_test_imag")

        if st.button("Analisar LGR",type="primary",use_container_width=True):
            try:
                data=_common_data()
                point=complex(st.session_state["input_test_real"],st.session_state["input_test_imag"])
                st.session_state["last_lgr_result"]=(data,point)
            except Exception as exc:
                st.error(f"Falha na entrada: {exc}")

        if "last_lgr_result" in st.session_state:
            st.markdown("---")
            data,point=st.session_state["last_lgr_result"]
            try:
                render_lgr_results(data,point)
            except Exception as exc:
                st.error("Falha durante a análise do LGR.")
                st.exception(exc)

    else:
        if parsed is not None:
            render_standalone_discretization(parsed)

    _cache_inputs()
