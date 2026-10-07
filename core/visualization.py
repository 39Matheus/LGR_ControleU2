"""Figuras de apoio ao projeto de controladores."""
import numpy as np
import plotly.graph_objects as go

try:
    import control as ctrl
except ImportError:
    ctrl=None

from .lgr import AnalisadorLGR
from .simulation import closed_loop_poles,closed_loop_transfer


def _finite_points(values):
    return [complex(v) for v in values if np.isfinite(complex(v).real) and np.isfinite(complex(v).imag)]


def _focus_ranges(points, padding=7.0):
    """Faixa inicial simétrica, ignorando ramos do LGR que tendem ao infinito."""
    pts=_finite_points(points)
    if not pts:
        return [-10.0,10.0],[-10.0,10.0]
    max_real=max(abs(p.real) for p in pts)
    max_imag=max(abs(p.imag) for p in pts)
    x_bound=max(10.0,max_real+padding)
    y_bound=max(10.0,max_imag+padding)
    return [-x_bound,x_bound],[-y_bound,y_bound]


def geometry_figure(design,padding=7.0):
    """Diagrama geométrico dos ângulos usados no ponto desejado."""
    fig=go.Figure()
    sd=complex(design.desired_pole)

    allc=design.plant_contributions+design.controller_contributions
    poles=[c for c in allc if c.kind=="pole"]
    zeros=[c for c in allc if c.kind=="zero"]

    if poles:
        fig.add_trace(go.Scatter(
            x=[c.singularity.real for c in poles],
            y=[c.singularity.imag for c in poles],
            mode="markers+text",
            marker=dict(symbol="x",size=12,line=dict(width=2)),
            text=[c.label for c in poles],
            textposition="top center",
            name="Polos",
        ))
    if zeros:
        fig.add_trace(go.Scatter(
            x=[c.singularity.real for c in zeros],
            y=[c.singularity.imag for c in zeros],
            mode="markers+text",
            marker=dict(symbol="circle-open",size=12,line=dict(width=2)),
            text=[c.label for c in zeros],
            textposition="top center",
            name="Zeros",
        ))

    for c in allc:
        fig.add_trace(go.Scatter(
            x=[c.singularity.real,sd.real],
            y=[c.singularity.imag,sd.imag],
            mode="lines",
            line=dict(width=1.5,dash="dot"),
            name=f"{c.label}: {c.angle_deg:.2f}°",
            hovertemplate=(
                f"{c.label}<br>ΔRe={c.dx:.5f}<br>ΔIm={c.dy:.5f}"
                f"<br>ângulo={c.angle_deg:.5f}°<extra></extra>"
            ),
        ))

    fig.add_trace(go.Scatter(
        x=[sd.real,sd.real],
        y=[sd.imag,-sd.imag],
        mode="markers",
        marker=dict(symbol="star",size=15,line=dict(width=1)),
        name="Polos desejados",
    ))
    fig.add_hline(y=0,line_width=1,line_color="black",opacity=.45)
    fig.add_vline(x=0,line_width=1,line_color="black",opacity=.45)
    focus_points=[c.singularity for c in allc]+[sd,sd.conjugate()]
    xrange,yrange=_focus_ranges(focus_points,padding=padding)
    fig.update_layout(
        title="Geometria da condição de ângulo",
        xaxis_title="Re(s)",yaxis_title="Im(s)",
        height=500,hovermode="closest",
        xaxis=dict(range=xrange),
        yaxis=dict(range=yrange,scaleanchor="x",scaleratio=1),
        margin=dict(l=20,r=20,t=50,b=20),
    )
    return fig


def root_locus_figure(num_g,den_g,num_h,den_h,desired_pole=None,controller_num=None,controller_den=None,title="LGR",padding=7.0):
    """LGR da planta original ou do sistema compensado."""
    ng=np.asarray(num_g,dtype=float)
    dg=np.asarray(den_g,dtype=float)
    if controller_num is not None:
        ng=np.polymul(np.asarray(controller_num,dtype=float),ng)
        dg=np.polymul(np.asarray(controller_den,dtype=float),dg)

    a=AnalisadorLGR(ng,dg,num_h,den_h)
    a.calcular_polos_zeros()
    a.calcular_segmentos_eixo_real()
    a.calcular_assintotas()
    a.adicionar_elementos_geometricos()
    a.calcular_lgr_exato()

    focus_points=list(a.polos)+list(a.zeros)
    if desired_pole is not None:
        sd=complex(desired_pole)
        focus_points.extend([sd,sd.conjugate()])
        a.fig.add_trace(go.Scatter(
            x=[sd.real,sd.real],
            y=[sd.imag,-sd.imag],
            mode="markers",
            marker=dict(symbol="star",size=14,line=dict(width=1)),
            name="Polo desejado",
        ))

    # O LGR possui ramos que podem tender ao infinito; usar todos os pontos no
    # autorange torna a região de interesse ilegível. O enquadramento inicial
    # usa apenas polos, zeros e polo desejado, com margem configurável (7 por padrão).
    # Zoom, pan e autoscale do Plotly continuam disponíveis manualmente.
    xrange,yrange=_focus_ranges(focus_points,padding=padding)
    a.fig.update_layout(
        title=title,
        height=500,
        xaxis=dict(range=xrange),
        yaxis=dict(range=yrange,scaleanchor="x",scaleratio=1),
    )
    return a.fig


def step_response_figure(num_g,den_g,num_h,den_h,design):
    """Compara a resposta ao degrau antes e depois da compensação quando possível."""
    if ctrl is None:
        return None

    original_poles=closed_loop_poles(num_g,den_g,num_h,den_h,[1.0],[1.0])
    controlled_poles=design.closed_loop_poles

    stable_rates=[-p.real for p in controlled_poles if p.real< -1e-8]
    if stable_rates:
        slow=min(stable_rates)
        t_end=min(max(10.0/slow,8.0),60.0)
    else:
        t_end=10.0
    t=np.linspace(0,t_end,2500)

    fig=go.Figure()

    def add_curve(num_c,den_c,label):
        num,den=closed_loop_transfer(num_g,den_g,num_h,den_h,num_c,den_c)
        sys=ctrl.TransferFunction(num,den)
        tout,yout=ctrl.step_response(sys,T=t)
        y=np.real(np.asarray(yout).squeeze())
        finite=np.isfinite(y)
        if not finite.any():
            return
        # Evita destruir a escala do gráfico quando a resposta original é instável.
        if np.nanmax(np.abs(y[finite]))>50:
            return
        fig.add_trace(go.Scatter(x=tout,y=y,mode="lines",name=label))

    if original_poles and all(p.real< -1e-8 for p in original_poles):
        add_curve([1.0],[1.0],"Sem controlador")
    add_curve(design.controller_num,design.controller_den,"Com controlador")

    fig.add_hline(y=1.0,line_width=1,line_dash="dash",opacity=.55,annotation_text="referência = 1")
    fig.update_layout(
        title="Resposta ao degrau em malha fechada",
        xaxis_title="Tempo (s)",yaxis_title="Saída",
        height=450,hovermode="x unified",
        margin=dict(l=20,r=20,t=50,b=20),
    )
    return fig
