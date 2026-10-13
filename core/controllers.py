"""Projeto de controladores pelo LGR com apresentação trigonométrica.

Além de PD/PI/PID, o núcleo aceita uma estrutura genérica formada por ganho,
zeros reais e polos reais. Uma singularidade real pode ser deixada como
incógnita e é determinada pela condição de ângulo; o ganho vem da condição de
módulo. Isso cobre, por exemplo, Gc(s)=a/(s+b) e PID com z1 != z2.
"""
import math

import numpy as np
import sympy as sp

from .geometry import phase_deg,singularity_contributions
from .models import ControllerDesign,TrigContribution
from .simulation import closed_loop_poles,meets_specs,step_metrics
from .specs import desired_pole,resolve_specs


def _as_float_list(values):
    if values is None:
        return []
    return [float(v) for v in values]


def _poly_from_locations(locations):
    """Polinômio monico com raízes em -location."""
    vals=_as_float_list(locations)
    if not vals:
        return np.asarray([1.0])
    return np.asarray(np.poly([-v for v in vals]),dtype=float)


def _fmt_factor(kind,value):
    symbol="z" if kind=="zero" else "p"
    return f"(s+{value:.8g})"


class ControllerDesigner:
    def __init__(self,num_g,den_g,num_h=(1.0,),den_h=(1.0,)):
        self.num_g=np.asarray(num_g,dtype=float)
        self.den_g=np.asarray(den_g,dtype=float)
        self.num_h=np.asarray(num_h,dtype=float)
        self.den_h=np.asarray(den_h,dtype=float)

    def loop_value(self,s):
        num=np.polyval(self.num_g,s)*np.polyval(self.num_h,s)
        den=np.polyval(self.den_g,s)*np.polyval(self.den_h,s)
        if abs(den)<1e-14:
            raise ZeroDivisionError("Ponto desejado coincide com um polo de malha aberta.")
        return complex(num/den)

    @staticmethod
    def _singularity_parameter(sd,angle_deg):
        """Para singularidade em -q: tan(theta)=wd/(q-sigma)."""
        sigma=-complex(sd).real
        wd=abs(complex(sd).imag)
        tangent=math.tan(math.radians(float(angle_deg)))
        if abs(tangent)<1e-12:
            raise ValueError("Ângulo da singularidade degenerado.")
        return sigma+wd/tangent

    @staticmethod
    def _controller_contributions(sd,zeros,poles,solved_kind=None,solved_count=0):
        out=[]
        zero_fixed_count=len(zeros)- (solved_count if solved_kind=="zero" else 0)
        pole_fixed_count=len(poles)- (solved_count if solved_kind=="pole" else 0)

        for i,z in enumerate(zeros,1):
            d=sd-complex(-z,0)
            role="solved" if solved_kind=="zero" and i>zero_fixed_count else "fixed"
            out.append(TrigContribution(
                f"z_c{i}",complex(-z,0),float(d.real),float(d.imag),
                phase_deg(d),"zero",role,
            ))
        for i,p in enumerate(poles,1):
            d=sd-complex(-p,0)
            role="solved" if solved_kind=="pole" and i>pole_fixed_count else "fixed"
            label="p_c=0" if abs(p)<1e-12 and i==1 else f"p_c{i}"
            out.append(TrigContribution(
                label,complex(-p,0),float(d.real),float(d.imag),
                phase_deg(d),"pole",role,
            ))
        return out

    def _design_template(self,controller_type,sd,config,settling_band=0.02):
        kind=str(controller_type)
        sd=complex(sd)
        if sd.real>=0 or sd.imag<=0:
            raise ValueError("Use o polo desejado do semiplano esquerdo com Im > 0.")

        fixed_zeros=_as_float_list(config.get("fixed_zeros",[]))
        fixed_poles=_as_float_list(config.get("fixed_poles",[]))
        unknown_kind=config.get("unknown_kind")
        unknown_count=int(config.get("unknown_count",1 if unknown_kind else 0))
        require_lhp=bool(config.get("require_lhp",False))

        if unknown_kind not in {None,"zero","pole"}:
            raise ValueError("A singularidade incógnita deve ser zero, polo ou nenhuma.")
        if unknown_kind is not None and unknown_count<1:
            raise ValueError("A quantidade de singularidades incógnitas deve ser positiva.")

        l0=self.loop_value(sd)
        fixed_factor=1+0j
        for z in fixed_zeros:
            fixed_factor*=sd+z
        for p in fixed_poles:
            if abs(sd+p)<1e-14:
                raise ValueError("O polo fixo do controlador coincide com o polo desejado.")
            fixed_factor/=sd+p

        base=l0*fixed_factor
        base_phase=phase_deg(base)

        solved_angle=0.0
        solved_parameter=None
        required_phase=0.0

        if unknown_kind=="zero":
            required_phase=(180.0-base_phase)%360.0
            solved_angle=required_phase/unknown_count
            if not 0.0<solved_angle<180.0:
                raise ValueError(
                    "A fase necessária não pode ser produzida pelos zeros reais "
                    "na configuração selecionada."
                )
            solved_parameter=self._singularity_parameter(sd,solved_angle)
            if require_lhp and solved_parameter<=0:
                raise ValueError(
                    "A solução exige zero fora da forma usual s=-z, z>0. "
                    "Revise a estrutura escolhida."
                )
            zeros=fixed_zeros+[solved_parameter]*unknown_count
            poles=fixed_poles
        elif unknown_kind=="pole":
            required_phase=(base_phase-180.0)%360.0
            solved_angle=required_phase/unknown_count
            if not 0.0<solved_angle<180.0:
                raise ValueError(
                    "A fase necessária não pode ser produzida pelos polos reais "
                    "na configuração selecionada."
                )
            solved_parameter=self._singularity_parameter(sd,solved_angle)
            if require_lhp and solved_parameter<=0:
                raise ValueError(
                    "A solução exige polo do controlador fora da forma s=-p, p>0."
                )
            zeros=fixed_zeros
            poles=fixed_poles+[solved_parameter]*unknown_count
        else:
            zeros=fixed_zeros
            poles=fixed_poles
            phase_mod=(base_phase%360.0+360.0)%360.0
            distance_to_lgr=min(abs(phase_mod-180.0),abs(phase_mod+180.0),abs(phase_mod-540.0))
            if distance_to_lgr>1e-5:
                raise ValueError(
                    f"Com todas as singularidades fixas, o ponto desejado não satisfaz "
                    f"a condição de ângulo (fase={base_phase:.5f}°)."
                )

        unit_num=_poly_from_locations(zeros)
        unit_den=_poly_from_locations(poles)
        unit_value=np.polyval(unit_num,sd)/np.polyval(unit_den,sd)
        kc=1.0/abs(unit_value*l0)

        num_c=np.asarray(unit_num,dtype=float)*kc
        den_c=np.asarray(unit_den,dtype=float)
        if abs(den_c[0])>1e-15:
            num_c=num_c/den_c[0]
            den_c=den_c/den_c[0]

        # Conversão para ganhos clássicos quando a estrutura for PID/PI/PD.
        kp=ki=kd=0.0
        upper=kind.upper()
        if upper=="PD" and len(zeros)==1 and not poles:
            kd=kc
            kp=kc*zeros[0]
        elif upper=="PI" and len(zeros)==1 and len(poles)==1 and abs(poles[0])<1e-12:
            kp=kc
            ki=kc*zeros[0]
        elif upper=="PID" and len(zeros)==2 and len(poles)==1 and abs(poles[0])<1e-12:
            kd=kc
            kp=kc*(zeros[0]+zeros[1])
            ki=kc*zeros[0]*zeros[1]
        elif upper=="P" and not zeros and not poles:
            kp=kc

        plant=singularity_contributions(
            np.polymul(self.num_g,self.num_h),
            np.polymul(self.den_g,self.den_h),
            sd,
        )
        ctrl_contrib=self._controller_contributions(
            sd,zeros,poles,unknown_kind,unknown_count
        )
        poles_cl=closed_loop_poles(
            self.num_g,self.den_g,self.num_h,self.den_h,num_c,den_c
        )
        metrics=step_metrics(
            self.num_g,self.den_g,self.num_h,self.den_h,
            num_c,den_c,settling_band,
        )

        s=sp.symbols("s")
        num_expr=sp.Poly.from_list(
            [sp.N(v,10) for v in num_c],gens=s
        ).as_expr()
        den_expr=sp.Poly.from_list(
            [sp.N(v,10) for v in den_c],gens=s
        ).as_expr()
        # Preserve a forma monica usada no projeto, em vez de deixar o SymPy
        # reescalar numerador/denominador por uma constante arbitrária.
        controller_latex=rf"\frac{{{sp.latex(sp.N(sp.expand(num_expr),8))}}}{{{sp.latex(sp.N(sp.expand(den_expr),8))}}}"

        if zeros:
            num_text="".join(_fmt_factor("zero",v) for v in zeros)
        else:
            num_text="1"
        if poles:
            den_text="".join(_fmt_factor("pole",v) for v in poles)
        else:
            den_text="1"
        controller_text=f"{kc:.8g}*{num_text}/{den_text}"

        zero_parameter=(
            float(solved_parameter)
            if unknown_kind=="zero" and solved_parameter is not None
            else (float(zeros[0]) if zeros else 0.0)
        )

        return ControllerDesign(
            controller_type=kind,
            desired_pole=sd,
            zero_parameter=zero_parameter,
            zero_locations=[-float(v) for v in zeros],
            kc=float(kc),
            kp=float(kp),ki=float(ki),kd=float(kd),
            controller_num=[float(v) for v in num_c],
            controller_den=[float(v) for v in den_c],
            base_phase_deg=float(base_phase),
            required_phase_deg=float(required_phase),
            zero_angle_deg=float(solved_angle if unknown_kind=="zero" else 0.0),
            plant_contributions=plant,
            controller_contributions=ctrl_contrib,
            closed_loop_poles=poles_cl,
            metrics=metrics,
            pole_locations=[-float(v) for v in poles],
            solved_singularity_kind=unknown_kind or "gain",
            solved_singularity_parameter=(
                float(solved_parameter) if solved_parameter is not None else None
            ),
            solved_singularity_angle_deg=(
                float(solved_angle) if unknown_kind is not None else None
            ),
            solved_singularity_count=unknown_count,
            controller_latex=controller_latex,
            controller_text=controller_text,
            controller_config=dict(config),
        )

    @staticmethod
    def _standard_config(controller_type,controller_config=None):
        cfg=dict(controller_config or {})
        kind=str(controller_type).upper()
        if kind=="P":
            return {"fixed_zeros":[],"fixed_poles":[],"unknown_kind":None,"unknown_count":0,"require_lhp":True}
        if kind=="PD":
            return {"fixed_zeros":[],"fixed_poles":[],"unknown_kind":"zero","unknown_count":1,"require_lhp":True}
        if kind=="PI":
            return {"fixed_zeros":[],"fixed_poles":[0.0],"unknown_kind":"zero","unknown_count":1,"require_lhp":True}
        if kind=="PID":
            mode=str(cfg.get("pid_zero_mode","equal")).lower()
            if mode in {"equal","iguais"}:
                return {"fixed_zeros":[],"fixed_poles":[0.0],"unknown_kind":"zero","unknown_count":2,"require_lhp":True}
            if mode in {"distinct","distintos"}:
                if "fixed_zero" not in cfg:
                    raise ValueError(
                        "Para PID com zeros distintos, informe um dos zeros. "
                        "O outro será calculado pela condição de ângulo."
                    )
                return {
                    "fixed_zeros":[float(cfg["fixed_zero"])],
                    "fixed_poles":[0.0],
                    "unknown_kind":"zero","unknown_count":1,"require_lhp":True,
                }
            raise ValueError("Modo de zeros do PID desconhecido.")
        if kind in {"GENÉRICO","GENERICO","GENERIC"}:
            return {
                "fixed_zeros":_as_float_list(cfg.get("fixed_zeros",[])),
                "fixed_poles":_as_float_list(cfg.get("fixed_poles",[])),
                "unknown_kind":cfg.get("unknown_kind"),
                "unknown_count":int(cfg.get("unknown_count",1 if cfg.get("unknown_kind") else 0)),
                "require_lhp":bool(cfg.get("require_lhp",False)),
            }
        raise ValueError("Use P, PD, PI, PID ou Genérico.")

    def design(self,controller_type,specs,auto_refine=False,controller_config=None):
        pole,zeta,wn,sigma=resolve_specs(specs)
        template=self._standard_config(controller_type,controller_config)

        if specs.pole is not None:
            source="pole"
        elif specs.zeta is not None and specs.omega_n is not None:
            source="zeta_wn"
        else:
            source="mp_ts"

        summary={
            "source":source,
            "initial_pole":complex(pole),
            "initial_zeta":float(zeta),
            "initial_omega_n":float(wn),
            "initial_sigma":float(sigma),
            "initial_omega_d":float(abs(complex(pole).imag)),
        }
        if specs.mp_max is not None:
            summary["mp_max"]=float(specs.mp_max)
        if specs.ts_max is not None:
            summary["ts_max"]=float(specs.ts_max)
            summary["settling_band"]=float(specs.settling_band)
        if specs.zeta is not None:
            summary["zeta_input"]=float(specs.zeta)
        if specs.omega_n is not None:
            summary["omega_n_input"]=float(specs.omega_n)
        if specs.pole is not None:
            summary["pole_input"]=complex(specs.pole)

        design=self._design_template(
            controller_type,pole,template,specs.settling_band
        )

        if auto_refine and specs.mp_max is not None and specs.ts_max is not None and design.metrics is not None:
            history=[]
            candidate=design
            factor=1.0
            for _ in range(60):
                m=candidate.metrics
                history.append({
                    "sigma":-candidate.desired_pole.real,
                    "kc":candidate.kc,
                    "singularity":(
                        candidate.solved_singularity_parameter
                        if candidate.solved_singularity_parameter is not None
                        else math.nan
                    ),
                    "mp":float(m.overshoot_percent) if m and m.overshoot_percent is not None else math.nan,
                    "ts":float(m.settling_time) if m and m.settling_time is not None else math.nan,
                })
                if meets_specs(m,specs.mp_max,specs.ts_max):
                    design=candidate
                    break
                factor*=1.02
                sigma_try=sigma*factor
                candidate=self._design_template(
                    controller_type,
                    desired_pole(zeta,sigma_try/zeta),
                    template,
                    specs.settling_band,
                )
            design.refinement_history=history

        final_pole=complex(design.desired_pole)
        final_wn=abs(final_pole)
        final_sigma=-final_pole.real
        final_zeta=final_sigma/final_wn
        summary.update({
            "zeta":float(final_zeta),
            "omega_n":float(final_wn),
            "sigma":float(final_sigma),
            "omega_d":float(abs(final_pole.imag)),
            "adjusted":bool(abs(final_pole-complex(pole))>1e-10),
        })
        design.specification_summary=summary
        design.controller_config=dict(controller_config or {})
        return design
