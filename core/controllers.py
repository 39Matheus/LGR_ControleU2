"""Projeto de PD, PI e PID pelo LGR com apresentação trigonométrica."""
import math
import numpy as np
from .geometry import phase_deficiency,phase_deg,real_zero_parameter,singularity_contributions,split_zero_phase
from .models import ControllerDesign,TrigContribution
from .simulation import closed_loop_poles,meets_specs,step_metrics
from .specs import desired_pole,resolve_specs

class ControllerDesigner:
    def __init__(self,num_g,den_g,num_h=(1.0,),den_h=(1.0,)):
        self.num_g=np.asarray(num_g,dtype=float); self.den_g=np.asarray(den_g,dtype=float)
        self.num_h=np.asarray(num_h,dtype=float); self.den_h=np.asarray(den_h,dtype=float)

    def loop_value(self,s):
        num=np.polyval(self.num_g,s)*np.polyval(self.num_h,s)
        den=np.polyval(self.den_g,s)*np.polyval(self.den_h,s)
        if abs(den)<1e-14: raise ZeroDivisionError("Polo de malha aberta no ponto desejado.")
        return complex(num/den)

    def _design_at_pole(self,controller_type,sd,settling_band=0.02):
        kind=controller_type.upper()
        if kind not in {"PD","PI","PID"}: raise ValueError("Use PD, PI ou PID.")
        sd=complex(sd)
        if sd.real>=0 or sd.imag<=0: raise ValueError("Polo desejado inválido.")
        l0=self.loop_value(sd)
        origin_pole=kind in {"PI","PID"}
        nzeros=2 if kind=="PID" else 1
        base=l0/sd if origin_pole else l0
        base_phase=phase_deg(base)
        required=phase_deficiency(base_phase)
        try:
            phi=split_zero_phase(required,nzeros)
        except ValueError as exc:
            raise ValueError(
                f"Fase incompatível com zeros reais nesta configuração: "
                f"fase base={base_phase:.3f}°, contribuição necessária={required:.3f}°, "
                f"{nzeros} zero(s)."
            ) from exc
        z=real_zero_parameter(sd,phi)
        if z<=0: raise ValueError("A solução exige zero fora da forma s=-z com z>0.")
        factor=(sd+z)**nzeros
        if origin_pole: factor/=sd
        kc=1.0/abs(factor*l0)
        if kind=="PD":
            num_c=[kc,kc*z]; den_c=[1.0]; kp,ki,kd=kc*z,0.0,kc
        elif kind=="PI":
            num_c=[kc,kc*z]; den_c=[1.0,0.0]; kp,ki,kd=kc,kc*z,0.0
        else:
            num_c=[kc,2*kc*z,kc*z*z]; den_c=[1.0,0.0]; kp,ki,kd=2*kc*z,kc*z*z,kc
        plant=singularity_contributions(np.polymul(self.num_g,self.num_h),np.polymul(self.den_g,self.den_h),sd)
        ctrl_contrib=[]
        if origin_pole:
            d=sd
            ctrl_contrib.append(TrigContribution("p_c=0",0j,float(d.real),float(d.imag),phase_deg(d),"pole"))
        for i in range(nzeros):
            loc=complex(-z,0); d=sd-loc
            ctrl_contrib.append(TrigContribution(f"z_c{i+1}",loc,float(d.real),float(d.imag),phase_deg(d),"zero"))
        poles=closed_loop_poles(self.num_g,self.den_g,self.num_h,self.den_h,num_c,den_c)
        metrics=step_metrics(self.num_g,self.den_g,self.num_h,self.den_h,num_c,den_c,settling_band)
        return ControllerDesign(kind,sd,float(z),[-float(z)]*nzeros,float(kc),float(kp),float(ki),float(kd),
            [float(v) for v in num_c],[float(v) for v in den_c],float(base_phase),float(required),float(phi),
            plant,ctrl_contrib,poles,metrics)

    def design(self,controller_type,specs,auto_refine=False):
        pole,zeta,wn,sigma=resolve_specs(specs)

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

        design=self._design_at_pole(controller_type,pole,specs.settling_band)

        if auto_refine and specs.mp_max is not None and specs.ts_max is not None and design.metrics is not None:
            history=[]
            candidate=design
            factor=1.0
            for _ in range(60):
                m=candidate.metrics
                history.append({
                    "sigma":-candidate.desired_pole.real,
                    "kc":candidate.kc,
                    "z":candidate.zero_parameter,
                    "mp":float(m.overshoot_percent) if m and m.overshoot_percent is not None else math.nan,
                    "ts":float(m.settling_time) if m and m.settling_time is not None else math.nan,
                })
                if meets_specs(m,specs.mp_max,specs.ts_max):
                    design=candidate
                    break
                factor*=1.02
                sigma_try=sigma*factor
                candidate=self._design_at_pole(
                    controller_type,
                    desired_pole(zeta,sigma_try/zeta),
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
        return design
