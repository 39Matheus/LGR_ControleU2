"""Fechamento de malha, polos e verificação temporal."""
import numpy as np
from .models import SimulationMetrics
try:
    import control as ctrl
except ImportError:
    ctrl=None

def _pad_add(a,b):
    a=np.asarray(a,dtype=float); b=np.asarray(b,dtype=float)
    if len(a)<len(b): a=np.pad(a,(len(b)-len(a),0))
    elif len(b)<len(a): b=np.pad(b,(len(a)-len(b),0))
    return a+b

def closed_loop_polynomial(num_g,den_g,num_h,den_h,num_c,den_c):
    d=np.polymul(np.polymul(den_c,den_g),den_h)
    n=np.polymul(np.polymul(num_c,num_g),num_h)
    return _pad_add(d,n)

def closed_loop_poles(num_g,den_g,num_h,den_h,num_c,den_c):
    return [complex(x) for x in np.roots(closed_loop_polynomial(num_g,den_g,num_h,den_h,num_c,den_c))]

def closed_loop_transfer(num_g,den_g,num_h,den_h,num_c,den_c):
    fnum=np.polymul(num_c,num_g)
    den=closed_loop_polynomial(num_g,den_g,num_h,den_h,num_c,den_c)
    num=np.polymul(fnum,den_h)
    return np.trim_zeros(num,"f"),np.trim_zeros(den,"f")

def step_metrics(num_g,den_g,num_h,den_h,num_c,den_c,settling_band=0.02):
    poles=closed_loop_poles(num_g,den_g,num_h,den_h,num_c,den_c)
    stable=bool(poles) and all(p.real < -1e-8 for p in poles)
    if not stable:
        return SimulationMetrics(stable=False,settling_band=settling_band)
    if ctrl is None:
        return SimulationMetrics(stable=True,settling_band=settling_band)
    num,den=closed_loop_transfer(num_g,den_g,num_h,den_h,num_c,den_c)
    sys=ctrl.TransferFunction(num,den)
    slowest=min(abs(p.real) for p in poles if p.real < -1e-8)
    t_end=min(max(12.0/slowest,10.0),180.0)
    t=np.linspace(0,t_end,7000)
    tout,yout=ctrl.step_response(sys,T=t)
    y=np.real(np.asarray(yout).squeeze())
    try:
        final=float(np.real(ctrl.dcgain(sys)))
        if not np.isfinite(final): final=float(y[-1])
    except Exception:
        final=float(y[-1])
    peak=float(np.max(y))
    if abs(final)>1e-10:
        mp=max(0.0,(peak-final)/abs(final)*100.0)
        tol=settling_band*abs(final)
    else:
        mp=None; tol=settling_band
    outside=np.flatnonzero(np.abs(y-final)>tol)
    if len(outside)==0: ts=0.0
    elif outside[-1]+1<len(tout): ts=float(tout[outside[-1]+1])
    else: ts=float("inf")
    return SimulationMetrics(True,final,peak,mp,ts,settling_band)

def meets_specs(metrics,mp_max=None,ts_max=None):
    if metrics is None or not metrics.stable: return False
    if mp_max is not None and (metrics.overshoot_percent is None or metrics.overshoot_percent>mp_max+1e-6):
        return False
    if ts_max is not None and (metrics.settling_time is None or not metrics.settling_time<ts_max):
        return False
    return True
