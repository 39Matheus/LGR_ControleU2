from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class TrigContribution:
    label: str
    singularity: complex
    dx: float
    dy: float
    angle_deg: float
    kind: str
    role: str="plant"


@dataclass
class SimulationMetrics:
    stable: bool
    final_value: Optional[float]=None
    peak_value: Optional[float]=None
    overshoot_percent: Optional[float]=None
    settling_time: Optional[float]=None
    settling_band: Optional[float]=None


@dataclass
class ControllerDesign:
    controller_type: str
    desired_pole: complex
    zero_parameter: float
    zero_locations: List[float]
    kc: float
    kp: float
    ki: float
    kd: float
    controller_num: List[float]
    controller_den: List[float]
    base_phase_deg: float
    required_phase_deg: float
    zero_angle_deg: float
    plant_contributions: List[TrigContribution]=field(default_factory=list)
    controller_contributions: List[TrigContribution]=field(default_factory=list)
    closed_loop_poles: List[complex]=field(default_factory=list)
    metrics: Optional[SimulationMetrics]=None
    refinement_history: List[Dict[str,float]]=field(default_factory=list)
    specification_summary: Dict[str,Any]=field(default_factory=dict)

    # Extensões para controladores genéricos / PID com zeros distintos.
    pole_locations: List[float]=field(default_factory=list)
    solved_singularity_kind: str="zero"
    solved_singularity_parameter: Optional[float]=None
    solved_singularity_angle_deg: Optional[float]=None
    solved_singularity_count: int=1
    controller_latex: str=""
    controller_text: str=""
    controller_config: Dict[str,Any]=field(default_factory=dict)


@dataclass
class DiscretizationResult:
    method: str
    sample_time: float
    target_label: str
    continuous_num: List[float]
    continuous_den: List[float]
    discrete_num: List[float]
    discrete_den: List[float]
    continuous_latex: str
    substitution_latex: str
    substituted_latex: str
    discrete_latex: str
    difference_equation_latex: str
    poles: List[complex]=field(default_factory=list)
    zeros: List[complex]=field(default_factory=list)
    notes: List[str]=field(default_factory=list)
