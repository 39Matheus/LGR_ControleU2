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
