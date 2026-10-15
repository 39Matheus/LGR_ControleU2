"""Camada de domínio da calculadora."""
from .lgr import AnalisadorLGR, VetorAnalise
from .controllers import ControllerDesigner
from .discretization import discretize_project,discretize_transfer
from .specs import DesignSpecs

__all__=[
    "AnalisadorLGR","VetorAnalise","ControllerDesigner","DesignSpecs",
    "discretize_transfer","discretize_project",
]
