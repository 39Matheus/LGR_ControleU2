"""Camada de domínio da calculadora."""
from .lgr import AnalisadorLGR, VetorAnalise
from .controllers import ControllerDesigner
from .specs import DesignSpecs

__all__=["AnalisadorLGR","VetorAnalise","ControllerDesigner","DesignSpecs"]
