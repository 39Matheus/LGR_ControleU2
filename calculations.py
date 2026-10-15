"""Compatibilidade com a versão anterior.

O núcleo do LGR foi movido para core.lgr. Este arquivo evita quebrar imports
antigos como: from calculations import AnalisadorLGR.
"""
from core.lgr import AnalisadorLGR,VetorAnalise

__all__=["AnalisadorLGR","VetorAnalise"]
