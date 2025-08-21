r"""
Module defines model parameters.
At the moment, two main models are supported:
    * cosmo-kr

    *cosmo-urry.

Other model parameter can be easily defined here such as M2 from Lindorff-larsen

Tesei, G., Schulze, T. K., Crehuet, R., & Lindorff-larsen, K. (2021).
Accurate model of liquid-liquid phase behaviour of intrinsically-disordered proteins from data-driven optimization of
single-chain properties. BioRxiv, 1–9.
"""

from .model_parameters import parameters
from .globular_model_parameters import globular_parameters