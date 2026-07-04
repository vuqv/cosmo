r"""Model parameters for the COSMO coarse-grained models.

At the moment the main supported models are:

* ``hps_kr`` -- HPS with the Kapcha-Rossy hydropathy scale.
* ``hps_urry`` -- HPS with the Urry hydropathy scale.

Other model parameters can be easily defined here (e.g. the M2 model of Tesei,
Schulze, Crehuet & Lindorff-Larsen, *Nat. Comput. Sci.* 2021, from data-driven
optimization of single-chain properties).
"""

from .model_parameters import parameters
