.. COSMO citation / how-to-cite page.

How to cite COSMO
=================

If COSMO contributed to work you are publishing, please cite it. Citing the
software lets others find the exact model and version you used, and credits the
people who developed the underlying force fields.

.. tip::

   The repository ships a :download:`CITATION.cff <../CITATION.cff>` file, so on
   GitHub you can click **“Cite this repository”** to export formatted APA and
   BibTeX automatically. Keep that file and this page in step whenever the
   citation details change.


Primary citation
----------------

Cite the paper that introduces this coarse-grained toolkit:

   Vu, Q. V.; Sitarik, I.; Li, M. S.; O'Brien, E. P. *Noncovalent Lasso
   Entanglements are Common in Experimentally Derived Intrinsically Disordered
   Protein Ensembles and Strongly Influenced by Protein Length and Charge.*
   **J. Phys. Chem. B** 129, 4682–4691 (2025).
   https://doi.org/10.1021/acs.jpcb.5c01260

.. code-block:: bibtex

   @article{cosmo,
     author  = {Vu, Quyen V. and Sitarik, Ian and Li, Mai Suan and O'Brien, Edward P.},
     title   = {Noncovalent Lasso Entanglements are Common in Experimentally Derived
                Intrinsically Disordered Protein Ensembles and Strongly Influenced by
                Protein Length and Charge},
     journal = {The Journal of Physical Chemistry B},
     year    = {2025},
     volume  = {129},
     pages   = {4682--4691},
     doi     = {10.1021/acs.jpcb.5c01260}
   }

The software itself is archived on `Zenodo <https://zenodo.org/>`_. Cite the
**concept DOI** (always the latest release) alongside the paper:

   Vu, Q. (2026). *COSMO: COarse-grained Simulation of intrinsically disordered
   prOteins with OpenMM* (Version 2026.1) [Computer software]. Zenodo.
   https://doi.org/10.5281/zenodo.21361272

To pin the exact version you ran, use the version-specific DOI for 2026.1 instead
(https://doi.org/10.5281/zenodo.21361273), and record the version
(``import cosmo; print(cosmo.__version__)`` or ``pyproject.toml``) and commit
(``git rev-parse --short HEAD``) — different versions can produce different numbers.


The force fields COSMO implements
---------------------------------

COSMO runs **sequence-based, one-bead-per-residue** intrinsically-disordered-protein
force fields. **Cite the model(s) you actually run** (set by the ``model`` key).

HPS hydropathy-scale family (``hps_urry`` / ``hps_kr`` / ``hps_ss``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Dignon, G. L.; Zheng, W.; Kim, Y. C.; Best, R. B.; Mittal, J. *Sequence
  Determinants of Protein Phase Behavior from a Coarse-Grained Model.* PLoS
  Comput. Biol. **14** (1), e1005941 (2018).
  https://doi.org/10.1371/journal.pcbi.1005941
* Regy, R. M.; Thompson, J.; Kim, Y. C.; Mittal, J. *Improved Coarse-Grained Model
  for Studying Sequence-Dependent Phase Separation of Disordered Proteins.*
  Protein Sci. **30** (7), 1371–1379 (2021). https://doi.org/10.1002/pro.4094
* Rizuan, A.; Jovic, N.; Phan, T. M.; Kim, Y. C.; Mittal, J. *Developing Bonded
  Potentials for a Coarse-Grained Model of Intrinsically Disordered Proteins.*
  J. Chem. Inf. Model. **62** (18), 4474–4485 (2022).
  https://doi.org/10.1021/acs.jcim.2c00450

Mpipi model (``mpipi``)
~~~~~~~~~~~~~~~~~~~~~~~~

* Joseph, J. A.; Reinhardt, A.; Aguirre, A.; Chew, P. Y.; Russell, K. O.;
  Espinosa, J. R.; Garaizar, A.; Collepardo-Guevara, R. *Physics-Driven
  Coarse-Grained Model for Biomolecular Phase Separation with Near-Quantitative
  Accuracy.* Nat. Comput. Sci. **1** (11), 732–743 (2021).
  https://doi.org/10.1038/s43588-021-00155-3


Co-translational synthesis (CSP)
--------------------------------

If you use the ``cosmo-csp`` / ``cosmo-cylinder`` runners, also cite the O'Brien-lab
per-codon, three-stage continuous-synthesis protocol they reproduce:

* Jiang, Y. *et al.* *How synonymous mutations alter enzyme structure and function
  over long timescales.* Nat. Chem. **15**, 308–318 (2023).
  https://doi.org/10.1038/s41557-022-01091-z

**Codon dwell-time datasets.** When you run synthesis with a species dwell-time
table (see :doc:`usage/codon_dwell_times`), cite the source dataset:

* *E. coli* — Fluitt, A.; Pienaar, E.; Viljoen, H. Comput. Biol. Chem. **31**,
  335–346 (2007). https://doi.org/10.1016/j.compbiolchem.2007.07.003
* Yeast — Gardin, J. *et al.* eLife **3**, e03735 (2014).
  https://doi.org/10.7554/eLife.03735
* Human — Gobet, C. *et al.* PNAS **117** (17), 9630–9641 (2020).
  https://doi.org/10.1073/pnas.1918145117
* *N. crassa* — Yang, Q. *et al.* Nucleic Acids Res. **47** (17), 9243–9258
  (2019). https://doi.org/10.1093/nar/gkz710


MD engine (always)
------------------

All dynamics run on OpenMM:

* Eastman, P. *et al.* *OpenMM 7: Rapid development of high performance algorithms
  for molecular dynamics.* PLoS Comput. Biol. **13** (7), e1005659 (2017).
  https://doi.org/10.1371/journal.pcbi.1005659


Questions
---------

For anything not covered here — collaboration, a preprint DOI, or how to cite a
specific analysis — open an issue on the
`GitHub repository <https://github.com/vuqv/cosmo/>`_.
