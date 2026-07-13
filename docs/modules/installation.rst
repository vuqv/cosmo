How to install
==============

COSMO is a Python package built on `OpenMM <https://openmm.org/>`_. It depends on
OpenMM plus NumPy, ParmEd, MDAnalysis, mdtraj, and pandas, which are best installed
from ``conda-forge``. The steps below target Linux.

Requirements
------------

* **OpenMM >= 7.7** (choose a CUDA toolkit compatible with your NVIDIA driver).
  ``getStepCount()`` is unreliable before 7.7 and is required to restart simulations;
  OpenMM 8.2 is recommended for better performance.
* **ParmEd**, and the analysis dependencies listed in ``requirements.txt``.

Steps
-----

1. **Create and activate a fresh conda/mamba environment** with a recent Python
   (e.g. 3.9+)::

       conda create -n cosmo python=3.10 && conda activate cosmo

2. **Install OpenMM 7.7 or later**::

       conda install -c conda-forge openmm=7.7 cudatoolkit=10.2

   Conda may pick a newer CUDA toolkit by default. Choose a version compatible with
   your NVIDIA driver.

3. **Install the remaining dependencies.** The full list (NumPy, ParmEd, MDAnalysis,
   mdtraj, pandas, plus analysis-only extras used by ``scripts/``) is in
   ``requirements.txt``::

       conda install --file requirements.txt

   If you use the pip install in step 5, COSMO's runtime dependencies are pulled in
   automatically from ``pyproject.toml``; ``requirements.txt`` still adds the
   analysis-only packages.

4. **Download this repository** to a target path, for example ``PATH_TO_CODE/cosmo/``.

5. **Install COSMO.** There are two main ways.

Two ways to install
-------------------

**(I) Add to** ``PYTHONPATH`` **(no install)**
    Add the repository root to your Python path (persist it in ``.bashrc``)::

        export PYTHONPATH=$PYTHONPATH:PATH_TO_CODE/cosmo/

    This exposes ``import cosmo`` and the module form of every tool
    (``python -m cosmo.mdrun``, etc.), but **not** the console commands.

**(II) Install with pip**
    From the repository root (the directory with ``pyproject.toml``). This
    additionally registers the ``cosmo-mdrun`` and other console commands on your CLI:

    * editable (recommended for development), reflects source edits immediately::

          pip install -e .

    * regular install, copies the package into ``site-packages`` (source edits
      require a reinstall)::

          pip install .

Remember to replace ``PATH_TO_CODE`` with your actual path.

Console commands
----------------

``pip install`` registers these entry points (each also has a module form,
``python -m <module>``):

.. list-table::
   :header-rows: 1
   :widths: 22 26 52

   * - Command
     - Module
     - Purpose
   * - ``cosmo-mdrun``
     - ``cosmo.mdrun``
     - Run an IDP / condensate simulation from an ``md.ini``.
   * - ``cosmo-csp``
     - ``cosmo.csp.protocol``
     - Co-translational synthesis on an explicit CG ribosome.
   * - ``cosmo-cylinder``
     - ``cosmo.csp.cylinder``
     - Co-translational synthesis through an analytic tunnel.
   * - ``cosmo-csp-movie``
     - ``cosmo.csp.movie``
     - Stitch per-residue/-stage synthesis trajectories into a VMD movie.
   * - ``cosmo-make-mrna``
     - ``cosmo.csp.synth_mrna``
     - Pre-generate a fastest/slowest synonymous-codon mRNA.

(``cosmo-simulation`` is a back-compat alias for ``cosmo-mdrun``.)

Verify
------

::

    cosmo-mdrun -h     # prints help if the console command is installed
    python -c "import cosmo; print(cosmo.__version__)"
