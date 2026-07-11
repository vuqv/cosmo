Synthesis control options
=========================

Both synthesis runners are configured from a ``.ini`` control file with a single
``[OPTIONS]`` section:

* **``cosmo-csp``** — the coarse-grained-ribosome runner (``csp.ini``), read by
  :func:`cosmo.csp.protocol.read_csp_config`;
* **``cosmo-cylinder``** — the analytic-tunnel runner (``cylinder.ini``), read by
  :func:`cosmo.csp.cylinder.read_cylinder_config`.

Both return a populated :class:`~cosmo.csp.core.RunParams`, so **most keys are shared**.
**This page is the single reference for every control key** — grouped into *shared*,
*coarse-grained-ribosome-only*, and *cylinder-only* sets below; the
:doc:`continuous_synthesis` and :doc:`cylinder_synthesis` pages cover the *physics* and
point back here for the keys.

* Comments: inline or on their own line, starting with ``;`` or ``#``.
* Keyword and value are separated by ``=`` or ``:``.
* Every option has a default **except** the ones marked *required*.
* Units are OpenMM defaults — nm, ps, kJ/mol, K, kJ/mol/nm² — and **dwell times are in
  seconds**. Integers may use ``_`` digit separators (``200_000``).

Running a synthesis
-------------------

Once cosmo is installed (``pip install -e .`` from the repo root):

.. code-block:: bash

    cosmo-csp      -f csp.ini         # coarse-grained-ribosome runner
    cosmo-cylinder -f cylinder.ini    # analytic-tunnel runner

    # stitch the per-stage trajectories into one VMD movie afterwards
    cosmo-csp-movie -o synth_out --ribosome 4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb

A GPU is recommended for ``cosmo-csp``: the truncated ribosome adds several thousand
rigid beads.

Example control files
---------------------

``csp.ini`` (coarse-grained ribosome):

.. code-block::

        [OPTIONS]
        ; --- inputs (pdb_file, ribosome are required) ---
        pdb_file = asyn.pdb              ; native PDB (all-atom OR Cα-only CG); CG model built from it
        ribosome = 4v9d_..._cg_trunc.pdb ; truncated CG ribosome (P-/A-anchors + rigid scenery)
        model    = hps_urry              ; nascent IDP force field (any model; hps_kr is the default)

        ; --- length schedule ---
        L0 = 1 ; L_max =                 ; final length blank -> full residue count

        ; --- codon-resolved kinetics ---
        mrna = mrna.txt ; codon_times = trans_times.txt ; scale_factor = 4331293
        time_stage_1 = 0.000340 ; time_stage_2 = 0.004201  ; [CSP only]
        random_seed = 1

        ; --- integrator / ribosome mechanics ---
        dt = 0.01 ; ref_t = 300 ; tau_t = 0.01 ; nstout = 50
        constraints = None ; restraint_k = 83680 ; minimize = yes
        ; tunnel_wall = yes                ; one-sided tunnel floor (default on)

        ; --- post-synthesis + resume + hardware ---
        ejection_steps = 20000 ; dissociation_steps = 0 ; resume = auto
        device = GPU ; ppn = 4 ; outdir = synth_out

``cylinder.ini`` (analytic tunnel) — the shared keys are identical; it drops ``ribosome``
and adds the ``tunnel_*`` geometry:

.. code-block::

        [OPTIONS]
        pdb_file = asyn.pdb ; model = hps_kr        ; (no `ribosome` PDB)
        L0 = 1 ; L_max = ; mrna = mrna.txt ; codon_times = trans_times.txt
        scale_factor = 4331293 ; random_seed = 1
        constraints = None ; restraint_k = 83680 ; minimize = yes
        dt = 0.01 ; ref_t = 300 ; tau_t = 0.01 ; nstout = 50

        ; --- analytic exit tunnel (cylinder only) ---
        tunnel_radius = 0.9 ; tunnel_length = 10.0 ; tunnel_x_lo = 0.0
        tunnel_center = 0.0, 0.0 ; tunnel_k = 8368 ; tunnel_mouth_round = 0.2

        ejection_steps = 300000 ; dissociation_steps = 0 ; resume = auto
        device = GPU ; ppn = 4 ; outdir = synth_out


Option reference
++++++++++++++++

"Required = yes" means the run cannot proceed without it. ``—`` in the *Default* column
means there is no default. The keys are grouped into three sets — all in the same single
``[OPTIONS]`` section:

* **Shared** — accepted (with identical meaning) by *both* runners.
* **Coarse-grained ribosome runner only** (``cosmo-csp`` / ``csp.ini``).
* **Cylinder runner only** (``cosmo-cylinder`` / ``cylinder.ini``).

Shared options (both runners)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Inputs & length schedule.* A length-``L`` model is ``buildCoarseGrainModel`` on residues
``1..L`` of ``pdb_file`` — cosmo's chain is sequence-based, so there is **no** ``domain_def``
/ ``stride_output_file`` / native-contact step (those are topo's Gō-model inputs).

.. list-table::
   :header-rows: 1
   :widths: 20 14 12 16 38

   * - Option
     - Type
     - Required
     - Default
     - Description
   * - ``pdb_file``
     - str
     - **yes**
     - ``—``
     - Native PDB of the target protein — **all-atom *or* a Cα-only CG structure both work** (only the Cα positions and residue names are read; no STRIDE / native-contact step). The chain at length ``L`` uses residues ``1..L``.
   * - ``model``
     - str
     - no
     - ``hps_kr``
     - Nascent IDP force field. **Any model works** (``hps_kr`` / ``hps_urry`` / ``mpipi``). The ribosome↔nascent 12-10-6 excluded volume uses the **model-independent** O'Brien ``Rmin/2`` radii (decoupled from the force field); the IDP↔IDP interaction is whatever the model provides (Ashbaugh–Hatch or Wang–Frenkel).
   * - ``L0``
     - int
     - no
     - ``1``
     - First nascent length to synthesize.
   * - ``L_max``
     - int
     - no
     - full length
     - Final nascent length. Omit/blank for the whole chain. Must satisfy ``1 ≤ L0 ≤ L_max ≤ N_full``.
   * - ``outdir``
     - str
     - no
     - ``synth_out``
     - Output root; each residue writes one ``L_<L>/`` folder (CSP: shared ``traj.psf`` + per-stage ``traj_s<1,2,3>.dcd``; cylinder: a single ``traj.dcd``).

*Codon kinetics & schedule length.*

.. list-table::
   :header-rows: 1
   :widths: 20 12 12 14 42

   * - Option
     - Type
     - Required
     - Default
     - Description
   * - ``mrna``
     - str
     - for per-codon timing
     - ``—``
     - mRNA sequence file (one codon per residue + one stop), **or** ``fastest``/``slowest``/``median`` to auto-build a synonymous-codon mRNA. Required unless ``codon_times`` is a number. A real filename must not be ``fastest``/``slowest``/``median``. See the notes below.
   * - ``codon_times``
     - str or float
     - for per-codon timing
     - ``—``
     - A **path** to a ``CODON  seconds  amino_acid`` table selects **per-codon** timing (no bundled default — pick one under ``assets/csp/codon_dwell_times/``); a **positive number of seconds** selects **uniform** timing (no ``mrna`` needed). A table filename must not be a bare number.
   * - ``scale_factor``
     - float
     - no
     - ``4331293``
     - In-vivo-seconds → in-silico-ns compression (``t_sim_ns = t_s · 1e9 / scale_factor``). Larger ⇒ fewer MD steps per residue ⇒ faster, preserving the relative fast/slow-codon timing.
   * - ``random_seed``
     - int
     - no
     - ``—`` (nondet.)
     - Seed for the first-passage-time sampler — makes the whole dwell schedule reproducible.
   * - ``ribosome_traffic``
     - bool
     - no
     - ``no``
     - Apply the ribosome-traffic (polysome) dwell-time correction, if available, on top of the per-codon kinetics. Off by default (single-ribosome timing).
   * - ``initiation_rate``
     - float
     - no
     - ``0.083333``
     - Translation initiation rate (1/s). Consumed **only** when ``ribosome_traffic = yes``.
   * - ``max_steps_per_stage``
     - int
     - no
     - ``—`` (uncapped)
     - **Testing only.** Upper clamp on each stage's MD step count. See the warning below.
   * - ``min_steps_per_stage``
     - int
     - no
     - ``1``
     - **Testing only.** Lower clamp on each stage's MD step count.

.. warning::

   **``max_steps_per_stage`` / ``min_steps_per_stage`` are testing-only knobs.** They
   clamp the MD step count so examples finish quickly, which **breaks the physical
   timescale mapping**. In production leave them **unset** so step counts come entirely
   from the kinetics (``scale_factor``, the codon times, ``dt``). The sampled dwell
   **times in seconds** are always written to ``dwell_times.dat``.

*Integrator & MD mechanics* (the shared per-length engine,
:class:`~cosmo.csp.core.RunParams`).

.. list-table::
   :header-rows: 1
   :widths: 22 14 16 48

   * - Option
     - Type
     - Default
     - Description
   * - ``dt``
     - float [ps]
     - ``0.01``
     - Integration timestep. cosmo uses flexible harmonic bonds and soft HPS/mpipi potentials.
   * - ``ref_t``
     - float [K]
     - ``300``
     - Langevin temperature. Set it to match your ``codon_times`` table's temperature.
   * - ``tau_t``
     - float [ps⁻¹]
     - ``0.01``
     - Langevin friction coefficient.
   * - ``nstout``
     - int
     - ``50``
     - Trajectory (DCD) and log output interval, in steps.
   * - ``constraints``
     - str
     - ``None``
     - Bond treatment: ``None`` (flexible harmonic bonds, the default — backbone flexibility is physically meaningful for disordered chains) or ``AllBonds`` (rigid distance constraints on the CA/P pseudo-bonds, allowing a larger timestep). Constraints act only on the pseudo-bonds; the non-bonded potentials are unaffected. The equilibrium-bond seeding keeps ``AllBonds`` stable too.
   * - ``restraint_k``
     - float [kJ/mol/nm²]
     - ``83680``
     - Stiffness of the C-terminus harmonic restraint to the PTC target (= 200 kcal/mol/Å²). Its ``k`` is a per-particle parameter so it coexists with the tunnel wall's global ``k``.
   * - ``minimize``
     - bool
     - ``yes``
     - Energy-minimize each seeded structure before running that stage's MD.
   * - ``device``
     - str
     - ``CPU``
     - Compute platform: ``CPU`` or ``GPU`` (CUDA). GPU recommended for ``cosmo-csp``.
   * - ``ppn``
     - int
     - ``1``
     - Number of CPU threads (``device = CPU``).

*Post-synthesis phases & resume.*

.. list-table::
   :header-rows: 1
   :widths: 22 10 12 56

   * - Option
     - Type
     - Default
     - Description
   * - ``ejection_steps``
     - int
     - ``0``
     - Post-synthesis ejection phase length (steps); ``0`` = skip. Releases the C-terminus restraint so the chain diffuses out (+x).
   * - ``dissociation_steps``
     - int
     - ``0``
     - Post-synthesis dissociation phase length (steps); ``0`` = skip. A further free run away from the ribosome.
   * - ``resume``
     - str
     - ``auto``
     - Resume policy for interrupted runs: ``auto`` (resume iff a ``progress.log`` is present under ``outdir``, else fresh), ``yes`` (require a resumable run, else error), ``no`` (always fresh). CLI ``--fresh`` / ``--no-resume`` forces ``no``. See :doc:`synthesis_resume`.

Coarse-grained ribosome runner only (``cosmo-csp``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These keys are read only from ``csp.ini``. (``time_stage_1`` / ``time_stage_2`` are
*accepted* by the cylinder INI but have **no effect** there — with a single MD segment
per residue the whole codon dwell is one segment, not a three-way split.)

.. list-table::
   :header-rows: 1
   :widths: 22 12 14 52

   * - Option
     - Type
     - Default
     - Description
   * - ``ribosome``
     - str
     - **required**
     - Truncated CG ribosome PDB — the P-/A-anchors and the rigid (mass-0) scenery. Supplying it *is* the signal to load it as rigid (no ``rigid_ribosome`` key). Must carry tRNA beads under the fixed names (segids ``PtR``/``AtR``, resid 76, beads ``R``/``P``/``BR2``).
   * - ``time_stage_1``
     - float [s]
     - ``0.00034``
     - Mean peptidyl-transfer (stage 1) dwell; the actual dwell is an exponential draw with this mean.
   * - ``time_stage_2``
     - float [s]
     - ``0.004201``
     - Mean translocation (stage 2) dwell. Stage 3's mean is the remainder ``τ(next codon) − time_stage_1 − time_stage_2``.
   * - ``tunnel_wall``
     - bool
     - ``yes``
     - Apply the one-sided half-harmonic tunnel wall (a floor below the synthesis point). The plane ``x₀`` is **auto-derived** from the ribosome (``min(A-target.x, P-target.x)``) and the stiffness is a fixed constant — neither is a key; only this on/off toggle.

Cylinder runner only (``cosmo-cylinder``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The analytic exit tunnel — a bore of radius ``r`` drilled through an infinite wall — that
replaces the explicit ribosome beads. Read only from ``cylinder.ini``.

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Option
     - Default
     - Description
   * - ``tunnel_radius``
     - ``0.9``
     - Bore radius ``r`` (nm); ~3 CG beads wide.
   * - ``tunnel_length``
     - ``10.0``
     - Bore length (nm). The exit face sits at ``x_exit = tunnel_x_lo + tunnel_length`` (derived, not a key).
   * - ``tunnel_x_lo``
     - ``0.0``
     - PTC / closed end of the bore (nm); the C-terminus is restrained on-axis here.
   * - ``tunnel_center``
     - ``0.0, 0.0``
     - Tunnel axis ``(y0, z0)`` (nm). The axis runs along +x.
   * - ``tunnel_k``
     - ``8368``
     - Wall stiffness (kJ/mol/nm² = 20 kcal/mol/Å²).
   * - ``tunnel_mouth_round``
     - ``0.2``
     - Mouth-corner fillet radius ``rho`` (nm); rounds the 90° inner corner so the potential is continuous.

.. note::

   Boolean options accept ``yes``/``no``, ``true``/``false``, ``1``/``0``. The
   ``trna_tether`` field exists on ``RunParams`` but is **forced off** by the CSP runner —
   CSP needs the switchable A↔P **position** restraint, so the O'Brien tRNA tether (which
   targets only the P-site) is incompatible with the 3-stage translocation.


Notes on individual options
+++++++++++++++++++++++++++

Required inputs (``pdb_file`` / ``ribosome``)
    ``pdb_file`` supplies the sequence (the length-``L`` model is
    ``buildCoarseGrainModel`` on residues ``1..L`` — no STRIDE, no contact map). **It may
    be an all-atom PDB or a Cα-only CG structure** — only the Cα positions and residue
    names are read (unlike the sibling ``topo`` Gō model, which needs an all-atom
    structure for STRIDE + heavy-atom native contacts). ``ribosome`` (``cosmo-csp`` only)
    provides the P-/A-anchors and rigid scenery; it must carry tRNA beads under the fixed
    names (segids ``PtR``/``AtR``, resid 76, beads ``R``/``P``/``BR2``) or anchor lookup —
    and the always-on PTC-geometry optimization — fail.

Per-codon vs. uniform timing (``mrna`` / ``codon_times``)
    ``codon_times`` selects the timing mode by its **value type**: a **path** → per-codon
    timing (``mrna`` required); a **positive number of seconds** → uniform timing (no
    ``mrna``). There is no bundled default — per-codon timing needs an explicit table path
    (pick one under ``assets/csp/codon_dwell_times/``). A codon-time table filename must
    not be a bare number. Setting ``mrna = fastest``, ``slowest`` or ``median`` auto-builds
    a synonymous-codon mRNA from the protein + table (see :doc:`codon_dwell_times`).

Kinetics and step counts (``scale_factor`` / ``time_stage_1`` / ``time_stage_2``)
    In ``cosmo-csp`` each residue is added over three sub-stages whose dwell times are
    exponential draws with means ``time_stage_1`` (peptidyl transfer), ``time_stage_2``
    (translocation), and ``τ(next codon) − time_stage_1 − time_stage_2`` (the decoding
    wait). Those seconds map to MD steps through ``scale_factor`` and ``dt``.
    **``cosmo-cylinder`` uses a single segment per residue**, so the whole codon dwell
    ``τ`` is one MD segment and ``time_stage_1`` / ``time_stage_2`` are ignored.

Ribosome ↔ nascent excluded volume
    The rigid ribosome interacts with the chain via the O'Brien **12-10-6** excluded
    volume (sum rule ``R = Rmin/2ᵢ + Rmin/2ⱼ``, ``ε = 0.000132`` kcal/mol) plus
    Debye–Hückel electrostatics. The per-bead ``Rmin/2`` values are **model-independent**
    steric radii (standalone ``OBRIEN_RMIN_2_NM`` / ``OBRIEN_RNA_RMIN_2_BEADS`` tables in
    ``cosmo.parameters.model_parameters``), decoupled from the force field, so the excluded
    volume is identical for every nascent model.

Tunnel wall vs. analytic tunnel (``tunnel_wall`` / ``tunnel_*``)
    ``cosmo-csp`` uses a one-sided half-harmonic **wall** (``tunnel_wall``) whose plane is
    auto-placed and stiffness fixed, so only the on/off toggle is exposed. ``cosmo-cylinder``
    instead uses a fully **analytic tunnel** — a cylindrical bore in an infinite wall —
    whose geometry is set by the ``tunnel_radius`` / ``tunnel_length`` / ``tunnel_x_lo`` /
    ``tunnel_center`` / ``tunnel_k`` / ``tunnel_mouth_round`` keys (see
    :doc:`cylinder_synthesis`).

Post-synthesis phases (``ejection_steps`` / ``dissociation_steps``)
    After the last residue, ``ejection_steps > 0`` runs a phase with the C-terminus
    restraint **released** (ribosome/tunnel still present), so the finished chain diffuses
    out. ``dissociation_steps > 0`` continues the free protein away. Both write their own
    output folders; ``0`` skips the phase.

Output layout, resume and movies
    Each residue writes one ``<outdir>/L_<L>/`` folder (CSP: shared ``traj.psf`` +
    per-stage ``traj_s<1,2,3>.dcd`` and a single ``traj_final.pdb``; cylinder: a single
    ``traj.dcd`` + ``traj_final.pdb``), plus a per-residue ``dwell_times.dat`` (the
    schedule, with a ``#PTC`` header for CSP) and an append-only ``progress.log``.
    Re-invoking a runner on an interrupted ``outdir`` **resumes** from the last completed
    residue (``resume = auto``); see :doc:`synthesis_resume`. Stitch a movie with
    ``cosmo-csp-movie -o <outdir>``; set ``COSMO_CSP_VERBOSE=1`` for the full per-stage
    banners.
