Synthesis control options
=========================

Continuous-synthesis (co-translational) runs are configured from a separate
``.ini`` control file (e.g. ``csp.ini``), read by
:func:`cosmo.csp.protocol.read_csp_config`, which returns a
:class:`~cosmo.csp.protocol.CSPConfig` (inputs, length schedule, and a populated
:class:`~cosmo.csp.core.RunParams`). This is the synthesis analogue of the
single-chain :doc:`simulation_control` page: the same ``[OPTIONS]`` /
``key = value`` grammar, a different set of options.

* Comments: inline or on their own line, starting with ``;`` or ``#``.
* Keyword and value are separated by ``=`` or ``:``.
* Every option below has a default **except** the ones marked *required*.
* Units are OpenMM defaults — nm, ps, kJ/mol, K, kJ/mol/nm² — and **dwell times
  are in seconds**. Integers may use ``_`` digit separators (``200_000``).

For the *physics* behind these knobs (the three-stage elongation cycle, the
codon-resolved kinetics, the O'Brien 12-10-6 ribosome wall), see
:doc:`continuous_synthesis`. This page is the parameter reference.

Running a synthesis
-------------------

The runner lives in the package as :mod:`cosmo.csp`. Once cosmo is installed
(``pip install -e .`` from the repo root) either of these works:

.. code-block:: bash

    cosmo-csp -f csp.ini              # installed console command
    python -m cosmo.csp -f csp.ini    # module form

    # stitch the per-stage trajectories into one VMD movie afterwards
    cosmo-csp-movie -o synth_out --ribosome 4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb

A GPU is recommended: the truncated ribosome adds several thousand rigid beads.

Example ``csp.ini``:

.. code-block::

        [OPTIONS]
        ; --- inputs (pdb_file, ribosome are required) ---
        pdb_file = asyn.pdb              ; full native PDB; the CG model is built from it
        ribosome = 4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb   ; truncated CG ribosome (P-/A-anchors + rigid scenery)
        model    = hps_urry              ; nascent IDP force field (any model works; hps_kr is the default)

        ; --- length schedule ---
        L0    = 5            ; first nascent length (default 1)
        L_max = 8            ; final nascent length (default: full residue count)

        ; --- codon-resolved kinetics ---
        mrna         = mrna.txt        ; one codon per residue (required for per-codon timing);
                                       ; or "fastest"/"slowest"/"median" to auto-build a synonymous-codon mRNA
        codon_times  = trans_times.txt ; table path = per-codon (required, no bundled default);
                                       ; a positive number of s = uniform instead
        scale_factor = 4331293         ; in-vivo seconds -> in-silico ns compression (larger = faster)
        time_stage_1 = 0.000340        ; mean peptidyl-transfer dwell (s)
        time_stage_2 = 0.004201        ; mean translocation dwell (s)
        random_seed  = 1               ; reproducible first-passage-time schedule

        ; --- test clamps (leave UNSET in production; see the warning below) ---
        max_steps_per_stage = 40
        min_steps_per_stage = 20

        ; --- integrator / output ---
        dt     = 0.01        ; timestep (ps)
        ref_t  = 300         ; temperature (K)
        tau_t  = 0.01        ; Langevin friction (1/ps)
        nstout = 20          ; trajectory/log output interval (steps)

        ; --- ribosome / PTC mechanics ---
        ; PTC geometry is always optimized (A/P targets one peptide bond, ~0.38 nm,
        ; the selected model's bond length, apart and EV-clear) -- no knob. Bonds default
        ; to flexible (constraints = None;
        ; set constraints = AllBonds for rigid bonds).
        restraint_k = 83680  ; C-terminus harmonic restraint (kJ/mol/nm^2)
        minimize    = yes    ; energy-minimize each seeded structure before its MD
        ; tunnel_wall = yes  ; one-sided tunnel floor (default on; plane auto-derived)

        ; --- post-synthesis phases (steps; 0 = skip) ---
        ejection_steps     = 20000   ; release the restraint; chain diffuses out (+x)
        dissociation_steps = 0       ; free run away from the ribosome

        ; --- hardware / output ---
        device = GPU
        ppn    = 4
        outdir = synth_out


Parameter summary
+++++++++++++++++

"Required = yes" means the run cannot proceed without it. ``—`` in the *Default*
column means there is no default. Options are grouped by role; they all live in
the same single ``[OPTIONS]`` section.

Inputs & length schedule
~~~~~~~~~~~~~~~~~~~~~~~~~~

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
     - Full native PDB of the target protein. cosmo builds the one-bead-per-residue CG model from it; the nascent chain at length ``L`` uses residues ``1..L``.
   * - ``ribosome``
     - str
     - **yes**
     - ``—``
     - Truncated CG ribosome PDB. Source of the P-/A-anchors and the rigid (mass-0) scenery. Supplying it *is* the signal to load it as rigid — there is no ``rigid_ribosome`` key. Must carry the tRNA beads under the expected names (segids ``PtR``/``AtR``, resid 76, beads ``R``/``P``/``BR2``).
   * - ``model``
     - str
     - no
     - ``hps_kr``
     - Nascent IDP force field. **Any model works** (``hps_kr`` / ``hps_urry`` / ``mpipi``); ``hps_kr`` is merely the default. The ribosome↔nascent 12-10-6 excluded volume uses the **model-independent** O'Brien ``Rmin/2`` collision radii (standalone ``OBRIEN_RMIN_2_NM`` / ``OBRIEN_RNA_RMIN_2_BEADS`` tables in ``cosmo.parameters.model_parameters``), decoupled from the force field; the nascent IDP↔IDP interaction is whatever the selected model provides (Ashbaugh–Hatch or Wang–Frenkel).
   * - ``mrna``
     - str
     - for per-codon timing
     - ``—``
     - mRNA sequence file (raw nucleotides; one codon per residue plus one stop), **or** the keyword ``fastest`` / ``slowest`` / ``median`` to auto-build a synonymous-codon mRNA (each residue's fastest/slowest/median-dwell-time codon per the ``codon_times`` table, written next to the PDB). Required unless ``codon_times`` is a number (uniform timing). A real filename must not be ``fastest``/``slowest``/``median``.
   * - ``codon_times``
     - str or float
     - for per-codon timing
     - ``—``
     - Codon-timing key. A **path** to a ``CODON  seconds  amino_acid`` table = **per-codon** timing (required, no bundled default -- pick one under ``assets/csp/codon_dwell_times/``); a **positive number of seconds** = **uniform** timing (no ``mrna`` needed). A table filename must **not** be a bare number.
   * - ``L0``
     - int
     - no
     - ``1``
     - First nascent length to synthesize.
   * - ``L_max``
     - int
     - no
     - full length
     - Final nascent length. Must satisfy ``1 <= L0 <= L_max <= N_full``.
   * - ``outdir``
     - str
     - no
     - ``synth_out``
     - Output root; each residue writes one ``L_<L>/`` folder (per-stage ``traj_s<1,2,3>.dcd`` + a shared ``traj.psf``).

.. note::

   There is **no** ``domain_def``, ``stride_output_file``, or ``nascent_ev_radii``
   option: those are topo's Gō-model inputs. cosmo's chain is sequence-based, so a
   length-``L`` model is just ``buildCoarseGrainModel`` on residues ``1..L`` — no
   STRIDE, no native-contact map.

Kinetics & schedule length
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 12 12 14 42

   * - Option
     - Type
     - Required
     - Default
     - Description
   * - ``scale_factor``
     - float
     - no
     - ``4331293``
     - In-vivo-seconds → in-silico-ns compression (``t_sim_ns = t_s · 1e9 / scale_factor``). Larger ⇒ fewer MD steps per residue ⇒ faster, while preserving the relative fast/slow-codon timing.
   * - ``time_stage_1``
     - float [s]
     - no
     - ``0.00034``
     - Mean peptidyl-transfer dwell. The per-residue dwell is an exponential draw with this mean.
   * - ``time_stage_2``
     - float [s]
     - no
     - ``0.004201``
     - Mean translocation dwell. Stage 3's mean is the remainder ``τ(next codon) − time_stage_1 − time_stage_2``.
   * - ``random_seed``
     - int
     - no
     - ``—`` (nondet.)
     - Seed for the first-passage-time sampler (reproducible dwell schedule).
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
     - **Testing only.** Upper clamp on each stage's MD step count. Breaks the physical timescale mapping — leave **unset** in production.
   * - ``min_steps_per_stage``
     - int
     - no
     - ``1``
     - **Testing only.** Lower clamp on each stage's MD step count.
   * - ``ejection_steps``
     - int
     - no
     - ``0``
     - Post-synthesis ejection phase (steps); ``0`` = skip. Releases the C-terminus restraint so the chain diffuses out (+x).
   * - ``dissociation_steps``
     - int
     - no
     - ``0``
     - Post-synthesis dissociation phase (steps); ``0`` = skip. A further free run away from the ribosome.

.. warning::

   **``max_steps_per_stage`` / ``min_steps_per_stage`` are testing-only knobs.**
   They clamp the MD step count so examples finish quickly, which **breaks the
   physical timescale mapping**. In production leave them **unset** so step counts
   come entirely from the kinetics (``scale_factor``, the codon times, ``dt``). The
   sampled dwell **times in seconds** are always written to ``dwell_times.dat``.

Integrator, ribosome & PTC mechanics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These configure the shared per-length engine (:class:`~cosmo.csp.core.RunParams`,
consumed by :func:`~cosmo.csp.core.run_length`).

.. list-table::
   :header-rows: 1
   :widths: 22 12 16 50

   * - Option
     - Type
     - Default
     - Description
   * - ``dt``
     - float [ps]
     - ``0.01``
     - Integration timestep. (cosmo uses flexible harmonic bonds and soft HPS/mpipi potentials — there is no rigid-bond / dt-halving stability guard.)
   * - ``ref_t``
     - float [K]
     - ``300``
     - Langevin temperature.
   * - ``tau_t``
     - float [ps⁻¹]
     - ``0.01``
     - Langevin friction coefficient.
   * - ``nstout``
     - int
     - ``50``
     - Trajectory (DCD) and log output interval, in steps, for every stage.
   * - ``constraints``
     - str
     - ``None``
     - Bond treatment: ``None`` (flexible harmonic bonds, the default) or ``AllBonds`` (rigid distance constraints). Both are supported. ``None`` is the default because backbone flexibility is physically meaningful for intrinsically disordered chains; ``AllBonds`` pins every CA/P bond at its equilibrium length, removing the fast bond-stretch mode so a larger timestep can be used. Constraints act only on the pseudo-bonds — the non-bonded potentials (Ashbaugh-Hatch / Wang-Frenkel and, with a ribosome, the 12-10-6 excluded volume) are unaffected, so ``AllBonds`` does not by itself prevent a stiff-EV blow-up (the per-stage dt-halving guard handles that).
   * - ``restraint_k``
     - float [kJ/mol/nm²]
     - ``83680``
     - Stiffness of the C-terminus harmonic position restraint to the A/P target point (= 200 kcal/mol/Å²). Switching the target A→P is how translocation is reproduced. (Its ``k`` is a per-particle parameter so it coexists with the tunnel wall's global ``k``.)
   * - ``minimize``
     - bool
     - ``yes``
     - Energy-minimize each seeded structure before running that stage's MD.
   * - ``tunnel_wall``
     - bool
     - ``yes``
     - Apply the one-sided half-harmonic tunnel wall (a floor below the synthesis point). The plane ``x₀`` is **auto-derived** from the ribosome (``min(A-target.x, P-target.x)``) and the stiffness is a fixed model constant — neither is a key; only this on/off toggle.
   * - ``device``
     - str
     - ``CPU``
     - Compute platform: ``CPU`` or ``GPU`` (CUDA). A GPU is recommended given the rigid-ribosome bead count.
   * - ``ppn``
     - int
     - ``1``
     - Number of CPU threads (``device = CPU``).

.. note::

   Boolean options accept ``yes``/``no``, ``true``/``false``, ``1``/``0``. The
   ``trna_tether`` field exists on ``RunParams`` but is **forced off** by the CSP
   runner — CSP needs the switchable A↔P **position** restraint, so the O'Brien
   tRNA tether (which targets only the P-site) is incompatible with the 3-stage
   translocation.


Notes on individual options
+++++++++++++++++++++++++++

Required inputs (``pdb_file`` / ``ribosome``)
    A synthesis run cannot start without both. ``pdb_file`` supplies the sequence
    (the length-``L`` model is ``buildCoarseGrainModel`` on residues ``1..L`` — no
    STRIDE, no contact map). ``ribosome`` is the truncated CG ribosome providing
    the P-/A-anchors and the rigid scenery; it must carry tRNA beads under the fixed
    names (segids ``PtR``/``AtR``, resid 76, beads ``R``/``P``/``BR2``) or anchor
    lookup — and the always-on PTC-geometry optimization — fail.

Per-codon vs. uniform timing (``mrna`` / ``codon_times``)
    ``codon_times`` selects the timing mode by its **value type**: a **path** →
    per-codon timing (``mrna`` required); a **positive number of seconds** →
    uniform timing (no ``mrna``). There is no bundled default — per-codon timing needs
    an explicit table path (pick one under ``assets/csp/codon_dwell_times/``). A
    codon-time table filename must not be a bare number. Setting ``mrna = fastest``,
    ``slowest`` or ``median`` auto-builds a synonymous-codon mRNA from the protein + table
    (see :doc:`codon_dwell_times`).

Ribosome ↔ nascent excluded volume
    The rigid ribosome interacts with the chain via the O'Brien **12-10-6**
    excluded volume (sum combining rule ``R = Rmin/2ᵢ + Rmin/2ⱼ``, ``ε = 0.000132``
    kcal/mol) plus Debye–Hückel (Yukawa) electrostatics. The per-bead ``Rmin/2``
    values are inherited from topo and are **model-independent** steric radii: they
    live in the standalone ``OBRIEN_RMIN_2_NM`` (per-AA) and ``OBRIEN_RNA_RMIN_2_BEADS``
    (rRNA ``P``/``R``/``BR``) tables in ``cosmo.parameters.model_parameters``, decoupled
    from any force field, so the excluded volume is identical for every nascent model.

Post-synthesis phases (``ejection_steps`` / ``dissociation_steps``)
    After the last residue, ``ejection_steps > 0`` runs a phase with the C-terminus
    restraint **released** (rigid ribosome and tunnel wall still present), so the
    finished chain diffuses out of the tunnel. ``dissociation_steps > 0`` continues
    the free protein away from the ribosome. Both write their own output folders.

Output layout and progress log
    Every stage writes a standalone, **nascent-only** trajectory to
    ``<outdir>/L_<L>/`` (per-stage ``traj_s<s>.dcd`` + a shared ``traj.psf``). A per-residue ``<outdir>/dwell_times.dat``
    records the codon and the three sampled dwell times (seconds, ns, MD steps).
    Stitch the per-stage trajectories into one VMD movie with ``cosmo-csp-movie -o
    <outdir> --ribosome <ribosome>.pdb``. Set ``COSMO_CSP_VERBOSE=1`` for the full
    per-stage banners. See :doc:`continuous_synthesis` for the output tree.
