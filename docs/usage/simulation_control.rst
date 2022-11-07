Simulation control options
========================================================= 

An example of how config file of simulation looks like.

.. code-block::

        [OPTIONS]
        md_steps = 30_000 # number of steps
        dt = 0.01 ; time step in ps
        nstxout = 1000 ; number of steps to write checkpoint = nstxout
        nstlog = 1000 ; number of steps to print log
        nstcomm = 100 ; frequency for center of mass motion removal
        ; select HPS model, available options: hps_kr, hps_urry, or hps_ss
        model = hps_urry

        ; control temperature coupling
        tcoupl = yes
        ref_t = 310 ; Kelvin- reference temperature
        tau_t = 0.01 ; ps^-1

        ;pressure coupling
        pcoupl = yes
        ref_p = 1
        frequency_p = 25

        ; Periodic boundary condition: if pcoupl is yes then pbc must be yes.
        pbc = yes
        ; if pbc=yes, then use box_dimension option to specify box_dimension = x or [x, y, z], unit of nanometer
        box_dimension = 30 ; [30, 30, 60]

        ; input
        protein_code = FUS_100chains
        pdb_file = FUS_100chains.pdb
        ; output
        checkpoint = FUS_100chains.chk
        ;Use GPU/CPU
        device = GPU
        ; If CPU is specified, then use ppn variable
        ppn = 4
        ;Restart simulation
        restart = no
        minimize = yes ;if not restart, then minimize will be loaded, otherwise, minimize=False

==============================================================================================================================

General information
++++++++++++++++++++++++
Simulation parameters are input from `.ini` file which is loaded by `ConfigParser` module in Python.
The section title ``[OPTIONS]`` is required, do not change section's name.

* Comment can be inline or in new line, start with `;` or `#`
* Keyword and value can be separated by `=` or `:`

Run control
++++++++++++

::

    md_steps:   (long int)
                (1) Maximum number of steps to integrate or minimize
    ------------------------------------------------------------------------------------
    dt:         (double)
                (0.01)[ps] Time step for integration
    ------------------------------------------------------------------------------------
    nstxout:    (int)
                (1) [step] number of steps that elapse between writing coordinates to output trajectory file,
                      the last coordinates are always written
    ------------------------------------------------------------------------------------
    nstlog:     (int)
                (1) number of steps that elapse between writing energies to the log file
    ------------------------------------------------------------------------------------
    nstcomm:    (int)
                (100) frequency for center of mass motion removal

Model parameter
+++++++++++++++
There are three models supported now: `hps_kr`, `hps_urry` and `hps_ss`.

`hps_kr` has parameters for a wide range of residues, i.e RNA, phosphorylation residues ... but this model is less accurate

::

    model:      (string)
                hps_kr: Kapcha-Rossy hydropathy scale, parameterize from OPLS-AA forcefield

                hps_urry (default): Urry hydropathy scale, parameterize from experiment

                hps_ss: hps_urry with bonded potential (angle and torsion)


Temperature coupling
+++++++++++++++++++++

::

    tcoupl:     (bool)
                yes (default) : The only available option for now, we don't care about NVE ensemble.
    ------------------------------------------------------------------------------------
    ref_t:      (double)
                (300) [K] : Reference temperature in unit of Kelvin
    ------------------------------------------------------------------------------------
    tau_t:      (double)
                [ps^-1] : The friction coefficient which couples the system to the heat bath (in inverse picoseconds)

Pressure coupling
+++++++++++++++++++

::

    pcoupl      (bool)
                yes : Using pressure coupling

                no (default) : Run on NVT ensemble only
    ------------------------------------------------------------------------------------
    ref_p       (double)
                 (1) [bar] The default pressure acting on the system.
    ------------------------------------------------------------------------------------
    frequency_p (int)
                (25) [steps] the frequency at which Monte Carlo pressure changes should be attempted

Periodic boundary condition:
+++++++++++++++++++++++++++++
if pcoupl is yes then pbc must be yes.

::

    pbc         (bool)
                yes : Using periodic boundary condition.
                        If this option is chosen, then it will affect to non-bonded forces in the system,
                        and the coordinate writen in PDB and DCD file as well. No worries since I have handled these.

                no (default) : Without periodic boundary condition.
    ------------------------------------------------------------------------------------
    box_dimension   (float or list of float)
                [nm] An example of box dimension:
                If you want a cubic box of 30x30x30 nm^3, put: 30 or [30, 30, 30]
                If you want a rectangular box? Put:  [30, 30, 60]

File input/output
+++++++++++++++++++

::

    protein_code    (string)
                    String for output prefix, i.e {protein_code}.dcd, {protein_code}.log
    ------------------------------------------------------------------------------------
    pdb_file        (string)
                    [.pdb, .cif] Input structure for loading topology and initial coordinate
    ------------------------------------------------------------------------------------
    checkpoint      (string)
                    [.chk] Checkpoint file name, here I ask you to provide it explicitly since
                            because checkpoint can be used to load state or save state.
                            in case if you restart simulation with different name, you have to provide it.

Simulation platform
+++++++++++++++++++++
Simulation can be run on CPU with number of threads is control by `ppn` or using GPU.
If `device=CPU` then ppn need to be specify, otherwise simulation will run on 1 core

::

    device          (string)
                    GPU : Use gpu to run simulation

                    CPU (default) : use cpu to run simulation, if you specify cpu, you should modify ppn option, it control
                            how many cores will be used to run simulation, if not, default is 1.
    ------------------------------------------------------------------------------------
    ppn             (int)
                    (1) [threads] Number of threads used to run simulation on CPU. When using GPU,
                                performance is boosted a lot so ppn in that case is set to 1.

Restart simulation
++++++++++++++++++++

::

    restart         (bool)
                    yes : restart simulation from checkpoint file. This can be True, 1 or whatever are not (FALSE)
                            in python condition. If this option is selected, minimize will be force to False.

                    no (default) : Run simulation from beginning, if this option is selected, you can choose if you want to minimize your
                        system before running simulation.
    ------------------------------------------------------------------------------------

    minimize        (bool)
                    yes (default) : perform energy minimization before run molecular dynamics.

                    no : Not running energy minimization. This is default option when restart option is set to yes.