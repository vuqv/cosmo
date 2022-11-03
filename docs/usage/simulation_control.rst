Simulation control options
========================================================= 

This file contains an example of how config file of simulation looks like.

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

