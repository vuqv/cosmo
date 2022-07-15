Slab simulation protocol to scale box dimension.
Run with Monte Carlo algorithm to scale the box.
The general procedures are:

1) we don't know exactly what is the reasonable box dimention for the system, so initially we just set a big box that
   contains whole system.
2) run NPT simulation to scale all dimensions of box to the reasonable.
3) When the system reach target pressure, we keep, e.g xy dimension is as the final results of NPT simulation
   and extend the z-dimention such that the density of system is comparable as what we want to study.
   Turn off the NPT and switch to NVT ensemble only to study the density. 