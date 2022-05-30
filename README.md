# hpsOpenMM

 OpenMM codebase for IDP which uses HPS-Urry/HPS-KR force field

This is a repository contains the source code of hps-urry/kr model.

 ```hps: hydropathy scale. ```
 ```urry: Hydropathy according to Urry scale.```
```kr: Kapcha-Rossy scale```
## requirements:
- **OpenMM 7.x**  
- **Parmed** 

## How to use hpsOpenMM:  
- Download folder and place in target location, for example: `/home/qvuvan/work3/code/hpsOpenMM/`
- Add folder in Python path (in .bashrc file): `export PYTHONPATH=$PYTHONPATH:/home/qvuvan/work3/code/hpsOpenMM/`  
- The standard example can be found at ```example/standard_example```

- *Note that in cluster, when submit job, the environment may not load `.bashrc` so need to add export command in job file* 

## Bugs
- If you encounter any bugs, please report the issue to Quyen Vu (```vuqv.phys@gmail.com```)

## Acknowledgments

- We would like to thank the OpenMM team for developing an excellent and open source engine. 
- We would like to thank the SBMOpenMM's developer team for developing an excellent model and open sourcing the software. 
- We also would like to thank Sugita's group to make a clear tutorial to introduce many features of OpenMM

## Do you want to cite this work?
- We will make this source code open to everyone to use.
