# uBooNE e+e- Trident

This repo stores an updated version of the TEG c++ neutrino trident event generator, shared by [Altmannshofer et al. ](https://arxiv.org/abs/1902.06765). The core code remains unchanged, and the executable now handles bnb gsimple formatted flux histograms, and can be used to calculate experiment and fixed energy cross sections as well as generating neutrino trident events. 

This repo will be updated as new tools are developed, looking at the kinematics of the trident process and making comparisons to other SM/BSM models with similar final states. Currently the focus will be on the coherent scattering to an $e^+e^-$ final state.

## Installation

To build TEG_v2, you will need an existing `ROOT` install. If ups is availble, this can be done by running 
```
setup root v6_28_12 -q e26:p3915:prof
```
in an SL7 container. Once `ROOT` is setup, you can use `g++` and the `Makefile` to compile the executable by running `make`. This will create an executable called TEV_v2, which you can call in this directory by using `./TEG_v2` or add it to your path to call anywhere.
