# TupleToolSemiLeptonic
This project provides several additional classes for the LHCb
`Analysis/Phys/DecayTreeTuple` package for lepton flavor universality violation
analysis in the semileptonic channel.

All rights reserved.


## Coding style
All `.cpp` and `.h` source files have been reformatted by:
```
clang-format -i -style=file *.cpp *.h
```
for visual consistency.

The formating template is located under `.clang-format` in the root directory
of the project.


## Acknowledgment
* `TupleToolApplyIsolationVetoDst, TupleToolTagDiscardDstMu` are used in LHCb
  run-1 semileptonic analysis. We obtained these from Phoebe Hamilton.
* `TupleToolApplyIsolation, TupleToolTauMuDiscrVars` are obtained from Patrick
  Owen. These can be found at [1].


[1]: https://gitlab.cern.ch/lhcb-slb/B02DplusTauNu/tree/master/tuple_production/tuple_tools_src
