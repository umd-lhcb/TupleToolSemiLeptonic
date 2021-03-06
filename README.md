# TupleToolSemiLeptonic
This project provides several additional classes for the LHCb
`Analysis/Phys/DecayTreeTuple` package for some semileptonic analyses.

All rights reserved.


## Coding style
All `.cpp` and `.h` source files have been reformatted by:
```
clang-format -i -style=file *.cpp *.h
```
for visual consistency.

The formating template is located under `.clang-format` in the root directory
of the project.


## Installation on lxplus
Type these commands in your terminal:

```
lb-dev DaVinci/<version> && cd DaVinciDev_<version>
git lb-use TupleToolSemiLeptonic https://github.com/umd-lhcb/TupleToolSemiLeptonic.git
git lb-checkout TupleToolSemiLeptonic/master Phys/TupleToolSemiLeptonic
make configure && make
```


## Development
It is highly recommended to use a docker with `DaVinci` and `CLion` shipped. It
has smart auto-completion and detects errors on-the-fly. To use the docker:
```
docker run --rm -it -v $(pwd):/data \
    -v $HOME/.config/clion-java:/home/physicist/.java \
    -v $HOME/.config/clion:/home/physicist/.CLion2019.2 \
    -v $XAUTHORITY:/home/physicist/.Xauthority \
    -e DISPLAY -e UID=$(id -u) -e GID=$(id -g) \
    --net=host \
    umdlhcb/lhcb-stack-cc7:DaVinci-v45r4
```
Note that the `clion` folder stores IDE configuration, and `clion-java` folder
store the licensing ticket.


## Acknowledgment
This repository is a superset of [`SemileptonicCommonTools`](https://gitlab.cern.ch/lhcb-slb/SemileptonicCommonTools), with additional
tools coming from [`B02DplusTauNu`](https://gitlab.cern.ch/lhcb-slb/B02DplusTauNu/tree/master/tuple_production/tuple_tools_src) analysis.
We made some changes and ported them to `DaVinci/v45`, the recommended version
for all run 1 and 2 analyses.
