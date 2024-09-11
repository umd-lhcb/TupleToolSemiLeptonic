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


## Local development

If you want to modify this `TupleTool` and run a local ntupling job with `lhcb-ntuples-gen` you will need to
set up a new `DaVinci` inside `docker` with all the dependencies we use in our reco script. The workflow would be
```shell
cd lhcb-ntuples-gen
make docker-dv ## Enter Docker
export DAVINCI_VERSION=v45r6
export PHYS_VERSION=v26r6
export TUPLETOOL_SL_VERSION=0.2.4
export TRACKER_ONLY_EMU_VERSION=0.2.2

# Setup a DaVinci dev
source /usr/local/bin/envset.sh
git config --global user.name "Physicist"
git config --global user.email "lhcb@physics.umd.edu"
lb-dev DaVinci/${DAVINCI_VERSION}
git clone https://github.com/umd-lhcb/TrackerOnlyEmu.git --branch ${TRACKER_ONLY_EMU_VERSION} --depth 1
cd DaVinciDev_${DAVINCI_VERSION}
git lb-use Phys
git lb-checkout Phys/${PHYS_VERSION} Phys/LoKiPhys
git lb-checkout Phys/${PHYS_VERSION} Phys/DaVinciTypes
git lb-checkout Phys/${PHYS_VERSION} Phys/RelatedInfoTools
git lb-use TupleToolSemiLeptonic https://github.com/umd-lhcb/TupleToolSemiLeptonic.git
git lb-checkout TupleToolSemiLeptonic/${TUPLETOOL_SL_VERSION} Phys/TupleToolSemiLeptonic
cp -r ../TrackerOnlyEmu/davinci/* .

## Now you can edit the source code in lhcb-ntuples-gen/DaVinciDev_v45r6/Phys/TupleToolSemiLeptonic/src

# Compile DaVinci
make configure && make

# Run ntupling job
cd ../run2-rdx
..//DaVinciDev_v45r6/run gaudirun.py conds/cond-std-2016.py ./reco_Dst_D0.py
```
If you `git annex get` the required `.DST` files, this will produce a sample `std.root` in the `run2-rdx` folder.

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
