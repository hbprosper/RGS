# RGS
## Introduction
The Random Grid Search (RGS) algorithm is a simple, but surprisingly effective, way to find rectangular cuts. Developed by Harrison Prosper, Chip Stewart, Pushpa Bhat and generalized by Sezen Sekmen to include *two-sided* and *staircase* cuts. A two-sided cut is a cut of the for (x1 < x < x2), while a staircase cut is the OR of two or more one-sided cuts.

## Installation
This package depends on the package [Root](https://root.cern.ch/downloading-root) from CERN. To install Root follow the instructions at the [Root](https://root.cern.ch/downloading-root) website. Then do
```
git clone https://github.com/hbprosper/RGS.git
cd RGS
make
source setup.sh
```
The setup need be done only once per terminal session. (Note: the bash setup can be "sourced" from any directory, but currently the non-bash setup must be sourced from the RGS directory.) 

## Examples
There are two examples in the *examples* directory of RGS. These examples require Root data files, which can get obtained using the commands
```
cd examples/data
wget http://www.hep.fsu.edu/~harry/RGS/data/Higgs.tar.gz 
tar zxvf Higgs.tar.gz
```
or downloaded from the website via a web browser. Use a similar procedure for the SUSY data files _Susy.tar.gz_.

### Higgs
This example illustrates three RGS optimizations, HO1, HO2, and HO3, designed to enhance the ratio of Higgs vector boson fusion (VBF) events to Higgs gluon gluon fusion (ggF) events and di-Z boson events. Each optimization can be run by executing the _train.py_ program followed by _analysis.py_ program. For example,
HO1 can be run as follows
```
cd examples/Higgs/HO1
./train.py
```
which will run RGS and store its results in a file called *HO1.root*. To analyze the results of RGS do
```
./analysis.py
```
which will read the results from *HO1.root* and write the results to *r_HO1.txt* and also produce a couple of plots.

### SUSY
This example illustrates three optimizations, SO1, SO2, and SO3 that use staircase cuts to improve the search for SUSY events. Switch to the *SUSY* directory and proceed as in the Higgs example. 
