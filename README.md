# RGS
## Introduction
The Random Grid Search (RGS) algorithm is a simple, but surprisingly effective, way to find rectangular cuts. Developed by Harrison Prosper, Chip Stewart, Pushpa Bhat and generalized by Sezen Sekmen to include *box* and *ladder* cuts. A box cut is the AND of two-sided cuts, for example (x1 < x < x2) AND (y1 < y < y2), while a ladder cut is the OR of two or more one-sided cuts.

## Installation
This package depends on the package [Root](https://root.cern.ch/downloading-root) from CERN and a small package [histutil](https://github.com/hbprosper/histutil) built on top of Root. To install Root follow the instructions at the [Root](https://root.cern.ch/downloading-root) website. The histutil package may be installed in any convenient place using 
```
git clone https://github.com/hbprosper/histutil.git
cd histutil
source setup.sh
cd -
```
Using setup.csh if you are using a non-bash shell. Then do
```
git clone https://github.com/hbprosper/RGS.git
cd RGS
make
source setup.sh
```
Both setups need be done only once per terminal session.

## Examples
There are two examples in the *examples* directory of RGS
### example1
This example illustrates a search for the best box cut. To run it do
```
cd examples/example1
python train.py
```
which will run RGS and store the results in a file called *example1.root* and also in a simple format in the text file *example1.txt*. To analyze the results of RGS do
```
python analysis.py
```
which will read the results from *example1.root*. The other program reads from the text file instead.

### example2
This example illustrates a search for the best ladder cut. Switch to the *example2* directory and proceed as in example1. 
