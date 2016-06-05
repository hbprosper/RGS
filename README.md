# RGS
## Introduction
The Random Grid Search (RGS) algorithm is a simple, but surprisingly effective, way to find rectangular cuts. Developed by yours truly, Chip Stewart, Pushpa Bhat and generalized by Sezen Sekmen to include *box* and *ladder* cuts. A box cut is the AND of two-sided cuts, for example (_x_1 < _x_ < _x_2) AND (_y_1 < _y_ < _y_2), while a ladder cut is the OR of two or more one-sided cuts.

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
