#!/bin/bash

cd dDEMENTpy/input
rm climate.csv
mv base.csv climate.csv
cd ../..

cd dDEMENTpy1/input
rm climate.csv
mv scenario_2012.csv climate.csv
cd ../..

cd dDEMENTpy2/input
rm climate.csv
mv scenario_2013.csv climate.csv
cd ../..

cd dDEMENTpy3/input
rm climate.csv
mv scenario_1213.csv climate.csv
cd ../..