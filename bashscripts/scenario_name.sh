#!/bin/bash

cd dementpy0/input
rm climate.csv
mv base.csv climate.csv
cd ../..

cd dementpy1/input
rm climate.csv
mv scenario_2012.csv climate.csv
cd ../..

cd dementpy2/input
rm climate.csv
mv scenario_2013.csv climate.csv
cd ../..