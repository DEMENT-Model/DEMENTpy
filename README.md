 <---! <p align="center"> <font size="6"> <b> DEMENTpy </b> </font> </p> --->
 
# DEMENTpy

## A trait- and individual-based spatially explicit soil microbial systems modelling framework

![GitHub repo size](https://img.shields.io/github/repo-size/bioatmosphere/DEMENTpy)
![GitHub contributors](https://img.shields.io/github/contributors/bioatmosphere/DEMENTpy)
![GitHub stars](https://img.shields.io/github/stars/bioatmosphere/DEMENTpy?style=social)
![GitHub forks](https://img.shields.io/github/forks/bioatmosphere/DEMENTpy?style=social)
![Twitter Follow](https://img.shields.io/twitter/follow/bioatmo_sphere?style=social)

This model is spatially and mechanistically explicit in simulating a soil microbial system comprised of a large number of microbial taxa. As indicated by the py in its name, DEMENTpy is developed and programmed in Python, based on its predecesor DEMENT which is R-based.

**Vision**

DEMENTpy is devoted to longterm maintanence and development with continuous updates not only from ourselves but also, hopefully, from the communities that can be as broad as microbial ecology, systems biology, theoretical ecology, etc. Just because of community inputs, please read closely our statement on **policies and rules** of making contributions to DEMENTpy or applying it to your own research.

**Structure and Process**

This model is built upon its predecessor--DEMENT, an R-based framework initially developed by Allison (2012). Except for the programming language change, a series of changes have been made with an overarching goal of making it more readily accessible to the research and teaching community revolving around microbial ecology, theoretical ecology, and ecosystem ecology, as well as biology. The structure of this model and processes therein are detailed below:

1. Structure:



2. Processes

- 2.0 Initialization

- 2.1 Degrdation

- 2.2 Uptake

- 2.3 Metabolism

- 2.3 Mortality

- 2.4 Reproduction

**Running DEMENTpy**

File structure

./src: all source code

./input: data required to drive the model

./output: folder where the output object in .pickle will be residing

dementpy.sh: bash script for submitting job to HPC.

**Contributing Guide**

Please follow these rules if you want to contribute to this open source model:

**References**
- Allison 2012 Ecology Letters
- Wang and Allison 2019 Soil Biology and Biochemistry
