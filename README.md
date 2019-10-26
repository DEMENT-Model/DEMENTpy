 

<p align="center"> <b>DEMENTpy</b> </p>

----

**A trait- and individual-based spatially explicit soil microbial systems modelling framework -- DEMENTpy v1.0**

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

File structure:

./src: all source code

./input: data required to drive the model

./output: folder where the output object in .pickle will be residing

dementpy.sh: bash script for submitting job to HPC.

**References**
- Allison 2012 Ecology Letters
- Wang and Allison 2019 Soil Biology and Biochemistry
- Read this blog post on how to make this project open-sourced: https://dev.to/yvonnickfrin/preparing-your-project-being-open-sourced-5bdp?utm_source=digest_mailer&utm_medium=email&utm_campaign=digest_email
