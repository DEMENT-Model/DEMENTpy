 <!-- <p align="center"> <font size="6"> <b> DEMENTpy </b> </font> </p> -->

![alt text](documentation/animations/bacteria.gif "Bacterial Taxon Dynamics"){ width=30% } ![alt text](documentation/animations/fungi.gif "Fungal Taxon Dynamics"){ width=30% } [alt text](documentation/animations/cellulose.gif "Cellulose Dynamics"){ width=30% }

<p align="center">
<img src="documentation/animations/bacteria.gif" width="256" title="Bacterial Taxon Dynamics"> <img src="documentation/animations/fungi.gif" width="256" title="Fungal Taxon Dynamics"> <img src="documentation/animations/cellulose.gif" width="256" title="Cellulose Dynamics">
</p>

# DEMENTpy
## A trait- and individual-based spatially explicit soil microbial systems modelling framework

[**NOTE: still under active development !!!**; if interested, feel free to reach out to me via any media]

![GitHub repo size](https://img.shields.io/github/repo-size/bioatmosphere/DEMENTpy)
![GitHub contributors](https://img.shields.io/github/contributors/bioatmosphere/DEMENTpy)
![GitHub stars](https://img.shields.io/github/stars/bioatmosphere/DEMENTpy?style=social)
![GitHub forks](https://img.shields.io/github/forks/bioatmosphere/DEMENTpy?style=social)
![Twitter Follow](https://img.shields.io/twitter/follow/bioatmo_sphere?style=social)

This model is spatially and mechanistically explicit in simulating a microbial system comprised of a large number of microbial taxa in terrestial environments. As indicated by the 'py' in its name, DEMENTpy is developed and programmed in Python, based on its predecesor DEMENT which is R-based.

**Vision**

DEMENTpy, as the first step of practicing real **Open Science** myself, is devoted to longterm maintanence and development with inputs not only from ourselves but also, hopefully, from the communities that can be as broad as microbial ecology, systems biology, theoretical ecology, etc. Just because of community inputs, please read closely our statement on **policies and rules** of making contributions to DEMENTpy or applying it to your own research.

**Structure and Process**

This model is built upon its predecessor--DEMENT, an R-based framework initially developed by Allison back in 2012. Except for the programming language change, a series of changes have been made with an overarching goal of making it more readily accessible to the research and teaching communities as broad as microbial ecology, theoretical ecology, and ecosystem ecology, as well as biology. This model simulates processes ranging from degradation of substrates through microbial processes encompassing uptake, metabolism, mortality, reproduction, and dispersal in a spatially explicit, mechanistically explicit fashion. Here is the underlying conceptual structure of DEMENTpy. For more detailed information about DEMENTpy, we refer readers/users to the **Documentation**(in pdf) archived in the documentation/ folder. 


**Running DEMENTpy**

Get the code:
```shell
git clone https://github.com/bioatmosphere/DEMENTpy
```

File structure and run DEMENTpy

- src/: all source code

- input/: data required to drive the model

- output/: folder where the output object in .pickle will be residing

Run DEMENTpy on HPC

- a simple example of bash script, dementpy.sh submitting jobs to HPC is provided.

**Contributing Guide**

Please follow these rules if you want to contribute to this open source project:

It is coming soon!

**License**

MIT LICENSE
