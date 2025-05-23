<!-- <p align="center"> <font size="6"> <b> DEMENTpy </b> </font> </p> -->

<!-- ![alt text](documentation/animations/bacteria.gif "Bacterial Taxon Dynamics"){ width=30% } ![alt text](documentation/animations/fungi.gif "Fungal Taxon Dynamics"){ width=30% } [alt text](documentation/animations/cellulose.gif "Cellulose Dynamics"){ width=30% } -->

<!--
<p align="center">
<img src="documentation/animations/bacteria.gif" width="256" title="Bacterial Taxon Dynamics"> <img src="documentation/animations/fungi.gif" width="256" title="Fungal Taxon Dynamics"> <img src="documentation/animations/cellulose.gif" width="256" title="Cellulose Dynamics">
</p>
-->

<p align='center'> <img src="documentation/animations/DEMENTpy_animation.gif"> </p>

# DEMENTpy
## A trait- and individual-based spatially explicit soil microbial systems modelling framework

<!--
<span style="color: red;"> [**NOTE: still under active development without any formal release of any version !!!**; if interested, feel free to reach out to me via any media] </span>
-->

![GitHub repo size](https://img.shields.io/github/repo-size/bioatmosphere/DEMENTpy)
![GitHub contributors](https://img.shields.io/github/contributors/bioatmosphere/DEMENTpy)
![GitHub stars](https://img.shields.io/github/stars/bioatmosphere/DEMENTpy?style=social)
![GitHub forks](https://img.shields.io/github/forks/bioatmosphere/DEMENTpy?style=social)
<!--![Twitter Follow](https://img.shields.io/twitter/follow/bioatmo_sphere?style=social)-->

DEMENT (Decomposition Model of Enzymatic Traits) is spatially and mechanistically explicit in simulating a microbial system comprising a large number of hypothetical microbial taxa in terrestial environments. As indicated by the 'py' in its name, DEMENTpy is developed and programmed in Python, based on its predecesor [DEMENT](https://github.com/stevenallison/DEMENT) developed by Steven Allison (originally programmed in R).

<p align='center'> <img src="documentation/figures/DEMENTpy_conceptual_structure.jpg"> </p>

### Vision

DEMENTpy was developed with practicing **Real Open Science** in mind and is devoted to 100% OPEN SOURCE and longterm maintanence and development. We sincerely seek inputs from communities ranging from microbial ecology, systems biology, theoretical ecology, etc. Please read closely our statement on **policies and rules** of making contributions to DEMENTpy or applying it to your own research.

<p align='center'> <img src="documentation/figures/DEMENTpy_programming_structure.png"> </p>

### Structure and Process

This model is built upon its predecessor--DEMENT, an R-based framework initially developed by Steven Allison back in 2012. In addition to the programming language, a series of changes have been made with an overarching goal of making it more readily accessible to the research and teaching communities. This model simulates plant litter decomposition from enzymatic degradation of substrates to microbial processes encompassing uptake, metabolism, mortality, reproduction, and dispersal in a spatially explicit, mechanistically explicit fashion. For more detailed information about DEMENTpy, we refer readers/users temporarily to an [Appendix](https://github.com/bioatmosphere/microbiome-drought-legacy/tree/master/writing) of the 1st manuscript applying this model. An actual documentation is under conceiving.


### Run DEMENTpy

**Get the code**:
```shell
git clone https://github.com/bioatmosphere/DEMENTpy
```

**Directory structure**:

- src/: all source code

- input/: data required to drive the model

- output/: folder where the output object in .pickle will be residing

- **SPA/: Scripts of Processing and Analyzing (SPA) outputs**

**Run DEMENTpy**:

[UV](https://github.com/astral-sh/uv) can be used to manage dependencies in `pyproject.toml`. In order to create a virtual ennviroment with DEMENTpy dependencies, do the following: 

```
uv venv
uv sync --dev 
source .venv/bin/activate
```
if you are not running as developer then leave off the `--dev` switch.

A simple example of bash script, dementpy.sh, for running jobs on HPC is provided. For testing, a simple run can be carried out as follows:

```
python dementpy.py grassland output 20250402 scrubland
```

### Contribution

The space for contribution is HUGE and OPEN! For instance:

- New processes/algorithms are more than welcome. 
- Programming improvement is absolutely needed.
- Any bugs are possible.
- ...

Feel free to reach out, or directly fork, change, and create pull requests, or create new branches to contribute.

### License

[MIT LICENSE](https://github.com/bioatmosphere/DEMENTpy/blob/master/LICENSE)

### References

1. Wang, B., & Steven D. Allison. (2022) [Climate-driven legacies in simulated microbial communities alter litter decomposition rates](https://www.frontiersin.org/articles/10.3389/fevo.2022.841824). Frontiers in Ecology and Evolution 10,841824 [[GitHub Repo](https://github.com/bioatmosphere/microbiome-climate-gradient)]

2. Wang, B., & Steven D. Allison. (2021) [Drought Legacies Mediated by Trait Tradeoffs in Soil Microbiomes]( https://doi.org/10.1002/ecs2.3562). Ecosphere 12, e03562 [[GitHub Repo](https://github.com/bioatmosphere/microbiome-drought-legacy)]

3. Wang, B., & Allison, S. D. (2019).[Emergent properties of organic matter decomposition by soil enzymes](https://doi.org/10.1016/j.soilbio.2019.107522). Soil Biology and Biochemistry, 136, 107522 [[GitHub Repo](https://github.com/bioatmosphere/An_emergent_soil_enzyme_decomposition_model)]

4. Allison, S. D., & Goulden, M. L. (2017). [Consequences of drought tolerance traits for microbial decomposition in the DEMENT model](https://doi.org/10.1016/j.soilbio.2017.01.001). Soil Biology and Biochemistry, 107, 104-113.

5. Allison, S. D. (2012). [A trait‐based approach for modelling microbial litter decomposition](https://doi.org/10.1111/j.1461-0248.2012.01807.x). Ecology Letters, 15, 1058-1070. 
