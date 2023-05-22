# The High-resolution Intermediate Complexity Atmospheric Research Model (HICAR)


[![DOI](https://zenodo.org/badge/638935780.svg)](https://zenodo.org/badge/latestdoi/638935780)

HICAR is a variant of the Intermediate Complexity Atmospheric Research (ICAR) model developed for sub-kilometer resolutions. The model is developed for downscaling of kilometer-scale NWP model output to resolutions used for land-surface simulations. HICAR features physics parameterizations shared by traditional weather models such as WRF, but with massively simplified dynamics which enable for run times nearly 600x faster than WRF.

More information about the model features enabling this can be found in Reynolds et al., 2023.

#### Compilation Requirements
While being fast to run compared to traditional weather models, HICAR has still been developed and intended for use on High Performance Computing (HPC) machines, although it can also be run on a local machine. HICAR thus uses a few package requirements which are common to HPC environments. They are:

- Parallel NetCDF4
- FFTW
- PETSc

HICAR currently supports either the GNU Fortran or Cray compilers. The Cray compiler offers a better adoption of the Coarray-Fortran standard and faster optimization options, and is thus significantly faster than the GNU compiler. However, debugging is sometimes easier with the GNU compiler.

#### Static data requirements

HICAR uses a domain file which defines land-surface variables and some terrain-descriptors. To generate a HICAR domain file using an existing DEM and land use data, the  python script gen_HICAR_dom.py, located in helpers/ can be used. See below for more details on how to run this script.

Example static data for running a 1-day simulation can be found under [HICAR-model/Test-data](https://github.com/HICAR-Model/Test-Data)

#### Forcing data requirements

HICAR requires at least the following 5 fields to run, all of which are contained within one netCDF file:
- U winds
- V winds
- Pressure
- Temperature (potential or normal temperature)
- Humidity (mixing ratio or specific humidity)

When using the variational wind solver, providing W winds from forcing data can lead to better estimates of wind speeds.

To perform nested runs, HICAR can be forced with the output from a previous HICAR simulation. Thus HICAR also supports the forcing of all hydrometeors and all of their moments as according to the microphysics scheme chosen. It is recommended to specify these forcing variables when forcing HICAR with output from coarser resolution HICAR runs.

HICAR reads in forcing data from a forcing file list supplied to the model in the namelist. A shell script for generating a forcing file list from a given directory is found within helper/filelist_script.sh

Example forcing data for running a 1-day simulation can be found under: [HICAR-model/Test-data](https://github.com/HICAR-Model/Test-Data)

#### Namelist

An example of a test namelist can be found under run/HICAR_Test_Case.nml. This namelist is the same which can be run with the test case provided in [HICAR-model/Test-data](https://github.com/HICAR-Model/Test-Data). The complete namelist options found in run/namelists/complete_hicar_options.nml show all possible name list options, with comments describing their function and use. Full documentation on namelist options to run with the model are in development…

#### Supplementary data

Supplementary data for running HICAR, for example look-up tables needed for the ISHMAEL microphysics scheme, can be found in the following repo: https://github.com/NCAR/icar_supporting_files

#### Developing
If you plan to make any major additions to HICAR or ICAR, please get in touch, for minor changes feel free to just submit a pull request. The current workflow is to make changes and pull requests to the `develop` branch.

If your code changes are core to the model and would benefit both HICAR and ICAR, please make them over at the ICAR model repo: https://github.com/NCAR/icar.

For an outline of the basic code structure see the [ICAR code overview](docs/icar_code_overview.md)

For reference working with the model code and git, see the [ICAR and Git workflow](docs/howto/icar_and_git_howto.md).


####Generating Static Data

HICAR relies on pre-computed static data to speed up some of it’s online calculations. To generate a HICAR domain file, an existing netCDF file with lat, lon, DEM, and landuse categories is needed. The lat and lon variables must be named **lat** and **lon**, and the terrain variable must be named **topo**. Additionally, a larger extent DEM of the same resolution is needed to generate parameters for terrain-shading of radiation. I.e., if you have a 50m resolution domain, a larger DEM with an extent ~20km beyond the boundaries of the target domain is also needed.

Once you have these two netCDF files, you can use a python script to generate the rest of the variables used by HICAR.

First, install the conda environment located in the HICAR_dom.yml file found in helpers/

```bash
conda env create -f HICAR_dom.yml
```

Once this environment is installed, activate it with:
```bash
conda activate HICAR_dom
```

Now you will need to install HORAYZON, a python package to efficiently calculate the horizon line matrix and sky view factor (Steger et al., 2022). To do so, type:

```bash
git clone https://github.com/ChristianSteger/HORAYZON.git
cd HORAYZON
python -m pip install .
```

Your python environment should now be complete. To generate the domain file, open gen_HICAR_dom.py and edit the paths for the target domain, radiation domain, and output domain. HICAR_Domain.py and ProjHelpers.py, both contained in the helpers/ directory, must be in the same directory as gen_HICAR_dom.py.

Now run:

```bash
python gen_HICAR_dom.py
```

#### Reference

Reynolds, D. S., Gutmann, E., Kruyt, B., Haugeneder, M., Jonas, T., Gerber, F., Lehning, M., and Mott, R.: The High-resolution Intermediate Complexity Atmospheric Research (HICAR v1.0) Model Enables Fast Dynamic Downscaling to the Hectometer Scale, Geosci. Model Dev. Discuss. [preprint], https://doi.org/10.5194/gmd-2023-16, in review, 2023. 

Gutmann, E. D., I. Barstad, M. P. Clark, J. R. Arnold, and R. M. Rasmussen (2016), *The Intermediate Complexity Atmospheric Research Model*, J. Hydrometeor, doi:[10.1175/JHM-D-15-0155.1](http://dx.doi.org/10.1175/JHM-D-15-0155.1).

Steger, C. R., Steger, B. and Schär, C. (2022): HORAYZON v1.2: an efficient and flexible ray-tracing algorithm to compute horizon and sky view factor, Geosci. Model Dev., 15, 6817–6840, https://doi.org/10.5194/gmd-15-6817-2022
