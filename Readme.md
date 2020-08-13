```
--- -- -- -- ---._      UJTI - general geodesic ray-tracer             
                o `.                                   o     _         
                    \          Cinespa UCR 2019  *  ( ) GANDREOLIVA.org
=======================================================================
```

# Overview
Ujti (pronounced in English as 'OOhh-tee') is a software package that uses Maxima, Fortran and Python, and calculates null geodesics (as well as other applications) for arbitrary axisymmetric metrics.

The program has three main starting points:

1. Symbolic calculations (Maxima)
  * Generation of general geodesic equations for an axisymmetric metric (auxiliary program)
  * Potential definition: takes a metric and generates Fortran code to evaluate it in 'numerical calculations'
2. Numerical calculations (Fortran)
  * Ray-tracing: traces null geodesics
  * Gravitational redshift, which also provides data for surface_colors and polar_cap
3. Post-processing (Python)
  * Plotting
  * Filter an arbitrary polar cap (or spot)
  * Plot surface points of the neutron star

# Compilation
0. Check that the following (empty) directories exist: `bin/so/`, `data/zt/meta`, `gensrc/`
1. Symbolic calculations: (metric->geod. eqns.) just run `maxima -b [name of the metric file].mac`. This generates the Fortran source code in `gensrc/`
2. Solvers: `./make_ujti.sh -num [name of the solver file].f90`. The available solver files are in `numeric/`, and they provide modules and .o objects that are reused in all applications (no binary executables are generated). The auxiliary module in `coord.f90` must also be compiled this way. Important note: the file `symbolic/geq-sphax-base.f90` can be also compiled this way, and must be compiled before the next step.
3. Fortran module for the specific metric: `./make_ujti.sh -symb [code of the metric]`. This provides the metric modules for the numerical calculations.
4. Numerical applications, in general: they should be compiled by including the .o files for the base axisymmetric modules, the specific metric, the solver and the coordinate transformations module.


# Quick guides
## Light scattering (ray tracing in the midplane)
In this example, we want to plot the null geodesics for the Fru16 metric.

Setup Ujti from scratch:
Generate the Fortran source code for the metric (necessary only once)
`maxima -b metric-frutos.mac`
Compile the necessary solvers and auxiliary modules (necessary only once)
`./make_ujti.sh -num numeric/coord.f90`
`./make_ujti.sh -num numeric/indiv-ngeod-solver.f90`
`./make_ujti.sh -num numeric/nsl-solver.f90`
`./make_ujti.sh -num symbolic/geq-sphax-base.f90`
Compile the metric (necessary only once)
`./make_ujti.sh -symb frutos`

Ray-tracing steps:
Compile the ray tracing binary for the Fru16 metric (necessary only once)
`./raytr_make.sh frutos`
Run the ray tracing binary. Check the configuration in the file. Notice the dataid in the file.
`./raytr_run.sh frutos`
Plot the geodesics. Notice the dataid in the file.
`python raytr_plot.py`

## Light curves for hotspots
After the Ujti setup is done, do the following:
Compile and run the ztoa (redshift+time of arrival) binary. Check the configuration in the file. Notice the dataid in the file. Repeat for different metrics, neutron star parameters or inclinations.
`./ztoa_make_run.sh`
Compute the light curve for a specific hotspot, specify configuration in command line args. Repeat for different hotspot configurations (location, size, shape). In this example, a circular hotspot of colatitude 45° and ang. radius 10°
`python lchotspot_get.py SHFT-161-frutos-circ_45_10-i0 SHFT-161-frutos-i0 45 10 0`
Plot the light curve. Notice the dataid in the file.
`python lchotspot_plot.py`


> Last update: 13.08.2020
