# Ujti

```
--- -- -- -- ---._      UJTI - general geodesic ray-tracer             
                o `.                                   o     _         
                    \          Cinespa UCR 2019  *  ( ) GANDREOLIVA.org
=======================================================================
```

## Overview
Ujti (pronounced in English as 'OOhh-tee', from Nawat language 'way') is a software package that uses Maxima, Fortran and Python, and calculates null geodesics (as well as other applications) for arbitrary axisymmetric metrics.

The program has three main starting points:

1. Symbolic calculations (Maxima)
    *    Generation of the geodesic equations from any (analytically defined) axisymmetric metric
2. Numerical calculations (Fortran)
    *    Ray-tracing: traces null geodesics
    *    General relativistic frequency shift and time of arrival (provides data for the post-processing applications)
3. Post-processing (Python)
    *    Thermal spectrum of the neutron star
    *    Light curves from a circular hotspot
    *    General relativistic redshift from the surface of the neutron star
    *    Shape of the neutron star as observer by a distant observer

## Setup instructions
0.    Check that the following (empty) directories exist: `bin/so/`, `data/zt/meta`, `gensrc/`
1.    Symbolic calculations: (metric->geod. eqns.) just run `maxima -b [name of the metric file].mac`. This generates the Fortran source code in `gensrc/`
2.    Solvers: `./make_ujti.sh -num [name of the solver file].f90`. The available solver files are in `numeric/`, and they provide modules and .o objects that are reused in all applications (no binary executables are generated). The auxiliary module in `coord.f90` must also be compiled this way. Important note: the file `symbolic/geq-sphax-base.f90` can be also compiled this way, and must be compiled before the next step.
3.    Fortran module for the specific metric: `./make_ujti.sh -symb [code of the metric]`. This provides the metric modules for the numerical calculations.
4.    Numerical applications, in general: they should be compiled by including the .o files for the base axisymmetric modules, the specific metric, the solver and the coordinate transformations module.


## Quick guides
### Light scattering (ray tracing in the midplane)
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

### Light curves for hotspots
After the Ujti setup is done, do the following:

Compile and run the ztoa (redshift+time of arrival) binary. Check the configuration in the file. Notice the dataid in the file. Repeat for different metrics, neutron star parameters or inclinations.

`./ztoa_make_run.sh`

Compute the light curve for a specific hotspot, specify configuration in command line args. Repeat for different hotspot configurations (location, size, shape). In this example, a circular hotspot of colatitude 45° and ang. radius 10°

`python lchotspot_get.py SHFT-161-frutos-circ_45_10-i0 SHFT-161-frutos-i0 45 10 0`

Plot the light curve. Notice the dataid in the file.

`python lchotspot_plot.py`


## Version
This is Ujti version 3.0. There was a previous version Ujti 2.0, still available from https://bitbucket.org/gandreoliva/ujti2 , that included the calculation of gravitational lenses. Ujti 2.0 was used in the article https://doi.org/10.15517/rmta.v22i2.20723 .

## License

---

Please acknowledge any use of this software by choosing any of the following means:

  *    Citation to Oliva-Mercado and Frutos-Alfaro 2020 (MNRAS, subm.) arXiv:2006.05948
  *    Acknowledgment including the URL of availability of the code (either http://cinespa.ucr.ac.cr or http://gandreoliva.org)

If you find this software useful, I encourage you to drop me an email.

---

BSD License

Copyright (c) 2019 Andre Oliva.
All rights reserved.


Redistribution and use in source and binary forms are permitted
provided that the above copyright notice and this paragraph are
duplicated in all such forms and that any documentation, advertising
materials, and other materials related to such distribution and use
acknowledge that the software was developed by Andre Oliva. The
name of the Andre Oliva may not be used to endorse or promote
products derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED “AS IS” AND WITHOUT ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
