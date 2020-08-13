# Makes and runs the computation of the gravitational redshift and time of arrival
# -----------------------------------------------
# make_ujti.sh --> (( ztoa_make_run.sh )) --> lchotspot_get.py --> lchotspot_plot.py
#                                         --> thspectrum_get_plot.py
#                                         --> surf_*_plot.py

# Usage: bash ztoa_make_run.sh

# Note that the data and info are stored by default in data/zt and data/zt/meta

# Metric name
metricname=frutos

# Parameter transformations (see TR in ztoa-metric.f90)
# {frutos}: 1; {mink}:2;  {htco, ht, gb}: 3
metriccat=1

# Coordinates: Boyer-Lindquist or spherical (in numeric/coord.f90)
coord="BL_COORD"

outmetricname=$metricname

# inclination angle
angle=0
# NS configuration name
objcode="SHFT"

# resolution ID for data id. To change resolution, go to ztoa-metric.f90
resid="161"

gfortran ztoa-metric.f90 -cpp \
  -DMETRIC=$metricname -DOBJ=$objcode -DTR=$metriccat -DUSECOORD=$coord \
  -O0 -Ibin/so bin/so/geq-sphax-base.o\
  bin/so/$metricname.o\
  bin/so/nsl-solver.o\
  bin/so/coord.o\
  -o bin/ztoa-$resid-${objcode}${outmetricname}${angle}


time bin/ztoa-${resid}-${objcode}${outmetricname}${angle} << eof
${objcode}-$resid-${outmetricname}-i${angle}
${angle}
eof
