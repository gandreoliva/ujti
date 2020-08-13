# Runs the ray tracing Fortran binary
# --------------------------------------
# make_ujti.sh --> raytr_make.sh --> (( raytr_run.sh ))
#
# The parameter values are specified here.

# Usage: bash raytr_run.sh metric_name


metricname=$1

m='0.2'
a='0.15'
q='0.01'
ar='0.6' # ratio of polar to equatorial radii.
         #if -1 < ar < 0, it's computed automatically assuming solid body rotation

view='1' # 1 is midplane, 2 is vertical plane

tol="1d-10"
dlmin="1d-8"

# (ray-tracing) data id (directory)
dataid="raytr_test1"

# create directory (notice root directory)
mkdir -p data/raytr/$dataid/meta


echo $dataid

# run the ray tracing binary
time bin/raytr-$metricname << eof
${dataid}
${metricname}
${m}
${a}
${q}
${ar}
${view}
${tol}
${dlmin}
eof
