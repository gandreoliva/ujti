# Compiles the ray tracing Fortran binary
# --------------------------------------
# make_ujti.sh --> (( raytr_make.sh )) --> raytr_run.sh
#
# This script adds the necessary parameter transformations for the general ray tracing
# Fortran source code in raytr-metric.f90
# Alternative: make the substitutions directly at the metric.mac files, but it may
# be slower! We only compute with RNS the values of m, a, q as defined in the HT metric,
# so, other metrics that use other definitions of q or parameters related to q
# should specify a transformation.
# If no transformation is found for a metric, none is inserted (a warning is displayed)


# Usage: bash raytr_make.sh metric_name

if [ -z $1 ]; then echo "$doc"; echo "**Error: No metric code to compile!"; exit 1; fi


paramtransf_tempfile="gensrc/raytr-paramtransf.f90"
metricname_tempfile="gensrc/raytr-metricname.f90"

# select metric module to be used
printf "use $1" > $metricname_tempfile

# Fortran code, parameter transformations
frutos_paramtransf="
q = -q
"
gb_paramtransf="
"
ht_paramtransf="
"
htco_paramtransf="
"
qm_paramtransf="
sigma = sqrt(m**2 - a**2)

q = -q/((sigma**2/m**2)**(3/2d0)*m**3)*15/2d0

if (a /= 0) then
  alpha = (sigma - m)/a
else
  alpha = 0
end if
"
mn2_paramtransf="
kk = sqrt(m**2-a**2)
if (a /= 0) then
  alpha = (kk - m)/a
else
  alpha = 0
end if

alpha2 = q/kk**3
"
# add new metrics here



# metric parameter transformation selector
# prints the transformations to a Fortran source file that is later included
# in raytr-metric.f90
case $1 in
"frutos")
  printf "$frutos_paramtransf" > $paramtransf_tempfile
  ;;
"ht")
  printf "$ht_paramtransf" > $paramtransf_tempfile
  ;;
"htco")
  printf "$htco_paramtransf" > $paramtransf_tempfile
  ;;
"qm")
  printf "$qm_paramtransf" > $paramtransf_tempfile
  ;;
"gb")
  printf "$gb_paramtransf" > $paramtransf_tempfile
  ;;
"mn2")
    printf "$mn2_paramtransf" > $paramtransf_tempfile
  ;;
# add new metrics here
*)
  printf "" > $paramtransf_tempfile
  echo " (!) Warning: no parameter transformation defined"
  ;;
esac

# compilation
gfortran raytr-metric.f90 -O0 -Ibin/so/ bin/so/geq-sphax-base.o\
  bin/so/$1.o\
  bin/so/indiv-ngeod-solver.o\
  bin/so/coord.o\
  -o bin/raytr-$1
