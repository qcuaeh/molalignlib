gfortran -fPIC -c messages.f90
gfortran -fPIC -c lapack.f90
gfortran -fPIC -c math.f90
gfortran -fPIC -c sorting.f90
gfortran -fPIC -c hungarian.f90
gfortran -fPIC -c strutils.f90
gfortran -fPIC -c maputils.f90
gfortran -fPIC -c printing.f90
gfortran -fPIC -c chemistry.f90
gfortran -fPIC -c options.f90
gfortran -fPIC -c random.f90
gfortran -fPIC -c decoding.f90
gfortran -fPIC -c assignment.f90
gfortran -fPIC -c assortment.f90
gfortran -fPIC -c translation.f90
gfortran -fPIC -c rotation.f90
gfortran -fPIC -c alignment.f90
gfortran -fPIC -c readwrite.f90
gfortran -fPIC -c biasing.f90
gfortran -fPIC -c remapping.f90
gfortran -fPIC -c superposition.f90
f2py3.4 -m ralign -h ralign.pyf options.f90 superposition.f90
f2py3.4 -L/usr/lib64/atlas -llapack -c ralign.pyf *.o
