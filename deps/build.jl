println("Building Penelope shared library from Fortran")

penelope_location="/Users/fowlerj/Software/penelope/fsource"
f77="gfortran"

run(`$(f77) -I$(penelope_location) -g -c -o wrappers.o wrappers.f`)
run(`$(f77) -shared -g -fPIC -o penelope.so wrappers.o`)
