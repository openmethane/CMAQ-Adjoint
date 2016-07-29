# make file to build bldmake program

 EXEC = bldmake

# CHANGE: specify local Fortran and C compilers here

### intel
# FC = /usr/local/apps/intel/composer_xe_2013.2.146/bin/intel64/ifort
 FC = ifort
# FC = /usr/local/intel/ictce/3.2.2.013/fc/bin/intel64/ifort
#F_FLAGS = -O2 -fixed -extend_source -WB 
 F_FLAGS = -g -fixed -extend_source -WB 

# CC = /usr/local/intel/ictce/3.2.2.013/fc/bin/intel64/icc
 CC = icc

### PGI
#FC = /usr/local/pgi-11.9/linux86-64/11.9/bin/pgf90
#F_FLAGS = -Bstatic -Mextend -Mbounds -Mfixed

### GFortran
#FC = /usr/local/gcc-4.6/bin/gfortran
#F_FLAGS = -ffixed-form -ffixed-line-length-132

 OBJS= cfg_module.o parser.o utils.o bldmake.o

.SUFFIXES: .f

$(EXEC): $(OBJS)
	$(FC) $(OBJS) -o $@

.f.o:
	$(FC) -c $(F_FLAGS) $<

clean:
	rm -f *.o *.mod bldmake
# make file to build bldmake program

