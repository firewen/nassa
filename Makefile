FC=ifort

target1=nassa

src=$(wildcard *.f90)
obj=$(patsubst %.f90, %.o, $(src))

subs=raytracing.o m_npy.o load_st.o get_par.o  spheredist.o fitness.o estimate.o NA.o
mod=raytracing.mod  m_npy.mod
Flags=-qopenmp

all : $(target1) 

$(target1) : $(subs) nassa.o
	$(FC) $^ $(Flags) -o $@

%.o : %.f90 
	$(FC) $^ $(Flags) -c 

%.mod : %.f90
	$(FC) $^ -c

.PHONY : clean
clean :
	rm *.o $(target1) $(target2)
