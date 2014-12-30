#Compilers
PCC = mpicc
FORTRAN = gfortran

#Flags C standard
FLAGS = -std=gnu99

#Flags for performance
OPTFLAGS = -O3 -march=nocona

#Libraries
LIBS = -lm -lpq

#Objects
SHAREDOBJS = rkmethods.o problems.o mathmethods.o riversys.o sort.o comm.o system.o processdata.o partition.o definetype.o misc.o \
	rainfall.o solvers.o io.o forcings.o compression.o date_manip.o asynch_interface.o modeloutputs.o data_types.o
LIBPYOBJS = asynch_interface_py.o
ASYNCHDISTOBJS = asynchdist.o

#Locations
LIBDIR = ./libs
OBJDIR = ./objects
WITHTAOOBJDIR = ./objects
FOBJDIR = ./objects
OBJDIR_LIB = ./objects
SHAREDOBJSLOC_LIB = $(addprefix $(OBJDIR_LIB)/,$(SHAREDOBJS))
SHAREDOBJSLOC_LIB_PY = $(addprefix $(OBJDIR_LIB)/,$(LIBPYOBJS))
SHAREDOBJSLOC = $(addprefix $(OBJDIR)/,$(SHAREDOBJS))
ASYNCHDISTOBJSLOC = $(addprefix $(OBJDIR)/,$(ASYNCHDISTOBJS))


#How to compile and link
$(OBJDIR)/%.o: %.c
	$(PCC) -c $*.c $(HEADLOC) $(FLAGS) $(OPTFLAGS) -o $(OBJDIR)/$*.o

$(WITHTAOOBJDIR)/%.o: %.c
	$(PCC) -c $*.c $(HEADLOC) $(FLAGS) $(OPTFLAGS) $(TAOFLAGS) -o $(WITHTAOOBJDIR)/$*.o

$(FOBJDIR)/%.o: %.f
	$(FORTRAN) -c $*.f $(HEADLOC) $(FLAGS) $(OPTFLAGS) $(PRINTFLAG) -o $(FOBJDIR)/$*.o

$(OBJDIR_LIB)/%.o: %.c
	$(PCC) -c -fPIC $*.c -I/usr/include/python2.4/ -L/usr/lib/python2.4/ $(HEADLOC) $(FLAGS) $(OPTFLAGS) -o $(OBJDIR)/$*.o

ASYNCHLIB: $(SHAREDOBJSLOC_LIB)
	$(PCC) $(SHAREDOBJSLOC_LIB) $(LIBS) $(FLAGS) $(OPTFLAGS) -shared -Wl,-soname,libasynch.so -o $(LIBDIR)/libasynch.so

ASYNCHLIB_PY: $(SHAREDOBJSLOC_LIB) $(SHAREDOBJSLOC_LIB_PY)
	$(PCC) $(SHAREDOBJSLOC_LIB) $(SHAREDOBJSLOC_LIB_PY) $(LIBS) $(FLAGS) $(OPTFLAGS) -lpython2.4 -shared -Wl,-soname,libasynch_py.so -o $(LIBDIR)/libasynch_py.so

ASYNCH: $(SHAREDOBJSLOC) $(ASYNCHDISTOBJSLOC)
	$(PCC) $(SHAREDOBJSLOC) $(ASYNCHDISTOBJSLOC) $(LIBS) $(FLAGS) $(OPTFLAGS) $(ALLTAUFLAGS) $(PRINTFLAG) -o ASYNCH


clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(LIBDIR)/*.so
	rm -f ASYNCH
