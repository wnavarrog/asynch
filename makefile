#Compilers
PCC = mpicc

#Flags C standard
FLAGS = -std=gnu99

#Flags for performance
OPTFLAGS = -O3 -march=nocona

#Locations
LIBDIR = ./libs
OBJDIR = ./objects

#Headers and libraries
LIBS = -lm -lpq

#Objects
SHAREDOBJS = $(addprefix $(OBJDIR)/,rkmethods.o problems.o mathmethods.o riversys.o sort.o comm.o system.o processdata.o partition.o definetype.o misc.o \
	rainfall.o solvers.o io.o forcings.o compression.o date_manip.o asynch_interface.o modeloutputs.o data_types.o)
LIBPYOBJS = $(addprefix $(OBJDIR)/,asynch_interface_py.o)
ASYNCHDISTOBJS = $(addprefix $(OBJDIR)/,asynchdist.o)

#How to compile and link
$(OBJDIR)/%.o: %.c
	$(PCC) -c $*.c $(FLAGS) $(OPTFLAGS) $(EXTRA_FLAGS) -o $(OBJDIR)/$*.o

ASYNCHLIB:
	$(MAKE) ASYNCHLIB_TARGET EXTRA_FLAGS=-fPIC

ASYNCHLIB_TARGET: $(SHAREDOBJS)
	$(PCC) $(SHAREDOBJS) $(LIBS) $(FLAGS) $(OPTFLAGS) -shared -Wl,-soname,libasynch.so -o $(LIBDIR)/libasynch.so

ASYNCHLIB_PY:
	$(MAKE) ASYNCHLIB_PY_TARGET EXTRA_FLAGS="-fPIC -I/usr/include/python2.4/ -L/usr/lib/python2.4/"

ASYNCHLIB_PY_TARGET: $(SHAREDOBJS) $(LIBPYOBJS)
	$(PCC) $(SHAREDOBJS) $(LIBPYOBJS) $(LIBS) $(FLAGS) $(OPTFLAGS) -lpython2.4 -shared -Wl,-soname,libasynch_py.so -o $(LIBDIR)/libasynch_py.so

ASYNCH: $(SHAREDOBJS) $(ASYNCHDISTOBJS)
	$(PCC) $(SHAREDOBJS) $(ASYNCHDISTOBJS) $(LIBS) $(FLAGS) $(OPTFLAGS) -o ASYNCH


clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(LIBDIR)/*.so
	rm -f ASYNCH

