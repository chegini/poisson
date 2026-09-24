KASKADE7 := /kaskade7

include $(KASKADE7)/Makefile.Local
include $(KASKADE7)/Makefile.Rules


OPTFLAGS = # -O3 -DNDEBUG
DEBUGFLAGS =  -g

default: bddc bddc-add

OBJ = bddc.o bddc-add.o

# Include the dependency files for the objects. This triggers the building of
# dependency files, and their presence ensures that the implicit compilation
# rule for creating objects (defined in Makefile.Rules) is found by make.
include $(patsubst %.o, %.d, $(OBJ))


bddc: bddc.o $(KASKADE7)/libs/*.o $(KASKADE7)/libs/*.a
	$(CXX) $< $(KASKADE7)/libs/*.o $(KASKADELIB) $(DUNELIB) $(UGLIB) $(BOOSTLIB) $(DIRECTSOLVERLIB) \
	$(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) $(DEBUGFLAGS) -o $@

bddc-add: bddc-add.o $(KASKADE7)/libs/*.a
	$(CXX) $< $(KASKADELIB) $(DUNELIB) $(UGLIB) $(BOOSTLIB) $(DIRECTSOLVERLIB) \
	$(BLASLIB) $(FTNLIB) $(NUMALIB) $(LINKFLAGS) $(DEBUGFLAGS) -o $@





clean:
	rm -f gccerr.txt $(OBJ) $(OBJ:.o=.d) bddc bddc-add *.vtu *.pvtu
