SRC := src
OBJ := build
MOD := modules
DEBS := precision.mod hung.mod fingerprint.mod
FC = gfortran

ifeq ($(FC), gfortran)
linking_flags = -Ofast -llapack -lblas -fopenmp -J$(OBJ)/$(MOD)
compile_flags = $(linking_flags)
else
linking_flags = -O3 -mkl -qopenmp -module $(OBJ)/$(MOD)
compile_flags = $(linking_flags)
endif

MODULE := $(SRC)/precision.f90 $(SRC)/hung.f90 $(SRC)/fingerprint.f90 $(SRC)/rot.f90 $(SRC)/simplex.f90
OBJECTS := $(OBJ)/$(MOD)/precision.o $(OBJ)/$(MOD)/hung.o $(OBJ)/$(MOD)/fingerprint.o $(OBJ)/$(MOD)/rot.o $(OBJ)/$(MOD)/simplex.o
SOURCES1 := $(SRC)/fp_distance.f90
SOURCES2 := $(SRC)/fp_for_clusters.f90
SOURCES3 := $(SRC)/fp_for_crystals.f90
SOURCES4 := $(SRC)/test_derivative.f90
SOURCES5 := $(SRC)/test_invariance.f90
SOURCES6 := $(SRC)/fp_simplex.f90

default: $(OBJECTS)

$(OBJ)/$(MOD)/%.o: $(SRC)/%.f90 | bin
	$(FC) -c $< -o $@ $(compile_flags)


fp_distance: $(SOURCES1) $(OBJECTS)
	$(FC) $(OBJ)/$(MOD)/fingerprint.o $(OBJ)/$(MOD)/hung.o $(OBJ)/$(MOD)/precision.o $(SOURCES1) -o $(OBJ)/fp_distance.x $(compile_flags)

fp_crystal: $(SOURCES2) $(OBJECTS)
	$(FC) $(OBJ)/$(MOD)/fingerprint.o $(OBJ)/$(MOD)/hung.o $(OBJ)/$(MOD)/precision.o $(SOURCES2) -o $(OBJ)/fp_crystal.x $(compile_flags)

fp_cluster: $(SOURCES3) $(OBJECTS)
	$(FC) $(OBJ)/$(MOD)/fingerprint.o $(OBJ)/$(MOD)/hung.o $(OBJ)/$(MOD)/precision.o $(SOURCES3) -o $(OBJ)/fp_cluster.x $(compile_flags)

test_derivative: $(SOURCES4) $(OBJECTS)
	cp $(SRC)/plot_results.gnuplot $(OBJ)/plot_results.gnuplot; $(FC) $(OBJ)/$(MOD)/fingerprint.o $(SOURCES4) -o $(OBJ)/test_derivative.x $(compile_flags)

test_invariance: $(SOURCES5) $(OBJECTS)
	$(FC) $(OBJ)/$(MOD)/fingerprint.o $(OBJ)/$(MOD)/rot.o $(OBJ)/$(MOD)/hung.o $(OBJ)/$(MOD)/precision.o $(SOURCES5) -o $(OBJ)/test_invariance.x $(compile_flags)

simplex: $(SOURCES6) $(OBJECTS)
	$(FC) $(OBJ)/$(MOD)/fingerprint.o $(OBJ)/$(MOD)/simplex.o $(SOURCES6) -o $(OBJ)/fp_simplex.x $(compile_flags)


bin:
	mkdir -p build
	mkdir -p build/modules

clean:
	rm -rf $(OBJ)
