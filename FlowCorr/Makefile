CXX := $(shell root-config --cxx)
CXX_FLAGS := $(shell root-config --cflags)
CXX_FLAGS += -Wall -Wextra -ggdb
LIBRARIES := $(shell root-config --glibs)

BIN := bin
SRC := src
INCLUDE := include
LIB := lib

CORR2D_PD0_MULT := corr2D_trg_pd0_mult
CORR2D_REF_MULT := corr2D_trg_ref_mult

all: $(BIN)/$(CORR2D_D0_MULT) $(BIN)/$(CORR2D_REF_MULT)

$(BIN)/$(CORR2D_PD0_MULT): $(SRC)/corr2D_trg_pd0_mult.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_REF_MULT): $(SRC)/corr2D_trg_ref_mult.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

clean:
	-rm $(BIN)/*

