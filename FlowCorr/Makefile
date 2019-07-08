CXX := $(shell root-config --cxx)
CXX_FLAGS := $(shell root-config --cflags)
CXX_FLAGS += -Wall -Wextra -ggdb
LIBRARIES := $(shell root-config --glibs)

BIN := bin
SRC := src
INCLUDE := include
LIB := lib

CORR2D_D0 := corr2D_trg_d0 
CORR2D_D0_NEW := corr2D_trg_d0_new
CORR2D_PD0 := corr2D_trg_pd0 
CORR2D_NPD0 := corr2D_trg_npd0 
CORR2D_REF := corr2D_trg_ref

all: $(BIN)/$(CORR2D_D0) $(BIN)/$(CORR2D_REF) $(BIN)/$(CORR2D_NPD0) $(BIN)/$(CORR2D_D0_NEW) $(BIN)/$(CORR2D_PD0)

$(BIN)/$(CORR2D_D0): $(SRC)/corr2D_trg_d0.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_D0_NEW): $(SRC)/corr2D_trg_d0_new.cc $(SRC)/Event_new.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_PD0): $(SRC)/corr2D_trg_pd0.cc $(SRC)/Event.cxx $(SRC)/functions.cxx $(INCLUDE)/myAnaConsts.h
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_NPD0): $(SRC)/corr2D_trg_npd0.cc $(SRC)/Event.cxx $(SRC)/functions.cxx 
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_REF): $(SRC)/corr2D_trg_ref.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

clean:
	-rm $(BIN)/*

