CXX := $(shell root-config --cxx)
CXX_FLAGS := $(shell root-config --cflags)
CXX_FLAGS += -Wall -Wextra -ggdb
LIBRARIES := $(shell root-config --glibs)

BIN := bin
SRC := src
INCLUDE := include
LIB := lib

CORR2D_PD0_MULT := corr2D_trg_pd0_mult
CORR2D_PD0_MULT_BOOST := corr2D_trg_pd0_mult_boost
CORR2D_PD0_MULT_NOPU := corr2D_trg_pd0_mult_noPU
CORR2D_PD0_MULT_LOOSE := corr2D_trg_pd0_mult_loose
CORR2D_PD0_MULT_TIGHT := corr2D_trg_pd0_mult_tight
CORR2D_REF_MULT := corr2D_trg_ref_mult
CORR2D_PD0      := corr2D_trg_pd0
CORR2D_NPD0     := corr2D_trg_npd0
CORR2D_NPD0_NEW := corr2D_trg_npd0_new
CORR2D_NPD0_NEW_MVA := corr2D_trg_npd0_new_mva
CORR2D_NPD0_NEW_LOOSE := corr2D_trg_npd0_new_loose
CORR2D_NPD0_NEW_TIGHT := corr2D_trg_npd0_new_tight
CORR2D_REF 		 := corr2D_trg_ref
COUNT := countEvt

all: $(BIN)/$(CORR2D_PD0_MULT) $(BIN)/$(CORR2D_REF_MULT) $(BIN)/$(CORR2D_PD0) $(BIN)/$(CORR2D_NPD0) $(BIN)/$(CORR2D_REF) $(BIN)/$(CORR2D_NPD0_NEW) $(BIN)/$(CORR2D_PD0_MULT_LOOSE) $(BIN)/$(CORR2D_PD0_MULT_TIGHT)  $(BIN)/$(CORR2D_NPD0_NEW_LOOSE) $(BIN)/$(CORR2D_NPD0_NEW_TIGHT) $(BIN)/$(COUNT) $(BIN)/$(CORR2D_PD0_MULT_NOPU) $(BIN)/$(CORR2D_PD0_MULT_BOOST) $(BIN)/$(CORR2D_NPD0_NEW_MVA) 

$(BIN)/$(COUNT): $(SRC)/countEvt.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_PD0_MULT): $(SRC)/corr2D_trg_pd0_mult.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_PD0_MULT_BOOST): $(SRC)/corr2D_trg_pd0_mult_boost.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_PD0_MULT_NOPU): $(SRC)/corr2D_trg_pd0_mult_noPU.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_PD0_MULT_LOOSE): $(SRC)/corr2D_trg_pd0_mult_loose.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_PD0_MULT_TIGHT): $(SRC)/corr2D_trg_pd0_mult_tight.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_REF_MULT): $(SRC)/corr2D_trg_ref_mult.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_PD0): $(SRC)/corr2D_trg_pd0.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_NPD0): $(SRC)/corr2D_trg_npd0.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_NPD0_NEW): $(SRC)/corr2D_trg_npd0_new.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_NPD0_NEW_LOOSE): $(SRC)/corr2D_trg_npd0_new_loose.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_NPD0_NEW_TIGHT): $(SRC)/corr2D_trg_npd0_new_tight.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_NPD0_NEW_MVA): $(SRC)/corr2D_trg_npd0_new_mva.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

$(BIN)/$(CORR2D_REF): $(SRC)/corr2D_trg_ref.cc $(SRC)/Event.cxx $(SRC)/functions.cxx
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

clean:
	-rm $(BIN)/*

