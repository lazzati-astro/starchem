# Makefile for starchem
#
#

# Project Name
PROJECT := starchem

# output directories
prefix ?= /usr/local
ROOT_BUILD_DIR := build

# local dir
LOCAL_DIR := /home/mauneyc/local

BOOST_DIR := $(LOCAL_DIR)/boost
PLOG_DIR  := $(LOCAL_DIR)/plog

# build dir
DEBUG ?= 0
ifeq ($(DEBUG),1)
	BUILD_DIR := $(ROOT_BUILD_DIR)/debug
else
	BUILD_DIR := $(ROOT_BUILD_DIR)/release
endif

NAME := $(BUILD_DIR)/bin/$(PROJECT)

# compiler
AR := ar rcs
CXX ?= g++

# libraries
LIBS := m netcdf_c++4 boost_program_options
LIB_DIR += $(BOOST_DIR)/lib
INC_DIR += ./include $(BOOST_DIR)/include $(PLOG_DIR)/include

# compiler flags
OMP_FLAGS := 
COMMON_FLAGS := $(addprefix -I, $(INC_DIR))
CXX_FLAGS := -std=c++1z $(OMP_FLAGS)
LINK_FLAGS := $(OMP_FLAGS)
LD_FLAGS := $(addprefix -l,$(LIBS)) $(addprefix -L, $(LIB_DIR))

ifeq ($(DEBUG),1)
	COMMON_FLAGS += -pg 
else
	COMMON_FLAGS += -O3 -march=native
endif

# source files
CXX_SRC := $(shell find src -name "*.cpp")
HXX_SRC := $(shell find include -name "*.h")

# object files
CXX_OBJ := $(addprefix $(BUILD_DIR)/objs/,$(CXX_SRC:.cpp=.o))

.PHONY: all clean

all: $(NAME)

$(NAME) : $(CXX_OBJ)
	@echo [ Linking ] $@
	@mkdir -p $(BUILD_DIR)/bin
	@$(CXX) -o $@ $(CXX_OBJ) $(COMMON_FLAGS) $(LD_FLAGS) $(LINK_FLAGS)

$(BUILD_DIR)/objs/%.o: %.cpp $(HXX_SRC)
	@echo [ CXX ] $<
	@$(foreach d, $(subst /, ,${@D}), mkdir -p $d && cd $d && ):
	@$(CXX) $(CXXFLAGS) $(COMMON_FLAGS) -c -o $@ $<

clean:
	@rm -rf $(ROOT_BUILD_DIR)

