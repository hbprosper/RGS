# ----------------------------------------------------------------------------
# Build libRGS.so
# Created 27 Feb 2013 HBP & SS
#         30 May 2015 HBP - standardize structure (src, lib, include) 
# ----------------------------------------------------------------------------
ifndef ROOTSYS
	$(error *** Please set up Root)
endif
ROOFIT	:= $(ROOTSYS)
# ----------------------------------------------------------------------------
NAME	:= RGS
incdir	:= include
srcdir	:= src
libdir	:= lib

# create lib directory if one does not exist
$(shell mkdir -p lib)

# get lists of sources

SRCS	:=  	$(srcdir)/RGS.cc

CINTSRCS:= $(wildcard $(srcdir)/*_dict.cc)

OTHERSRCS:= $(filter-out $(CINTSRCS) $(SRCS),$(wildcard $(srcdir)/*.cc))

# list of dictionaries to be created
DICTIONARIES	:= $(SRCS:.cc=_dict.cc)

# get list of objects
OBJECTS		:= $(SRCS:.cc=.o) $(OTHERSRCS:.cc=.o) $(DICTIONARIES:.cc=.o)
#say := $(shell echo "DICTIONARIES:     $(DICTIONARIES)" >& 2)
#say := $(shell echo "" >& 2)
#say := $(shell echo "SRCS: $(SRCS)" >& 2)
#say := $(shell echo "OBJECTS: $(OBJECTS)" >& 2)
#$(error bye)
# ----------------------------------------------------------------------------
ROOTCINT	:= rootcint
# check for clang++, otherwise use g++
COMPILER	:= $(shell which clang++)
ifneq ($(COMPILER),)
CXX		:= clang++
LD		:= clang++
else
CXX		:= g++
LD		:= g++
endif

CPPFLAGS	:= -I. -I$(incdir)
CXXFLAGS	:= -O -Wall -fPIC -g -ansi -Wshadow -Wextra \
$(shell root-config --cflags)
LDFLAGS		:= -g
# ----------------------------------------------------------------------------
# which operating system?
OS := $(shell uname -s)
ifeq ($(OS),Darwin)
	LDFLAGS += -dynamiclib
	LDEXT	:= .dylib
else
	LDFLAGS	+= -shared
	LDEXT	:= .so
endif
LDFLAGS += $(shell root-config --ldflags)
LIBS	:= $(shell root-config --libs) -lTreePlayer -lPyROOT
LIBRARY	:= $(libdir)/lib$(NAME)$(LDEXT)
# ----------------------------------------------------------------------------
all: $(LIBRARY)

$(LIBRARY)	: $(OBJECTS)
	@echo ""
	@echo "=> Linking shared library $@"
	$(LD) $(LDFLAGS) $^ $(LIBS)  -o $@

$(OBJECTS)	: %.o	: 	%.cc
	@echo "=> Compiling $<"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(DICTIONARIES)	: $(srcdir)/%_dict.cc	: $(incdir)/%.h
	@echo "=> Building dictionary $@"
	$(ROOTCINT) -f $@ -c $(CPPFLAGS) $^
	find $(srcdir) -name "*.pcm" -exec mv {} $(libdir) \;

clean:
	rm -rf $(libdir)/* $(srcdir)/*_dict*.* $(srcdir)/*.o 
