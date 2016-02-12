#
# shared.mk
#
# Common definitions for compiling C++ code with gcc
#
# Copyright Peter Andrews CSHL 2014
#

# What shell to use
SHELL	= bash

# What to build by default
.PHONY	: all
all	: $(EXEC) lint

# Third party packages to link to
THIRD			= /data/software
tcmalloc		= # defined
ifdef tcmalloc
	tclib		= $(THIRD)/gperf/2.1.90/lib
	unwindlib	= $(THIRD)/libunwind/0.99-beta/lib
	tclibs		= -ltcmalloc # -lprofiler
	extralibs	= :$(tclib):$(unwindlib)
	extralibdirs	= $(tclibs) -L$(tclib) -L$(unwindlib)
endif

# Compiler choice
newgcc			= defined
ifdef newgcc
	GCC_DIR		= /data/software/gcc/4.9.2/rtf
	PATH		= $(GCC_DIR)/bin:/bin:/usr/bin
	LD_LIBRARY_PATH	= $(GCC_DIR)/lib64$(extralibs)
	STD		= -std=c++11 -pedantic
endif

# Compilation of source files into object files
FAST	= -Ofast -march=native -m64
#DEBUG	= -ggdb -pg 
WARN	= -Wall -Wextra -Wc++11-compat -Woverloaded-virtual -Wsign-promo -Wdouble-promotion -Wformat-security -Winit-self -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wsync-nand -Wuninitialized -Wunknown-pragmas -Wfloat-equal -Wcast-qual -Wcast-align -Wzero-as-null-pointer-constant -Wlogical-op -Wpacked -Wredundant-decls -Wvla -Wdisabled-optimization -Wundef -Wshadow
CXXFLAGS	= $(DEBUG) $(STD) $(FAST) $(WARN) $(EXTRA_OPTS)

# Linking of object files into executables - assume c++ for all
CC     	=  g++
LDFLAGS	= -Xlinker -rpath=$(LD_LIBRARY_PATH) $(extralibdirs) $(EXTRA_LD)

# Linting of all code
CODE 	= $(wildcard *.h) $(wildcard *.cpp) $(wildcard *.cc)
LINT	= ./cpplint.py --filter=-readability/streams,-runtime/references,-runtime/printf,-build/header_guard
FILT	= grep -v -e 'Done process' -e 'Total err'
lint 	: $(CODE:=.lint)
%.lint	: % ; @ $(LINT) $< 2>&1 | tee >($(FILT) 1>&2) > $@

# Parallel remake without stdout or repetition in stderr, to see stderr better
UNIQ 	= perl -ne 'print unless $$seen{$$_}++'
warn	: ; @ $(MAKE) -s clean ; $(MAKE) -j -s 2> >($(UNIQ) 1>&2)

# Clean-up of directory
clean	: ; rm -Rf *.o{,m,l} *.d *.lint *~ $(EXEC)

# Dependency file inclusion and generation
DOPT	= $(STD) -MM -MT $(subst .d,.o,$@) -MF $@
%.d	: %.cpp ; @ $(CXX) $(STD) $(DOPT) $<
-include $(subst .cpp,.d,$(wildcard *.cpp))
%.d	: %.cc ; @ $(CXX) $(STD) $(DOPT) $<
-include $(subst .cc,.d,$(wildcard *.cc))

# Package code as zip in home directory for sending via email
DIR	:= $(notdir $(PWD))
zip	:  clean ; rm -f ~/$(DIR).zip ; (cd .. ; zip -r ~/$(DIR) $(DIR))
