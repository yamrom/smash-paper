# Copyright Peter Andrews CSHL 2013

# What to build by default
EXEC	= fastqs_to_sam mappability_tag mummer mummer-medium mummer-long 

EXTRA_OPTS	= -pthread
#
EXTRA_LD	= -pthread
# 
# Shared rules
include shared.mk

# Linking object files into executable for each int size
fastqs_to_sam	: fastqs_to_sam.o strings.o util.o
mappability_tag	: mappability_tag.o strings.o util.o
MUMMER	= mummer.o fasta.o locked.o longSA.o memsam.o qsufsort.o query.o util.o
mummer		: $(MUMMER)
mummer-medium	: $(MUMMER:.o=.om) ; $(CXX) $(LDFLAGS) -o $@ $^
mummer-long	: $(MUMMER:.o=.ol) ; $(CXX) $(LDFLAGS) -o $@ $^

# Compilation of source files into object files for each int size
%.om		: %.cpp	; $(CXX) $(CXXFLAGS)   -c -DSINTS -o $@ $<
%.ol		: %.cpp	; $(CXX) $(CXXFLAGS)   -c -DSINTS -DUINTS -o $@ $<

# Dependency file inclusion and generation - overrides DOPT from shared.mk
TRANS	= '$$d.=$$_;END{$$_=$$d;s/\.o:/.om:/;print;s/\.om:/.ol:/;print;}'
DOPT	= $(STD) -MM -MT $(subst .d,.o,$@) -MF >(perl -pe $(TRANS) > $@)

