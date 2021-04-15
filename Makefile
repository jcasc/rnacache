#--------------------------------------------------------------------
# htslib options
#--------------------------------------------------------------------
HTS_INCLUDE = dep/htslib/include/htslib
HTS_LIB = dep/htslib/lib/libhts.a


#--------------------------------------------------------------------
# names & flags
#--------------------------------------------------------------------
ARTIFACT     = rnacache

COMPILER     = $(CXX)
DIALECT      = -std=c++14
WARNINGS     = -Wall -Wextra -Wpedantic
OPTIMIZATION = -O3
#-march native -fomit-frame-pointer 

ifeq ($(RC_BAM), TRUE)
    DEP_CXXFLAGS = -I$(HTS_INCLUDE) -DRC_BAM
endif

REL_CXXFLAGS = $(INCLUDES) $(MACROS) $(DIALECT) $(OPTIMIZATION) $(WARNINGS) $(DEP_CXXFLAGS)
DBG_CXXFLAGS = $(INCLUDES) $(MACROS) $(DIALECT) -O0 -g $(WARNINGS) $(DEP_CXXFLAGS)
PRF_CXXFLAGS = $(INCLUDES) $(MACROS) $(DIALECT) $(OPTIMIZATION) -g $(WARNINGS) $(DEP_CXXFLAGS)

ifeq ($(RC_BAM), TRUE)
    STATIC_LIBS  = $(HTS_LIB)
    DEP_LDFLAGS  = -lz -llzma -lbz2
endif

REL_LDFLAGS  = -pthread -s $(DEP_LDFLAGS)
DBG_LDFLAGS  = -pthread $(DEP_LDFLAGS)
PRF_LDFLAGS  = -pthread $(DEP_LDFLAGS)

REL_DIR     = build_release
DBG_DIR     = build_debug
PRF_DIR     = build_profile

REL_ARTIFACT     = $(ARTIFACT)
DBG_ARTIFACT     = $(ARTIFACT)_debug
PRF_ARTIFACT     = $(ARTIFACT)_prf


#--------------------------------------------------------------------
# main targets
#--------------------------------------------------------------------
.PHONY: all release debug profile clean 
	
release: $(REL_DIR) $(REL_ARTIFACT)
debug:   $(DBG_DIR) $(DBG_ARTIFACT)
profile: $(PRF_DIR) $(PRF_ARTIFACT)

all: release debug profile

clean : 
	rm -rf build_*
	rm -f *.exe
	rm -f *.json
	rm -f 
	rm -f $(REL_ARTIFACT)
	rm -f $(DBG_ARTIFACT)
	rm -f $(PRF_ARTIFACT)


#--------------------------------------------------------------------
# dependencies
#--------------------------------------------------------------------
HEADERS = \
          src/alignment.h \
          src/batch_processing.h \
          src/bitmanip.h \
          src/candidates.h \
          src/chunk_allocator.h \
          src/classification.h \
          src/classification_statistics.h \
          src/classify_common.h \
          src/cmdline_utility.h \
          src/config.h \
          src/database.h \
          src/dna_encoding.h \
          src/filesys_utility.h \
          src/hash_dna.h \
          src/hash_int.h \
          src/hash_multimap.h \
          src/io_error.h \
          src/io_options.h \
          src/io_serialize.h \
          src/matches_per_target.h \
          src/modes.h \
          src/options.h \
          src/printing.h \
          src/querying.h \
          src/sequence_io.h \
          src/sequence_view.h \
          src/stat_combined.h \
          src/stat_confusion.h \
          src/stat_moments.h \
          src/string_utils.h \
          src/taxonomy.h \
          src/taxonomy_io.h \
          src/timer.h \
          src/version.h \
          dep/edlib.h

SOURCES = \
          src/classify_common.cpp \
          src/classify_rna.cpp \
          src/cmdline_utility.cpp \
          src/database.cpp \
          src/filesys_utility.cpp \
          src/main.cpp \
          src/mode_build.cpp \
          src/mode_help.cpp \
          src/mode_info.cpp \
          src/mode_merge.cpp \
          src/mode_query.cpp \
          src/options.cpp \
          src/printing.cpp \
          src/sequence_io.cpp \
          src/taxonomy_io.cpp \
          dep/edlib.cpp


#--------------------------------------------------------------------
# subtarget generator for out-of-place build
#--------------------------------------------------------------------
PLAIN_SRCS = $(notdir $(SOURCES))
PLAIN_OBJS = $(PLAIN_SRCS:%.cpp=%.o)

# $(1):    $(2)       $(3)     $(4)
# artifact build_dir cxxflags ldflags
define make_subtargets =
$(1): $$(PLAIN_OBJS:%=$(2)/%)
	$(COMPILER) -o $(1) $$(PLAIN_OBJS:%=$(2)/%) $(STATIC_LIBS) $(4)

$(2):
	mkdir $(2) 
    
$(2)/classify_rna.o : src/classify_rna.cpp $(HEADERS) 
	$(COMPILER) $(3) -c $$< -o $$@
	
$(2)/classify_common.o : src/classify_common.cpp $(HEADERS) 
	$(COMPILER) $(3) -c $$< -o $$@
	
$(2)/cmdline_utility.o : src/cmdline_utility.cpp src/cmdline_utility.h 
	$(COMPILER) $(3) -c $$< -o $$@
	
$(2)/filesys_utility.o : src/filesys_utility.cpp src/filesys_utility.h 
	$(COMPILER) $(3) -c $$< -o $$@
	
$(2)/main.o : src/main.cpp src/modes.h 
	$(COMPILER) $(3) -c $$< -o $$@

$(2)/mode_build.o : src/mode_build.cpp $(HEADERS) 
	$(COMPILER) $(3) -c $$< -o $$@

$(2)/mode_help.o : src/mode_help.cpp src/modes.h src/filesys_utility.h 
	$(COMPILER) $(3) -c $$< -o $$@
	
$(2)/mode_info.o : src/mode_info.cpp $(HEADERS) 
	$(COMPILER) $(3) -c $$< -o $$@

$(2)/mode_merge.o : src/mode_merge.cpp $(HEADERS) 
	$(COMPILER) $(3) -c $$< -o $$@
	
$(2)/mode_query.o : src/mode_query.cpp $(HEADERS) 
	$(COMPILER) $(3) -c $$< -o $$@
	
$(2)/printing.o : src/printing.cpp $(HEADERS) 
	$(COMPILER) $(3) -c $$< -o $$@
	
$(2)/options.o : src/options.cpp $(HEADERS)  
	$(COMPILER) $(3) -c $$< -o $$@

$(2)/sequence_io.o : src/sequence_io.cpp src/sequence_io.h src/io_error.h 
	$(COMPILER) $(3) -c $$< -o $$@
	
$(2)/taxonomy_io.o : src/taxonomy_io.cpp src/taxonomy_io.h src/cmdline_utility.h src/filesys_utility.h 
	$(COMPILER) $(3) -c $$< -o $$@

$(2)/database.o : src/database.cpp $(HEADERS)
	$(COMPILER) $(3) -c $$< -o $$@

$(2)/edlib.o : dep/edlib.cpp dep/edlib.h
	$(COMPILER) $(3) -c $$< -o $$@
	
endef


#--------------------------------------------------------------------
# generate all subtargets
#--------------------------------------------------------------------
$(eval $(call make_subtargets,$(REL_ARTIFACT),$(REL_DIR),$(REL_CXXFLAGS),$(REL_LDFLAGS)))
$(eval $(call make_subtargets,$(DBG_ARTIFACT),$(DBG_DIR),$(DBG_CXXFLAGS),$(DBG_LDFLAGS)))
$(eval $(call make_subtargets,$(PRF_ARTIFACT),$(PRF_DIR),$(PRF_CXXFLAGS),$(PRF_LDFLAGS)))
