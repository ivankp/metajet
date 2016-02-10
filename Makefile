CC  := gcc
CXX := g++

CXXFLAGS := -std=c++11 -Wall -g -Isrc
# $(shell root-config --cflags)
LIBS :=
# $(shell root-config --libs)

.PHONY: all clean

NODEPS := clean

SRCDIR := src
BLDDIR := .build
EXEDIR := bin

BLDDIRS := $(shell find $(SRCDIR) -type d | sed "s|^$(SRCDIR)|$(BLDDIR)|")
EXEDIRS := $(shell find $(SRCDIR) -type d | sed "s|^$(SRCDIR)|$(EXEDIR)|")

SRCS := $(shell find $(SRCDIR) -type f -name "*.cc")
OBJS := $(patsubst $(SRCDIR)/%.cc,$(BLDDIR)/%.o,$(SRCS))
DEPS := $(patsubst %.o,%.d,$(OBJS))

GREP_EXE := grep -r '^ *int *main *(' $(SRCDIR) | cut -d':' -f1
EXES := $(patsubst $(SRCDIR)/%.cc,$(EXEDIR)/%,$(shell $(GREP_EXE)))
EXE_DEPS := $(patsubst $(EXEDIR)/%,$(BLDDIR)/%.d,$(EXES))

all: $(EXES)

#Don't create dependencies when we're cleaning, for instance
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
-include $(DEPS)
endif

# object and executable dependencies
$(EXE_DEPS): $(BLDDIR)/%.d: $(SRCDIR)/%.cc
	@echo DEP $(notdir $@)
	@$(CXX) $(CXXFLAGS) -MM -MT '$(<:$(SRCDIR)/%.cc=$(BLDDIR)/%.o)' $< -MF $@
	@objs=""; \
	for f in `sed 's| \\\\$$||g;s| $(SRCDIR)/|\n|g' $@ \
	      | sed '/\.hh$$/!d;s|\.hh$$||'`; do \
	  [ -s "$(SRCDIR)/$$f.cc" ] && objs="$${objs} $(BLDDIR)/$$f.o"; \
	done; \
	if [ "$${objs}" ]; then \
	  echo "$(@:$(BLDDIR)/%.d=$(EXEDIR)/%):$${objs}" >> $@; \
	fi

# object dependencies
$(BLDDIR)/%.d: $(SRCDIR)/%.cc
	@echo DEP $(notdir $@)
	@$(CXX) $(CXXFLAGS) -MM -MT '$(<:$(SRCDIR)/%.cc=$(BLDDIR)/%.o)' $< -MF $@

# compile objects
$(BLDDIR)/%.o :
	@echo CXX $(notdir $@)
	@$(CXX) -c -I$(SRCDIR) $(CXXFLAGS) $< -o $@

# link executables
$(EXEDIR)/% : $(BLDDIR)/%.o
	@echo LD $(notdir $@)
	@$(CXX) $(filter %.o,$^) -o $@ $(LIBS)

# directories as order-only-prerequisites
$(OBJS) $(DEPS): | $(BLDDIRS)
$(EXES): | $(EXEDIRS)

# make directories
$(BLDDIRS) $(EXEDIRS):
	mkdir $@

clean:
	@rm -vfr $(BLDDIR) $(EXEDIR)
