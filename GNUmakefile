# --------------------------------------------------------------
# GNUmakefile
# Tarak Thakore, tarakstar@gmail.com, July 7, 2013
# --------------------------------------------------------------

APP      = anal_ical

SRCEXT   = cc
SRCDIR   = src
OBJDIR   = obj
BINDIR   = ./bin

SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)' )
SRCDIRS := $(shell find . -name '*.$(SRCEXT)' -exec dirname {} \; | uniq)
OBJS    := $(patsubst %.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))

G4_INC = $(shell geant4-config --cflags) -I${CLHEP_BASE_DIR}/include
G4_LIB = $(shell geant4-config --libs) -L/usr/lib64 -lconfig++ -L${CLHEP_BASE_DIR}/lib -lCLHEP

#EXTRALIBS += `root-config --glibs` -lGeom

DEBUG    = 
INCLUDES = -I./inc -I `root-config --incdir` -I/usr/local/include  $(G4_INC)
CFLAGS   = -c $(DEBUG) $(INCLUDES) -O3
LDFLAGS  = `root-config --libs` -lGeom -L/usr/local/lib -lCLHEP -lxerces-c -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint  -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lMinuit -L${CLHEP_BASE_DIR}/lib -lCLHEP -L/usr/lib64 -lc -lm -ldl -lcrypt -lpthread $(G4_LIB) # -lCint
#-L/home/apoorva/products/root61406/install/lib -lGeom -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic $(G4_LIB) 
#LDFLAGS  = `root-config --libs` -lxerces-c -lGeom $(G4_LIB)

# The library -lGeom has been added next to `root-config --libs`, wriiten originally by Tarak.
# g++ -c main.cxx -I$(root-config --incdir)
# g++ -O main.o $(root-config --libs) -lGeom -o main.exe

CC = g++
ifeq ($(CC),)
CC = g++
endif

.PHONY: all clean distclean

all: $(BINDIR)/$(APP)

$(BINDIR)/$(APP): buildrepo $(OBJS)
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) $(OBJS) $(LDFLAGS) -o $@

$(OBJDIR)/%.o: %.$(SRCEXT)
#	@echo "Generating dependencies for $<..."
#	@$(call make-depend,$<,$@,$(subst .o,.d,$@))
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

clean:
	$(RM) -r $(OBJDIR) $(APP) $(BINDIR)

distclean: clean
	$(RM) -r $(BINDIR)

buildrepo:
	@$(call make-repo)

define make-repo
   for dir in $(SRCDIRS); \
   do \
	mkdir -p $(OBJDIR)/$$dir; \
   done
endef

# usage: $(call make-depend,source-file,object-file,depend-file)
define make-depend
  $(CC) -MM       \
        -MF $3    \
        -MP       \
        -MT $2    \
        $(CFLAGS) \
        $1
endef
