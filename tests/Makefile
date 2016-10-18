#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic Licence.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#@file    Makefile
#@brief   Makefile for unit tests
#@author  Felipe Serrano

#-----------------------------------------------------------------------------
# Probably nothing below here should be changed
# TODO: use the $ORIGIN variable
# NOTE: currently, compilation with SHARED=false is not supported.
# One would need to add more flags (eg $(LDFLAGS)) when linking
#-----------------------------------------------------------------------------

UNITTESTSSRC	=	$(shell find src/ -name *.c)

# remove unittests in src/bugs if BUGS is not true
ifneq ($(BUGS),true)
	UNITTESTSSRC 	:= $(filter-out $(shell find src/bugs/ -name *.c),$(UNITTESTSSRC))
endif

# removes `src/` prefix and file.c postfix, then sorts, which removes duplicates!
TESTDIRS 	= $(sort $(dir $(patsubst src/%,%,$(UNITTESTSSRC))))
#$(warning testdirs: $(TESTDIRS))

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

SCIPDIR		=	..

#-----------------------------------------------------------------------------
# SET OPTIONS FROM SCIP BINARY
#-----------------------------------------------------------------------------

ifneq ("$(wildcard $(SCIPDIR)/bin/scip)","")
override ARCH		:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*ARCH=\(.*\)/\1/p')
override COMP		:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*COMP=\(.*\)/\1/p')
override EXPRINT	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*EXPRINT=\(.*\)/\1/p')
override GAMS		:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*GAMS=\(.*\)/\1/p')
override GMP		:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*GMP=\(.*\)/\1/p')
override IPOPT		:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*IPOPT=\(.*\)/\1/p')
override IPOPTOPT	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*IPOPTOPT=\(.*\)/\1/p')
override LPS		:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*LPS=\(.*\)/\1/p')
override LPSCHECK	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*LPSCHECK=\(.*\)/\1/p')
override LPSOPT 	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*LPSOPT=\(.*\)/\1/p')
override NOBLKBUFMEM	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*NOBLKBUFMEM=\(.*\)/\1/p')
override NOBLKMEM	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*NOBLKMEM=\(.*\)/\1/p')
override NOBUFMEM	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*NOBUFMEM=\(.*\)/\1/p')
override OPT		:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*OPT=\(.*\)/\1/p')
override OSTYPE		:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*OSTYPE=\(.*\)/\1/p')
override PARASCIP	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*PARASCIP=\(.*\)/\1/p')
override READLINE	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*READLINE=\(.*\)/\1/p')
override SHARED		:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*SHARED=\(.*\)/\1/p')
override USRARFLAGS	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*USRARFLAGS=\(.*\)/\1/p')
override USRCFLAGS	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*USRCFLAGS=\(.*\)/\1/p')
override USRCXXFLAGS	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*USRCXXFLAGS=\(.*\)/\1/p')
override USRDFLAGS	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*USRDFLAGS=\(.*\)/\1/p')
override USRFLAGS	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*USRFLAGS=\(.*\)/\1/p')
override USRLDFLAGS	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*USRLDFLAGS=\(.*\)/\1/p')
override USROFLAGS	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*USROFLAGS=\(.*\)/\1/p')
override VERSION	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*VERSION=\(.*\)/\1/p')
override ZIMPL		:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*ZIMPL=\(.*\)/\1/p')
override ZIMPLOPT	:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*ZIMPLOPT=\(.*\)/\1/p')
override ZLIB		:=$(shell $(SCIPDIR)/bin/scip -v | sed -n 's/^\s*ZLIB=\(.*\)/\1/p')
ifeq ($(SHARED), false)
$(error Please compile SCIP with SHARED=true)
endif
else
override OPT	=	dbg
override SHARED	=	true
endif

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------

include $(SCIPDIR)/make/make.project

BASESUFFIX	=	.$(BASE).$(LPS)$(EXEEXTENSION)

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

src2obj		=	$(patsubst $(SRCDIR)%,$(OBJDIR)%,$(1:.c=.o))
#$(warning obj2src: src/cons/quadratic/quadratic.c -> \
	$(call src2obj, src/cons/quadratic/quadratic.c))

#-----------------------------------------------------------------------------
# Unit test definitions
#-----------------------------------------------------------------------------


BINDIRS 	= 	$(addprefix $(BINDIR)/, $(TESTDIRS))
#$(warning bindirs: $(BINDIRS))

OBJDIRS 	= 	$(addprefix $(OBJDIR)/, $(TESTDIRS))
#$(warning objdirs: $(OBJDIRS))

# remove the .c from UNITTESTSSRC: src/lpi/bases, src/cons/quadratic/gauge, etc
UNITTESTS 	= 	$(basename $(UNITTESTSSRC))
#$(warning unittests: $(UNITTESTS))

# build object files from source files: obj/O.linux.../cons/expr/free.o, etc
UNITTESTSOBJ 	= 	$(call src2obj, $(UNITTESTSSRC))
#$(warning unittestobj: $(UNITTESTSOBJ))

# substitute src/ for $(OBJDIR) and append $(BASESUFFIX) to UNITTESTS: bin/cons/expr/free.linux..., etc
UNITTESTSBIN 	= 	$(addsuffix $(BASESUFFIX), $(patsubst $(SRCDIR)%,$(BINDIR)%,$(UNITTESTS)))
#$(warning unittestsbin: $(UNITTESTSBIN))

UNITTESTSEXECS	=	$(addsuffix $(BASESUFFIX), $(UNITTESTS))
UNITTESTSDEP	=	$(SRCDIR)/depend.$(OPT)

#-----------------------------------------------------------------------------
# Flags
# To make criterion compile without warnings
#-----------------------------------------------------------------------------

CFLAGS 		:=	$(filter-out -Wdeclaration-after-statement -Wstrict-prototypes, $(CFLAGS))
CFLAGS 		+=	-std=c99
# DFLAGS (flags for make depend) uses the -MM option which creates a make target with dependencies. -MT target-name specifies the name of target
DFLAGS 		+= 	-MT
# to find include/scip_test.h
FLAGS 		+= 	-I.

#-----------------------------------------------------------------------------
# Criterion Framework
#-----------------------------------------------------------------------------
CR_INC		=	$(wildcard Criterion/include)
CR_LIBDIR	=	$(wildcard Criterion/build)
CR_LIB		=	criterion

#-----------------------------------------------------------------------------
# Rules
# Small remainder:
# $@ refers to the target
# $^ refers to all the dependencies
# $< refers to all the first dependency
# $? refers to the dependencies that are newer than the target
# Newer dependencies are one reason why a rule is executed.
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(UNITTESTSBIN) $(UNITTESTSOBJ)
endif


.PHONY: all
all:            Criterion $(OBJDIR) $(OBJDIRS) $(BINDIRS) $(UNITTESTSBIN)
		make ctestrunner
		ctest --output-on-failure

.PHONY: lint
lint:		$(UNITTESTSSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) $(SCIPDIR)/lint/scip.lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

$(OBJDIR):
		@echo "-> Creating $@ directory"
		@-mkdir -p $@
$(OBJDIRS):
		@echo "-> Creating $@ directory"
		@-mkdir -p $@
$(BINDIRS):
		@echo "-> Creating $@ directory"
		@-mkdir -p $@

.PHONY: test
test:           $(UNITTESTSBIN)
		make ctestrunner
		ctest --output-on-failure

.PHONY: ctestrunner
ctestrunner: 	$(UNITTESTSBIN)
		@echo -e "#Automatic generated file\n" > "CTestTestfile.cmake"
		@echo -e $(foreach bin,$(UNITTESTSBIN), \
			"ADD_TEST($(shell echo $(bin:$(BASESUFFIX)=) | sed s/^bin/test/g | sed s/\\//_/g) $(bin))\n" >> "CTestTestfile.cmake")

.PHONY: clean
clean:
		@echo "-> remove main objective files"
		$(foreach dir, $(OBJDIRS) $(OBJDIR), \
			@-rm -f $(dir)/*.o)
		@echo "-> remove binaries"
		$(foreach bin, $(UNITTESTSBIN), \
			@-rm -f $(bin))

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*.c src/*.h $(SCIPDIR)/src/scip/*.c $(SCIPDIR)/src/scip/*.h;

.PHONY: depend
depend:		$(SCIPDIR)
		@- rm -f $(UNITTESTSDEP)
		$(foreach src, $(UNITTESTSSRC), \
			$(SHELL) -ec '$(DCC) $(FLAGS) $(src) \
			$(DFLAGS) $(call src2obj, $(src))' >> $(UNITTESTSDEP);)



-include	$(UNITTESTSDEP)

$(BINDIR)/%$(BASESUFFIX): $(OBJDIR)/%.o $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE)
		@echo "-> linking $@"
		$(LINKCXX) $< \
		$(LINKCXX_L)$(SCIPDIR)/lib \
		$(LINKCXX_l)$(SCIPLIB) $(LINKCXX_l)$(OBJSCIPLIB) $(LINKCXX_l)$(LPILIB) $(LINKCXX_l)$(NLPILIB) \
		$(addprefix $(LINKCXX_L),$(CR_LIBDIR)) \
		$(LINKCXX_l)$(CR_LIB) \
		$(LINKRPATH)$(realpath $(SCIPDIR)/lib) \
		$(addprefix $(LINKRPATH),$(realpath $(CR_LIBDIR))) \
		$(SANITZERFLAGS) \
		$(LINKCXX_o)$@

$(LPILIBFILE):
		$(error $@ not available. SCIP must be compiled with OPT=dbg SHARED=true IPOPT=$(IPOPT) LPS=$(LPS))
$(SCIPLIBFILE):
		$(error $@ not available. SCIP must be compiled with OPT=dbg SHARED=true IPOPT=$(IPOPT) LPS=$(LPS))
$(NLPILIBFILE):
		$(error $@ not available. SCIP must be compiled with OPT=dbg SHARED=true IPOPT=$(IPOPT) LPS=$(LPS))

$(OBJDIR)/%.o : $(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(addprefix -I,$(CR_INC)) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

# TODO: checkout v2.3.0 when it becomes available
Criterion:
	@echo -n "Criterion not found, downloading and building it in"
	@for (( count=3; $$count ; count-- )) do echo -n " $$count" ; sleep 1 ; done ; echo
	git clone https://github.com/Snaipe/Criterion; \
	cd Criterion; \
	git reset --hard 3bad3c6aaa779d9ac580ae74cfcaa245aee0ee65; \
	mkdir build; \
	cd build; \
	cmake ..; \
	make
#---- EOF --------------------------------------------------------------------
