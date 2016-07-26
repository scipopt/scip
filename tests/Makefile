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
# OPTIONS
#-----------------------------------------------------------------------------

override OPT	=	dbg
override SHARED	=	true

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

SCIPDIR		=	..

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------

include $(SCIPDIR)/make/make.project
#$(warning OPT=$(OPT))

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
CR_INC		=	Criterion/include
CR_LIBDIR	=	Criterion/build/
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
		$(LINKCXX_L)$(CR_LIBDIR) \
		$(LINKCXX_l)$(CR_LIB) \
		$(LINKRPATH)$(realpath $(SCIPDIR)/lib) \
		$(LINKRPATH)$(realpath $(CR_LIBDIR)) \
		$(LINKCXX_o)$@

$(LPILIBFILE):
		$(error SCIP must be compiled with OPT=dbg IPOPT=true SHARED=true)
$(SCIPLIBFILE):
		$(error SCIP must be compiled with OPT=dbg IPOPT=true SHARED=true)
$(NLPILIBFILE):
		$(error To compile the tests with IPOPT=true, SCIP must be compiled with OPT=dbg IPOPT=true SHARED=true)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) -I$(CR_INC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

Criterion:
	@echo "Criterion not found, building it"
	git clone -b v2.2.2 https://github.com/Snaipe/Criterion
	cd Criterion; \
	mkdir build; \
	cd build; \
	cmake ..; \
	make
#---- EOF --------------------------------------------------------------------
