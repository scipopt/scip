#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
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
# put an @ at the end of each line of scip -v to use as a marker for sed
SCIPVERSION			:=$(shell $(SCIPDIR)/bin/scip -v | sed -e 's/$$/@/')
override ARCH		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ARCH=\([^@]*\).*/\1/')
override COMP		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* COMP=\([^@]*\).*/\1/')
override EXPRINT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* EXPRINT=\([^@]*\).*/\1/')
override GAMS		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* GAMS=\([^@]*\).*/\1/')
override GMP		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* GMP=\([^@]*\).*/\1/')
override SYM		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SYM=\([^@]*\).*/\1/')
override IPOPT		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* IPOPT=\([^@]*\).*/\1/')
override IPOPTOPT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* IPOPTOPT=\([^@]*\).*/\1/')
override LPS		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPS=\([^@]*\).*/\1/')
override LPSCHECK	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPSCHECK=\([^@]*\).*/\1/')
override LPSOPT 	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* LPSOPT=\([^@]*\).*/\1/')
override NOBLKBUFMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBLKBUFMEM=\([^@]*\).*/\1/')
override NOBLKMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBLKMEM=\([^@]*\).*/\1/')
override NOBUFMEM	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* NOBUFMEM=\([^@]*\).*/\1/')
override OPT		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* OPT=\([^@]*\).*/\1/')
override OSTYPE		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* OSTYPE=\([^@]*\).*/\1/')
override PARASCIP	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* PARASCIP=\([^@]*\).*/\1/')
override READLINE	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* READLINE=\([^@]*\).*/\1/')
override SANITIZE	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SANITIZE=\([^@]*\).*/\1/')
override SHARED		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SHARED=\([^@]*\).*/\1/')
override USRARFLAGS	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* USRARFLAGS=\([^@]*\).*/\1/')
override USRCFLAGS	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* USRCFLAGS=\([^@]*\).*/\1/')
override USRCXXFLAGS	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* USRCXXFLAGS=\([^@]*\).*/\1/')
override USRDFLAGS	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* USRDFLAGS=\([^@]*\).*/\1/')
override USRFLAGS	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* USRFLAGS=\([^@]*\).*/\1/')
override USRLDFLAGS	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* USRLDFLAGS=\([^@]*\).*/\1/')
override USROFLAGS	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* USROFLAGS=\([^@]*\).*/\1/')
override VERSION	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* VERSION=\([^@]*\).*/\1/')
override ZIMPL		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZIMPL=\([^@]*\).*/\1/')
override ZIMPLOPT	:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZIMPLOPT=\([^@]*\).*/\1/')
override ZLIB		:=$(shell echo "$(SCIPVERSION)" | sed -e 's/.* ZLIB=\([^@]*\).*/\1/')
ifeq ($(SHARED), false)
$(error Please compile SCIP with SHARED=true)
endif
ifneq ($(OPT), dbg)
$(error Please compile SCIP with OPT=dbg)
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
# filter out some flags to make criterion compile without warnings with gcc
# and make the tests run with clang.
# TODO: in newer version, check if we still need this
#-----------------------------------------------------------------------------

CFLAGS 		:=	$(filter-out -Wdeclaration-after-statement -Wstrict-prototypes, $(CFLAGS))
ifeq ($(COMP), clang)
CFLAGS 		:=	$(filter-out -Wcast-align, $(CFLAGS))
endif
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
		setarch `uname -m` -R ctest --output-on-failure

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

# run the tests with address space layout randomization disabled, as this isn't
# compatible with Criterions parametrized tests (https://github.com/Snaipe/Criterion/issues/208)
.PHONY: test
test:           $(UNITTESTSBIN)
		make ctestrunner
		setarch `uname -m` -R ctest --output-on-failure

.PHONY: ctestrunner
ctestrunner: 	$(UNITTESTSBIN)
		@echo -e "#Automatic generated file\n" > "CTestTestfile.cmake"
		@echo -e $(foreach bin,$(UNITTESTSBIN), \
			"ADD_TEST($(shell echo $(bin:$(BASESUFFIX)=) | sed s/^bin/test/g | sed s/\\//_/g) $(bin))\n" >> "CTestTestfile.cmake")

.PHONY: clean
clean:
		@echo "-> remove main object files"
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

$(BINDIR)/%$(BASESUFFIX): $(OBJDIR)/%.o $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(TPILIBFILE)
		@echo "-> linking $@"
		$(LINKCXX) $< $(LINKCXXSCIPALL) \
		$(LINKCXX_l)$(TPILIB) \
		$(addprefix $(LINKCXX_L),$(CR_LIBDIR)) \
		$(LINKCXX_l)$(CR_LIB) \
		$(LINKRPATH)$(realpath $(SCIPDIR)/lib/shared) \
		$(addprefix $(LINKRPATH),$(realpath $(CR_LIBDIR))) \
		$(SANITZERFLAGS) \
		$(LINKCXX_o)$@

$(LPILIBFILE):
		$(error $@ not available. SCIP must be compiled with OPT=dbg SHARED=true IPOPT=$(IPOPT) LPS=$(LPS))
$(SCIPLIBFILE):
		$(error $@ not available. SCIP must be compiled with OPT=dbg SHARED=true IPOPT=$(IPOPT) LPS=$(LPS))
$(NLPILIBFILE):
		$(error $@ not available. SCIP must be compiled with OPT=dbg SHARED=true IPOPT=$(IPOPT) LPS=$(LPS))
$(TPILIBFILE):
		$(error $@ not available. SCIP must be compiled with OPT=dbg SHARED=true IPOPT=$(IPOPT) LPS=$(LPS) TPI=$(TPI))


$(OBJDIR)/%.o : $(SRCDIR)/%.c | $(OBJDIR) $(OBJDIRS) $(BINDIRS)
		@echo "-> compiling $@"
		$(CC) $(addprefix -I,$(CR_INC)) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

# TODO: checkout v2.3.0 when it becomes available
Criterion:
	@echo -n "Criterion not found, downloading and building it in"
	@for (( count=3; $$count ; count-- )) do echo -n " $$count" ; sleep 1 ; done ; echo
	git clone --depth 1 -b v2.2.2 https://github.com/Snaipe/Criterion
	cd Criterion; \
	mkdir build; \
	cd build; \
	cmake ..; \
	make
#---- EOF --------------------------------------------------------------------
