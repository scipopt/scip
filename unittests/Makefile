#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            *
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
#@author  Ambros Gleixner
#@author  Marc Pfetsch


#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

SCIPDIR         =       ..


#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------

include $(SCIPDIR)/make/make.project

BASESUFFIX	=	.$(BASE).$(LPS)$(EXEEXTENSION)


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

VERSION		=	unittests

#-----------------------------------------------------------------------------
# Unit test definitions
#-----------------------------------------------------------------------------

# add new unit tests here:
UNITTESTS	=	unittest-slack unittest-depthlevel unittest-boundchg

# UNITTESTS, EXECCOMMANDS, and TESTCASES depend on each other
# command which should be used to start the unittest, e.g. optimize or presolve
EXECCOMMANDS    =       testslack NOEXEC NOEXEC
# testsets or NOTEST if there is no testset
TESTCASES       =       short nofile nofile

UNITTESTSSRC	=	$(addprefix $(SRCDIR)/, $(addsuffix .c, $(UNITTESTS)))
UNITTESTSOBJ	=	$(addprefix $(OBJDIR)/, $(addsuffix .o, $(UNITTESTS)))
UNITTESTSEXECS	=	$(addsuffix $(BASESUFFIX), $(UNITTESTS))
UNITTESTSDEP	=	$(SRCDIR)/depend.$(OPT)

# number of test cases in the whole unittest
NUMBERS        :=       $(shell seq 1 $(words $(UNITTESTSEXECS)))

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(UNITTESTSEXECS) $(UNITTESTSOBJ)
endif

.PHONY: all
all:            $(SCIPDIR) $(UNITTESTSEXECS)

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
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: test
test:           $(UNITTESTSEXECS)
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/getlastprob.awk);

		@cd check; \
		$(foreach var,$(NUMBERS), \
			echo; \
			echo "--------------------------------------------------------------------"; \
			echo "running test $(var) with executable $(word $(var),$(UNITTESTSEXECS))"; \
			$(SHELL) -ec './check.sh $(word $(var),$(TESTCASES)) $(BINDIR)/$(word $(var),$(UNITTESTSEXECS)) $(SETTINGS) $(word $(var),$(UNITTESTSEXECS)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(VALGRIND) $(word $(var),$(EXECCOMMANDS))')


.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		@-(rm -f $(OBJDIR)/*.o && rmdir $(OBJDIR));
		@echo "-> remove main objective files"
endif
		@-rm -f $(BINDIR)/$(UNITTESTSEXECS)
		@echo "-> remove binary"

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*.c src/*.h $(SCIPDIR)/src/scip/*.c $(SCIPDIR)/src/scip/*.h;

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(UNITTESTSSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(UNITTESTSDEP)'

-include	$(UNITTESTSDEP)

$(UNITTESTSEXECS): $(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(UNITTESTSOBJ)
		@echo "-> linking $@"
		$(LINKCXX) $(OBJDIR)/$(notdir $(patsubst %$(BASESUFFIX),%,$@)).o \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$(BINDIR)/$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
