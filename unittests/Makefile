#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            *
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
UNITTESTS	=  unittest-boundchg unittest-setters unittest-relax unittest-cons unittest-presol unittest-depthlevel unittest-slack unittest-lpibases unittest-regression
EXTRASRC	=  unittest-relax/relax_unittest unittest-cons/cons_unittest unittest-presol/presol_unittest


# UNITTESTS, EXECCOMMANDS, and TESTCASES depend on each other
# command which should be used to start the unittest, e.g. optimize or presolve
EXECCOMMANDS    =       NOEXEC NOEXEC NOEXEC NOEXEC NOEXEC NOEXEC testslack NOEXEC NOEXEC
# testsets or NOTEST if there is no testset
TESTCASES       =       nofile nofile nofile nofile nofile nofile short nofile nofile

UNITTESTSSRC	=	$(addprefix $(SRCDIR)/, $(join $(addsuffix /,$(UNITTESTS)), $(addsuffix .c, $(UNITTESTS))))
UNITTESTSOBJ	=	$(addprefix $(OBJDIR)/, $(addsuffix .o, $(UNITTESTS)))
UNITTESTSOBJ	+=	$(addprefix $(OBJDIR)/, $(notdir $(addsuffix .o, $(EXTRASRC))))
UNITTESTSBIN	=	$(addprefix $(BINDIR)/, $(addsuffix $(BASESUFFIX), $(UNITTESTS)))
UNITTESTSEXECS	=	$(addsuffix $(BASESUFFIX), $(UNITTESTS))
UNITTESTSDEP	=	$(SRCDIR)/depend.$(OPT)

# number of test cases in the whole unittest
NUMBERS        :=       $(shell seq 1 $(words $(UNITTESTSEXECS)))

# where files should be looked in
VPATH = $(addprefix $(SRCDIR)/, $(UNITTESTS))

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(UNITTESTSBIN) $(UNITTESTSOBJ)
endif

.PHONY: all
all:            $(SCIPDIR) $(UNITTESTSBIN)

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

# if $(BINDIR) does *not* exist: create it and generate links
t_bin:
		-@test -d $(BINDIR) || { \
		echo "-> Creating $(BINDIR) directory"; \
		mkdir -p $(BINDIR); }

.PHONY: test
test:           $(UNITTESTSBIN)
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/getlastprob.awk);

		@cd check; \
		$(foreach var,$(NUMBERS), \
			echo; \
			echo "--------------------------------------------------------------------"; \
			echo "running test $(var) with executable $(word $(var),$(UNITTESTSBIN))"; \
			$(SHELL) -ec './check.sh $(word $(var),$(TESTCASES)) $(word $(var),$(UNITTESTSBIN)) $(SETTINGS) $(word $(var),$(UNITTESTSEXECS)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(VALGRIND) $(word $(var),$(EXECCOMMANDS))')


.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		@-(rm -f $(OBJDIR)/*.o && rmdir $(OBJDIR));
		@echo "-> remove main objective files"
endif
		@-rm -f $(UNITTESTSBIN)
		@echo "-> remove binaries"

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*.c src/*.h $(SCIPDIR)/src/scip/*.c $(SCIPDIR)/src/scip/*.h;

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(UNITTESTSSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(UNITTESTSDEP)'

-include	$(UNITTESTSDEP)

$(UNITTESTSBIN): $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(UNITTESTSOBJ) | t_bin
		@echo "-> linking $@"
		$(LINKCXX) $(OBJDIR)/$(notdir $(patsubst %$(BASESUFFIX),%,$@)).o \
		$(addprefix $(OBJDIR)/, $(addsuffix .o, $(notdir $(filter $(notdir $(patsubst %$(BASESUFFIX),%,$@))/%, $(EXTRASRC))))) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@


$(OBJDIR)/%.o:	%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
