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
#@brief   Makefile for SCIP callable library example project
#@author  Stefan Vigerske


#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

SCIPDIR		=	../..


#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAMES	=	string gastrans circle brachistochrone circpacking
MAINOBJ		=	$(addsuffix .o,$(MAINNAMES))
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.c))
MAINDEP		=	$(SRCDIR)/depend.$(OPT)

MAINS		=	$(addsuffix .$(BASE).$(LPS)$(EXEEXTENSION),$(MAINNAMES))
MAINFILES	=	$(addprefix $(BINDIR)/,$(MAINS))
MAINSHORTLINKS	=	$(addprefix $(BINDIR)/,$(MAINNAMES))
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILES) $(MAINOBJFILES) $(MAINSHORTLINKS)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILES) $(MAINSHORTLINKS)

.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		-ln -fs $(SCIPDIR)/lint
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) -I$(SCIPDIR) lint/main-gcc.lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

$(MAINSHORTLINKS):	$(MAINFILES)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $@.$(BASE).$(LPS)$(EXEEXTENSION)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: test
test: $(MAINFILES)
		for i in $(MAINFILES) ; do $$i || break 1 ; done

.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		-rm -f $(OBJDIR)/*.o
		-rmdir $(OBJDIR)
endif
		-rm -f $(MAINFILES)

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).cpp|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.cpp|g'\'' \
		>$(MAINDEP)'

-include	$(MAINDEP)

.PHONY: doc
doc:
		cd doc && doxygen scip.dxy

$(BINDIR)/%.$(BASE).$(LPS)$(EXEEXTENSION): $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(OBJDIR)/%.o | $(BINDIR)
		@echo "-> linking $@"
		$(LINKCXX) $(OBJDIR)/$*.o $(LINKCXXSCIPALL)	$(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c | $(OBJDIR)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

#---- EOF --------------------------------------------------------------------
