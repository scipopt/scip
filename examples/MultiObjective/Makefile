#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic Licence.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#@file    Makefile
#@brief   Makefile for Özpeynirci-Köksalan Algorithm using SCIP as a callable library
#@author  Timo Strunk

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

SCIPDIR		=	../../

SOLDIR		=	solutions

# define to be able to locate library files
ifeq ($(OSTYPE),mingw)
MULTOBJDIR		=	./
else
MULTOBJDIR		=	$(realpath .)
endif

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------

include $(SCIPDIR)/make/make.project

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	multiopt
MAINOBJ		=	main.o
SECOBJS		= 	reader_mop.o WeightedSolver.o LiftedWeightSpaceSolver.o\
			Skeleton.o WeightSpaceVertex.o Objectives.o

MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINDEP		=	$(SRCDIR)/depend.cppmain.$(OPT)

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ) $(SECOBJS) )

TESTSETDIR	=	data
TESTFILES       =       data/*.mop
#INCLUDE 	= 	-Ilib/lemon-src

VERBOSE 	= 	false
MAKELEMONLINKS  = 	true

FLAGS		+=	-I$(LIBDIR)/lemoninc
LEMONLINKS	=	$(LIBDIR)/lemoninc
LEMONLINKS	+=	$(LIBDIR)/liblemon.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
LEMONLINKS	+=	$(LIBDIR)/liblemon.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
MAKE		+= -s
endif

.PHONY: all
all:            $(LEMONLINKS) $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) lint/$(MAINNAME).lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

.PHONY: scip
scip:
		@$(MAKE) -C $(SCIPDIR) libs $^

.PHONY: doc
doc:
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/scip.css);
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/pictures/scippy.png);
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/pictures/miniscippy.png);
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/scipfooter.html footer.html);
		cd doc; $(DOXY) $(MAINNAME).dxy

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

$(SOLDIR):
		@-mkdir -p $(SOLDIR)

$(LIBDIR):
		@-mkdir -p $(LIBDIR)

$(TESTSETDIR):
		@-mkdir -p check/$(TESTSETDIR)

.PHONY: test
test:           $(MAINFILE) $(TESTSETDIR)
		@for i in $(TESTFILES) ; do $(MAINFILE) $$i || break 1 ; done

.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		-rm -f $(OBJDIR)/*.o
		-rmdir $(OBJDIR)
endif
		-rm -f $(MAINFILE)

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).cpp|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.cpp|g'\'' \
		>$(MAINDEP)'

-include	$(MAINDEP)

$(MAINFILE):	$(BINDIR) $(OBJDIR) $(LIBDIR) $(MAINOBJFILES)
		@$(MAKE) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES) $(OBJS) $(EXTOBS) $(SOLDIR)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) $(OBJS) $(EXTOBJS) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
		$(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) 	\
		$(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) -Wno-shadow $(LINKCXX_o) $@ \
		$(LINKCXX_L)$(LIBDIR) $(LINKCXX_l)lemon.$(OSTYPE).$(ARCH).$(COMP)


$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -Wno-shadow -c $< $(CXX_o)$@

info:
	@echo "$(MAINOBJFILES)"

.PHONY: $(LEMONLINKS)
ifeq ($(MAKELEMONLINKS), true)
		@$(SHELL) -ec 'if test ! -e $@ ; \
			then \
				DIRNAME=`dirname $@` ; \
				BASENAMEA=`basename $@ .$(STATICLIBEXT)` ; \
				BASENAMESO=`basename $@ .$(SHAREDLIBEXT)` ; \
				echo ; \
				echo "- preparing missing soft-link \"$@\":" ; \
				if test -e $$DIRNAME/$$BASENAMEA.$(SHAREDLIBEXT) ; \
				then \
					echo "* this soft-link is not necessarily needed since \"$$DIRNAME/$$BASENAMEA.$(SHAREDLIBEXT)\" already exists - press return to skip" ; \
				fi ; \
				if test -e $$DIRNAME/$$BASENAMESO.$(STATICLIBEXT) ; \
				then \
					echo "* this soft-link is not necessarily needed since \"$$DIRNAME/$$BASENAMESO.$(STATICLIBEXT)\" already exists - press return to skip" ; \
				fi ; \
				echo "> Enter soft-link target file or directory for \"$@\" (return if not needed): " ; \
				echo -n "> " ; \
				cd $$DIRNAME ; \
				eval $(READ) TARGET ; \
				cd - ; \
				if test "$$TARGET" != "" ; \
				then \
					echo "-> creating softlink \"$@\" -> \"$$TARGET\"" ; \
					rm -f $@ ; \
					$(LN_s) $$TARGET $@ ; \
				else \
					echo "* skipped creation of softlink \"$@\". Call \"make links\" if needed later." ; \
				fi ; \
				echo ; \
			fi'
endif

#---- EOF --------------------------------------------------------------------
