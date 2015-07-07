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

VERBOSE 	= 	false
MAKELEMONLINKS  = 	true

FLAGS		+=	-I$(LIBDIR)/lemon/include
LEMONLINKS	=	$(LIBDIR)/lemon
LINKSMARKERFILE	=	$(LIBDIR)/linkscreated
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
		@for i in $(TESTFILES) ; do $(MAINFILE) $$i || exit 1 ; done

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
		$(LINKCXX_L)$(LIBDIR)/lemon/lib $(LINKCXX_l)emon


$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -Wno-shadow -c $< $(CXX_o)$@

info:
	@echo "$(MAINOBJFILES)"

.PHONY: removemarker
removemarker:
		@rm -f $(LINKSMARKERFILE)

.PHONY: links
links:		|$(LIBDIR) removemarker echosoftlinks $(LEMONLINKS)

.PHONY: echosoftlinks
echosoftlinks:
		@echo
		@echo "- Current settings: LPS=$(LPS) OSTYPE=$(OSTYPE) ARCH=$(ARCH) COMP=$(COMP) SUFFIX=$(LINKLIBSUFFIX) ZIMPL=$(ZIMPL) ZIMPLOPT=$(ZIMPLOPT) IPOPT=$(IPOPT) IPOPTOPT=$(IPOPTOPT) EXPRINT=$(EXPRINT) GAMS=$(GAMS)"
		@echo
		@echo "*   This example code for multiobjective optimization requires an installation of the lemon";
		@echo "*   graph library:";
		@echo
		@echo "*   'lib/lemon' should point to the main directory of your lemon installation";
		@echo "*   containing the subdirectories 'include' and 'lib', i.e.";
		@echo "*   'lib/lemon/include/lemon/list_graph.h' should exist";
		@echo "*   The links will be installed in the 'lib' directory.";
		@echo "*   For more information and if you experience problems see the INSTALL file.";

.PHONY: $(LEMONLINKS)
$(LEMONLINKS): $(LIBDIR)
ifeq ($(MAKELEMONLINKS), true)
		@$(SHELL) -ec 'if test ! -e $(LINKSMARKERFILE); \
			then \
				echo ; \
				echo "- preparing missing soft-link \"$@\":" ; \
				DIRNAME=`dirname $@` ; \
				# test for a systemwide lemon library - grep returns 0 if it succesfully finds a match \
				if test "`ldconfig -p | grep libemon`" == "0" ; \
				then \
					echo "> Enter soft-link target file or directory for \"$@\" (not required because lemon installation was found in system): " ; \
				else \
					echo "> Enter soft-link target file or directory for \"$@\" (required because no lemon installation was found in system): " ; \
				fi ; \
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
		@rm -f $(LINKSMARKERFILE)
		@echo "this is only a marker" > $(LINKSMARKERFILE)
endif

#---- EOF --------------------------------------------------------------------
