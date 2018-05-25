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
#@brief   Makefile for GAMS interface project
#@author  Stefan Vigerske


#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

SCIPDIR		=	../..
GAMSDIR		=	gams

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------

include $(SCIPDIR)/make/make.project

ifneq ($(SHARED), true)
$(error Please compile SCIP and this interface with SHARED=true)
endif

ifneq ($(READLINE), false)
$(error Please compile SCIP and this interface with READLINE=false)
endif

GAMSIODIR	=	$(GAMSDIR)/apifiles/C/api

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	gamsscip
# C sources
MAINCOBJ	=	reader_gmo.o event_solvetrace.o GamsSolveTrace.o
MAINCSRC	=	$(addprefix $(SRCDIR)/,$(MAINCOBJ:.o=.c))
# C++ sources
MAINCXXOBJ	=	GamsScip.o GamsSolver.o GamsOptions.o
MAINCXXSRC	=	$(addprefix $(SRCDIR)/,$(MAINCXXOBJ:.o=.cpp))
# C sources in GAMS I/O directory
MAINGOBJ	=	gdxcc.o gevmcc.o gmomcc.o optcc.o
MAINGSRC	=	$(addprefix $(GAMSIODIR)/,$(MAINGOBJ:.o=.c))
MAINDEP		=	$(SRCDIR)/depend.$(OPT)

MAIN		=	lib$(MAINNAME).$(BASE).$(LPS).$(LIBEXT)
MAINFILE	=	$(LIBDIR)/$(MAIN)
MAINSHORTLINK	=	$(LIBDIR)/lib$(MAINNAME).$(LIBEXT)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINCOBJ) $(MAINCXXOBJ) $(MAINGOBJ))

#-----------------------------------------------------------------------------
# GMO library
#-----------------------------------------------------------------------------

# where to find header files for GAMS I/O files
FLAGS += -I$(GAMSIODIR)

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: install
install:
		@./gamsinst.sh "$(GAMSDIR)"

.PHONY: lint
lint:		$(MAINCSRC) $(MAINCXXSRC)
		-rm -f lint.out
		-ln -fs $(SCIPDIR)/lint
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) lint/scip.lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):	
		@-mkdir -p $(OBJDIR)

$(LIBDIR):
		@-mkdir -p $(LIBDIR)

.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		-rm -f $(OBJDIR)/*.o
		-rmdir $(OBJDIR)
endif
		-rm -f $(MAINFILE)

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(MAINCSRC) $(MAINCXXSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).cpp|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.cpp|g'\'' \
		| sed '\''s|$(GAMSIODIR).*h||g'\'' \
		>$(MAINDEP)'

.PHONY: test
test:
		./test.sh $(MAINFILE)

-include	$(MAINDEP)

$(MAINFILE):	$(LIBDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ \
		$(MAINOBJFILES) $(LINKCXXSCIPALL)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif


$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

$(OBJDIR)/%.o:	$(GAMSIODIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS:-pedantic=) -Wno-shadow -c $< $(CC_o)$@

#---- EOF --------------------------------------------------------------------
