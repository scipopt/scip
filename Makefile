# $Id: Makefile,v 1.46 2003/10/15 13:09:25 bzfpfend Exp $
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2002 Tobias Achterberg                              *
#*                            Thorsten Koch                                  *
#*                  2002-2002 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the SCIP Academic Licence.        *
#*                                                                           *
#*  You should have received a copy of the SCIP Academic License             *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#@file    Makefile
#@brief   SCIP Makefile
#@author  Thorsten Koch
#@author  Tobias Achterberg

# Several symlinks are nesseccary:
#
# lib/cpxinc                                        -> directory with cplex.h
# lib/spxinc                                        -> directory with SoPlex's *.h
# lib/libcplex.$(OSTYPE).$(ARCH).a                  -> libcplex.a
# lib/libsoplex.$(OSTYPE).$(ARCH).$(COMP).$(OPT).a  -> libsoplex.a
#
#
.PHONY:		depend clean lint doc test testcplex

ARCH            :=      $(shell uname -m | \
                        sed \
			-e s/sun../sparc/ \
			-e s/i.86/x86/ \
			-e s/IP../mips/ \
			-e s/9000..../hppa/ \
			-e s/00........../pwr4/)
OSTYPE		:=	$(shell uname -s | tr A-Z a-z)


#-----------------------------------------------------------------------------
# default settings
#-----------------------------------------------------------------------------

OPT		=	dbg
LPS		=	cpx
PAR		=	seq

COMP		=	gnu
CC		=	gcc
CXX		=	g++
DCC		=	gcc
DCXX		=	g++
LINT		=	flexelint
AR		=	ar
RANLIB		=	ranlib
DOXY		=	doxygen

SETTINGS        =       hybridhistory

SRCDIR		=	src
BINDIR		=	bin
LIBDIR		=	lib

CPPFLAGS	=	-I$(SRCDIR)
CFLAGS		=	
CXXFLAGS	=	
LDFLAGS		=	-lpthread -lm
ARFLAGS		=	cr
DFLAGS		=	-MM

GCCWARN		=	-Wall -W -Wpointer-arith -Wcast-align -Wwrite-strings \
			-Wstrict-prototypes -Wmissing-prototypes \
			-Wmissing-declarations -Wno-unknown-pragmas \
			-Wno-unused

GXXWARN		=	-Wall -W -Wpointer-arith -Wbad-function-cast \
			-Wcast-align -Wwrite-strings -Wconversion \
			-Wstrict-prototypes -Wmissing-prototypes \
			-Wmissing-declarations -Wno-unknown-pragmas \
			-Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder \
			-Woverloaded-virtual -Wsign-promo -Wsynth -Wundef \
			-Wcast-qual -Wold-style-cast -Wshadow 

BASE		=	$(PAR).$(OSTYPE).$(ARCH).$(COMP).$(OPT)
OBJDIR		=	obj/O.$(BASE)


#-----------------------------------------------------------------------------
include make/make.$(BASE)
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Program Components
#-----------------------------------------------------------------------------

NAME		=	scip
OBJECTS		=

OBJXXX		=	$(addprefix $(OBJDIR)/,$(OBJECTS))
OBJSRC		=	$(addprefix $(SRCDIR)/,$(OBJECTS:.o=.c))
OBJDEP		=	src/depend.obj


#-----------------------------------------------------------------------------
# LP Solver Interface
#-----------------------------------------------------------------------------

ifeq ($(LPS),cpx)
CPPFLAGS	+=	-I$(LIBDIR)/cpxinc
MAINOBJ		=	cmain.o
LPSLIB		=	cplex.$(OSTYPE).$(ARCH)
LPIOBJ		=	lpi_cpx.o bitencode.o
endif

ifeq ($(LPS),spx)
CPPFLAGS	+=	-I$(LIBDIR)/spxinc 
MAINOBJ		=	cppmain.o
LPSLIB		=	soplex.$(OSTYPE).$(ARCH).$(COMP).$(OPT)
LPIOBJ		=	lpi_spx.o bitencode.o
endif

LPILIB		=	$(LPS)lp.$(BASE)
LPILIBFILE	=	$(LIBDIR)/lib$(LPILIB).a
LPIXXX		=	$(addprefix $(OBJDIR)/,$(LPIOBJ))
LPISRC  	=	$(addprefix $(SRCDIR)/,$(LPIOBJ:.o=.c))
LPIDEP		=	src/depend.$(LPS)
MAINXXX		=	$(addprefix $(OBJDIR)/,$(MAINOBJ))
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.c))


#-----------------------------------------------------------------------------
# SCIP Library
#-----------------------------------------------------------------------------

SCIPOBJ		=	branch.o \
			buffer.o \
			clock.o \
			conflict.o \
			cons.o \
			cutpool.o \
			disp.o \
			event.o \
			heur.o \
			interrupt.o \
			lp.o \
			mem.o \
			memory.o \
			message.o \
			misc.o \
			nodesel.o \
			paramset.o \
			presol.o \
			price.o \
			primal.o \
			prob.o \
			reader.o \
			retcode.o \
			scip.o \
			scipdefplugins.o \
			sepa.o \
			sepastore.o \
			set.o \
			sol.o \
			solve.o \
			stat.o \
			tree.o \
			var.o \
			branch_fullstrong.o \
			branch_mostinf.o \
			branch_leastinf.o \
			cons_and.o \
			cons_binpack.o \
			cons_eqknapsack.o \
			cons_invarknapsack.o \
			cons_integral.o \
			cons_knapsack.o \
			cons_linear.o \
			cons_logicor.o \
			cons_setppc.o \
			cons_varlb.o \
			cons_varub.o \
			disp_default.o \
			heur_diving.o \
			heur_rounding.o \
			nodesel_bfs.o \
			nodesel_dfs.o \
			nodesel_restartdfs.o \
			presol_dualfix.o \
			presol_trivial.o \
			reader_cnf.o \
			reader_mps.o \
			sepa_gomory.o

SCIPLIB		=	$(NAME).$(BASE)
SCIPLIBFILE	=	$(LIBDIR)/lib$(SCIPLIB).a
SCIPDEP		=	src/depend.scip
SCIPXXX		=	$(addprefix $(OBJDIR)/,$(SCIPOBJ))
SCIPSRC		=	$(addprefix $(SRCDIR)/,$(SCIPOBJ:.o=.c))


#-----------------------------------------------------------------------------
# Target
#-----------------------------------------------------------------------------

TARGET		=	$(NAME).$(BASE).$(LPS)
BINARY		=	$(BINDIR)/$(TARGET)


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

$(BINARY):	$(OBJDIR) $(BINDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(OBJXXX) $(MAINXXX)
		$(CXX) $(CXXFLAGS) $(OBJXXX) $(MAINXXX) \
		-L$(LIBDIR) -l$(SCIPLIB) -l$(LPILIB) -l$(LPSLIB) \
		$(LDFLAGS) -o $@

lint:		$(SCIPSRC) $(OBJSRC) $(MAINSRC)
		$(LINT) lint/$(NAME).lnt -os\(lint.out\) \
		$(CPPFLAGS) -UNDEBUG $^

doc:		
		cd doc; $(DOXY) $(NAME).dxy

tmpcheck:
		cd check; \
		/bin/sh ./tmpcheck.sh $(TEST).test $(BINARY);

test:		
		cd check; \
		/bin/sh ./check.sh $(TEST).test $(BINARY) $(SETTINGS);

testcplex:		
		cd check; \
		/bin/sh ./check_cplex.sh $(TEST).test $(OSTYPE).$(ARCH);

$(OBJDIR):	
		-mkdir -p $(OBJDIR)

$(LIBDIR):
		-mkdir -p $(LIBDIR)

$(BINDIR):
		-mkdir -p $(BINDIR)

clean:
		-rm -rf $(OBJDIR)/* $(SCIPLIBFILE) $(LPILIBFILE) $(TARGET)

depend:
		$(SHELL) -ec '$(DCC) $(DFLAGS) $(CPPFLAGS) $(SCIPSRC) \
		| sed '\''s|^\([0-9A-z\_]\{1,\}\)\.o|$$\(OBJDIR\)/\1.o|g'\'' \
		>$(SCIPDEP)'
ifeq ($(LPS),cpx)
		$(SHELL) -ec '$(DCC) $(DFLAGS) $(CPPFLAGS) $(MAINSRC) $(LPISRC) \
		| sed '\''s|^\([0-9A-z\_]\{1,\}\)\.o|$$\(OBJDIR\)/\1.o|g'\'' \
		>$(LPIDEP)'
endif
ifeq ($(LPS),spx)
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(CPPFLAGS) $(MAINSRC) $(LPISRC) \
		| sed '\''s|^\([0-9A-z\_]\{1,\}\)\.o|$$\(OBJDIR\)/\1.o|g'\'' \
		>$(LPIDEP)'
endif
		$(SHELL) -ec '$(DCC) $(DFLAGS) $(CPPFLAGS) $(OBJSRC) \
		| sed '\''s|^\([0-9A-z\_]\{1,\}\)\.o|$$\(OBJDIR\)/\1.o|g'\'' \
		>$(OBJDEP)'

-include	$(SCIPDEP)
-include 	$(LPIDEP)
-include	$(OBJDEP)

$(SCIPLIBFILE):	$(OBJDIR) $(LIBDIR) $(SCIPXXX) 
		-rm -f $@
		$(AR) $(ARFLAGS) $@ $(SCIPXXX) 
		$(RANLIB) $@

$(LPILIBFILE):	$(OBJDIR) $(LIBDIR) $(LPIXXX)
		-rm -f $@
		$(AR) $(ARFLAGS) $@ $(LPIXXX)
		$(RANLIB) $@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# --- EOF ---------------------------------------------------------------------
