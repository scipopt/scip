#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2002 Tobias Achterberg                              *
#*                                                                           *
#*                  2002-2002 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the SCIP Academic Licence.        *
#*                                                                           *
#*  You should have received a copy of the SCIP Academic License             *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: Makefile,v 1.53 2003/11/26 16:08:57 bzfpfend Exp $

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
MACHINENAME	:=	$(shell uname -n | tr A-Z a-z)


#-----------------------------------------------------------------------------
# default settings
#-----------------------------------------------------------------------------

OPT		=	dbg
LPS		=	cpx
TIME     	=  	3600
TEST		=	miplib

COMP		=	gnu
CC		=	gcc
CXX		=	g++
DCC		=	gcc
DCXX		=	g++
LINT		=	flexelint
AR		=	ar
RANLIB		=	ranlib
DOXY		=	doxygen

SETTINGS        =       default

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

BASE		=	$(OSTYPE).$(ARCH).$(COMP).$(OPT)
OBJDIR		=	obj/O.$(BASE)


#-----------------------------------------------------------------------------
include make/make.$(BASE)
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Memory Management
#-----------------------------------------------------------------------------

#CFLAGS		+=	-DNOSAFEMEM
#CFLAGS		+=	-DNOBLOCKMEM


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
LPSLIB		=	cplex.$(OSTYPE).$(ARCH)
LPIOBJ		=	lpi_cpx.o bitencode.o
LPISRC  	=	$(addprefix $(SRCDIR)/,$(LPIOBJ:.o=.c))
endif

ifeq ($(LPS),spx)
CPPFLAGS	+=	-I$(LIBDIR)/spxinc 
LPSLIB		=	soplex.$(OSTYPE).$(ARCH)
LPIOBJ		=	lpi_spx.o bitencode.o
LPISRC  	=	lpi_spx.cpp bitencode.c
endif

ifeq ($(LPS),spxdbg)
CPPFLAGS	+=	-I$(LIBDIR)/spxinc 
LPSLIB		=	soplexdbg.$(OSTYPE).$(ARCH)
LPIOBJ		=	lpi_spxdbg.o bitencode.o
LPISRC  	=	lpi_spxdbg.cpp bitencode.c
endif

LPILIB		=	$(LPS)lp.$(BASE)
LPILIBFILE	=	$(LIBDIR)/lib$(LPILIB).a
LPIXXX		=	$(addprefix $(OBJDIR)/,$(LPIOBJ))
LPIDEP		=	src/depend.$(LPS)
MAINOBJ		=	cmain.o
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
			dialog.o \
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
			pricestore.o \
			pricer.o \
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
			dialog_default.o \
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

test:		
		cd check; \
		/bin/sh ./check.sh $(TEST) $(BINARY) $(SETTINGS) $(TARGET).$(MACHINENAME) $(TIME);

testcplex:		
		cd check; \
		/bin/sh ./check_cplex.sh $(TEST) $(OSTYPE).$(ARCH).$(MACHINENAME) $(TIME);

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
