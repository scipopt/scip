# $Id: Makefile,v 1.43 2003/09/03 17:13:53 bzfpfend Exp $
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
# lib/cpxinc                         -> directory with cplex.h
# lib/spxinc                         -> directory with SoPlex *.h
# lib/libcplex.$(OSTYPE).$(ARCH).a   -> libcplex.a
# lib/libsoplex.$(BASE).a            -> itself
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
NAME		=	scip

#-----------------------------------------------------------------------------

LIBOBJ		=	branch.o \
			branch_fullstrong.o \
			branch_mostinf.o \
			branch_leastinf.o \
			buffer.o \
			clock.o \
			conflict.o \
			cons.o \
			cons_and.o \
			cons_binpack.o \
			cons_bitarith.o \
			cons_bitvar.o \
			cons_eqknapsack.o \
			cons_integral.o \
			cons_linear.o \
			cons_logicor.o \
			cons_invarknapsack.o \
			cons_knapsack.o \
			cons_setppc.o \
			cons_varlb.o \
			cons_varub.o \
			cutpool.o \
			disp.o \
			disp_default.o \
			event.o \
			heur.o \
			heur_diving.o \
			heur_rounding.o \
			interrupt.o \
			lp.o \
			mem.o \
			memory.o \
			message.o \
			misc.o \
			nodesel.o \
			nodesel_bfs.o \
			nodesel_dfs.o \
			nodesel_restartdfs.o \
			paramset.o \
			presol.o \
			presol_dualfix.o \
			presol_trivial.o \
			price.o \
			primal.o \
			prob.o \
			reader.o \
			reader_cnf.o \
			reader_mps.o \
			reader_rtp.o \
			retcode.o \
			scip.o \
			sepa.o \
			sepa_gomory.o \
			sepastore.o \
			set.o \
			sol.o \
			solve.o \
			stat.o \
			tree.o \
			var.o


#-----------------------------------------------------------------------------
include make/make.$(BASE)
#-----------------------------------------------------------------------------

LIBRARY	=	$(LIBDIR)/lib$(NAME).$(BASE).a
DEPLIB	=	src/depend.scip

LIBXXX	=	$(addprefix $(OBJDIR)/,$(LIBOBJ))
LIBSRC	=	$(addprefix $(SRCDIR)/,$(LIBOBJ:.o=.c))

#-----------------------------------------------------------------------------
# LP Solver CPLEX
#-----------------------------------------------------------------------------
ifeq ($(LPS),cpx)
	CPPFLAGS+=	-I$(LIBDIR)/cpxinc
	OBJECT	=	cmain.o
	LPSLIB	=	cplex.$(OSTYPE).$(ARCH)
	LPIOBJ	=	lpi_cpx.o bitencode.o
	LPILIB	=	$(LIBDIR)/libcpxlp.$(BASE).a
	LPIXXX	=	$(addprefix $(OBJDIR)/,$(LPIOBJ))
	LPISRC  =	$(addprefix $(SRCDIR)/,$(LPIOBJ:.o=.c))
	TARGET	=	$(NAME).$(BASE).cpx
	BINARY	=	$(BINDIR)/$(TARGET)
	DEPLPS	=	src/depend.cpx
	OBJXXX	=	$(addprefix $(OBJDIR)/,$(OBJECT))
	OBJSRC	=	$(addprefix $(SRCDIR)/,$(OBJECT:.o=.c))

$(BINARY):	$(OBJDIR) $(BINDIR) $(OBJXXX) $(LIBRARY) $(LPILIB) 
		$(CXX) $(CXXFLAGS) $(OBJXXX) \
		-L$(LIBDIR) -l$(NAME).$(BASE) -lcpxlp.$(BASE) \
		-l$(LPSLIB) $(LDFLAGS) -o $@

endif
#-----------------------------------------------------------------------------
# LP Solver SoPlex
#-----------------------------------------------------------------------------
ifeq ($(LPS),spx)
	CPPFLAGS+=	-I$(LIBDIR)/spxinc 
	OBJECT	=	cppmain.o
	LPSLIB	=	soplex.$(OSTYPE).$(ARCH).$(COMP).$(OPT)
	LPIOBJ	=	lpi_spx.o bitencode.o
	LPILIB	=	$(LIBDIR)/libspxlp.$(BASE).a
	LPIXXX	=	$(addprefix $(OBJDIR)/,$(LPIOBJ))
	LPISRC  =	$(addprefix $(SRCDIR)/,$(LPIOBJ:.o=.cpp))
	TARGET	=	$(NAME).$(BASE).spx
	BINARY	=	$(BINDIR)/$(TARGET)
	DEPLPS	=	src/depend.spx
	OBJXXX	=	$(addprefix $(OBJDIR)/,$(OBJECT))
	OBJSRC	=	$(addprefix $(SRCDIR)/,$(OBJECT:.o=.cpp))

$(BINARY):	$(OBJDIR) $(BINDIR) $(OBJXXX) $(LIBRARY) $(LPILIB)
		$(CXX) $(CXXFLAGS) $(OBJXXX) \
		-L$(LIBDIR) -l$(NAME).$(BASE) -lspxlp.$(BASE) \
		-l$(LPSLIB) $(LDFLAGS) -o $@

endif

#-----------------------------------------------------------------------------
# SCIP Library
#-----------------------------------------------------------------------------
$(LIBRARY):	$(OBJDIR) $(LIBDIR) $(LIBXXX) 
		-rm -f $@
		$(AR) $(ARFLAGS) $@ $(LIBXXX) 
		$(RANLIB) $@

$(LPILIB):	$(OBJDIR) $(LIBDIR) $(LPIXXX)
		-rm -f $@
		$(AR) $(ARFLAGS) $@ $(LPIXXX)
		$(RANLIB) $@

#-----------------------------------------------------------------------------

lint:		$(LIBSRC) $(OBJSRC)
		$(LINT) lint/scip.lnt -os\(lint.out\) \
		$(CPPFLAGS) -UNDEBUG $^

doc:		
		cd doc; $(DOXY) scip.dxy

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
		-rm -rf $(OBJDIR)/* $(LIBRARY) $(TARGET)

depend:
		$(SHELL) -ec '$(DCC) $(DFLAGS) $(CPPFLAGS) $(LIBSRC) \
		| sed '\''s|^\([0-9A-z\_]\{1,\}\)\.o|$$\(OBJDIR\)/\1.o|g'\'' \
		>$(DEPLIB)'
ifeq ($(LPS),cpx)
		$(SHELL) -ec '$(DCC) $(DFLAGS) $(CPPFLAGS) $(OBJSRC) $(LPISRC) \
		| sed '\''s|^\([0-9A-z\_]\{1,\}\)\.o|$$\(OBJDIR\)/\1.o|g'\'' \
		>$(DEPLPS)'
endif
ifeq ($(LPS),spx)
		$(SHELL) -ec '$(DCXX) $(DFLAGS) $(CPPFLAGS) $(OBJSRC) $(LPISRC) \
		| sed '\''s|^\([0-9A-z\_]\{1,\}\)\.o|$$\(OBJDIR\)/\1.o|g'\'' \
		>$(DEPLPS)'
endif

-include	$(DEPLIB)
-include 	$(DEPLPS)

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# --- EOF ---------------------------------------------------------------------
