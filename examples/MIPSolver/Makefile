#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2004 Tobias Achterberg                              *
#*                                                                           *
#*                  2002-2004 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the SCIP Academic Licence.        *
#*                                                                           *
#*  You should have received a copy of the SCIP Academic License             *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: Makefile,v 1.2 2004/11/15 10:55:43 bzfpfend Exp $

#@file    Makefile
#@brief   Makefile for C++ example project using SCIP as a callable library
#@author  Thorsten Koch
#@author  Tobias Achterberg

.PHONY:		depend clean lint doc test

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
COMP		=	gnu
TIME     	=  	3600
NODES           =       2100000000
MEM		=	1024
TEST		=	
SETTINGS        =       default

LINKER  	=	CPP

CC		=	gcc
CXX		=	g++
DCC		=	gcc
DCXX		=	g++
AR		=	ar
RANLIB		=	ranlib
LINT		=	flexelint
DOXY		=	doxygen

FLAGS		=	-I$(SRCDIR) -I$(LIBDIR)/scip/src
CFLAGS		=	
CXXFLAGS	=	
LDFLAGS		=	-lpthread -lm -lz
ARFLAGS		=	cr
DFLAGS		=	-MM

GCCWARN		=	-Wall -W -Wpointer-arith -Wcast-align -Wwrite-strings \
			-Wstrict-prototypes -Wmissing-prototypes \
			-Wmissing-declarations -Wno-unknown-pragmas \
			-Wno-unused

GXXWARN		=	-Wall -W -Wpointer-arith \
			-Wcast-align -Wwrite-strings -Wconversion \
			-Wstrict-prototypes -Wmissing-prototypes \
			-Wno-unknown-pragmas \
			-Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder \
			-Woverloaded-virtual -Wsign-promo -Wsynth -Wundef \
			-Wcast-qual -Wold-style-cast -Wshadow 

BASE		=	$(OSTYPE).$(ARCH).$(COMP).$(OPT)
OBJDIR		=	obj/O.$(BASE)
SRCDIR		=	src
BINDIR		=	bin
LIBDIR		=	lib


#-----------------------------------------------------------------------------
include make/make.$(BASE)
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Memory Management
#-----------------------------------------------------------------------------

#FLAGS		+=	-DNOSAFEMEM
#FLAGS		+=	-DNOBLOCKMEM


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	example
MAINOBJ		=	cppmain.o

MAIN		=	$(MAINNAME).$(BASE).$(LPS)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINXXX		=	$(addprefix $(OBJDIR)/,$(MAINOBJ))
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINDEP		=	src/depend.main


#-----------------------------------------------------------------------------
# LP Solver Interface
#-----------------------------------------------------------------------------

LPILIBNAME	=	lpi$(LPS)

ifeq ($(LPS),cpx)
LPSLIB		=	cplex.$(OSTYPE).$(ARCH).$(COMP)
endif

ifeq ($(LPS),spx)
LINKER		=	CPP
LPSLIB		=	soplex.$(OSTYPE).$(ARCH).$(COMP)
endif

LPILIB		=	$(LPILIBNAME).$(BASE)
LPILIBFILE	=	$(LIBDIR)/scip/lib/lib$(LPILIB).a


#-----------------------------------------------------------------------------
# SCIP Library
#-----------------------------------------------------------------------------

SCIPLIBNAME	=	scip
SCIPLIB		=	$(SCIPLIBNAME).$(BASE)
SCIPLIBFILE	=	$(LIBDIR)/scip/lib/lib$(SCIPLIB).a


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

all:            $(MAINFILE)

lint:		$(MAINSRC)
		$(LINT) lint/$(MAINNAME).lnt -os\(lint.out\) \
		$(FLAGS) -UNDEBUG $^

doc:		
		cd doc; $(DOXY) $(MAINNAME).dxy

test:		
		cd check; \
		/bin/sh ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(MAIN).$(MACHINENAME) $(TIME) $(NODES) $(MEM);

$(OBJDIR):	
		-mkdir -p $(OBJDIR)

$(LIBDIR):
		-mkdir -p $(LIBDIR)

$(BINDIR):
		-mkdir -p $(BINDIR)

clean:
		-rm -rf $(OBJDIR)/* $(MAINFILE)

depend:
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-z\_]\{1,\}\)\.o|$$\(OBJDIR\)/\1.o|g'\'' \
		>$(MAINDEP)'

-include	$(MAINDEP)

$(MAINFILE):	$(OBJDIR) $(BINDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(MAINXXX)
ifeq ($(LINKER),C)
		$(CC) $(CFLAGS) $(MAINXXX) \
		-L$(LIBDIR)/scip/lib -l$(SCIPLIB) -l$(LPILIB) -l$(LPSLIB) -L$(LIBDIR) \
		$(LDFLAGS) -o $@
endif
ifeq ($(LINKER),CPP)
		$(CXX) $(CXXFLAGS) $(MAINXXX) \
		-L$(LIBDIR)/scip/lib -l$(SCIPLIB) -l$(LPILIB) -l$(LPSLIB) -L$(LIBDIR) \
		$(LDFLAGS) -o $@
endif

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		$(CC) $(FLAGS) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		$(CXX) $(FLAGS) $(CXXFLAGS) -c $< -o $@

# --- EOF ---------------------------------------------------------------------
