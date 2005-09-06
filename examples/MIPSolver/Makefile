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
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: Makefile,v 1.8 2005/09/06 17:28:59 bzfpfend Exp $

#@file    Makefile
#@brief   Makefile for C++ example project using SCIP as a callable library
#@author  Thorsten Koch
#@author  Tobias Achterberg

ARCH            :=      $(shell uname -m | \
                        sed \
			-e s/sun../sparc/ \
			-e s/i.86/x86/ \
			-e s/IP../mips/ \
			-e s/9000..../hppa/ \
			-e s/Power\ Macintosh/ppc/ \
			-e s/00........../pwr4/)
OSTYPE		:=	$(shell uname -s | tr '[:upper:]' '[:lower:]' | sed -e s/cygwin.*/cygwin/ -e s/irix../irix/ )
HOSTNAME	:=	$(shell uname -n | tr '[:upper:]' '[:lower:]')


#-----------------------------------------------------------------------------
# default settings
#-----------------------------------------------------------------------------

TIME     	=  	3600
NODES           =       2100000000
MEM		=	1024
TEST		=	
SETTINGS        =       default

VERBOSE		=	false
OPT		=	opt
LPS		=	cpx
COMP		=	gnu
LINK		=	static
LINKER  	=	CPP

READLINE	=	true
ZLIB		=	true
ZIMPL		=	false

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
LDFLAGS		=	-lpthread -lm
ARFLAGS		=	cr
DFLAGS		=	-MM

GCCWARN		=	-Wall -W -Wpointer-arith -Wcast-align -Wwrite-strings -Wshadow \
			-Wno-unknown-pragmas -Wno-unused-parameter \
			-Wredundant-decls -Wdisabled-optimization \
			-Wsign-compare -Wstrict-prototypes \
			-Wmissing-declarations -Wmissing-prototypes

GXXWARN		=	-Wall -W -Wpointer-arith -Wcast-align -Wwrite-strings -Wshadow \
			-Wno-unknown-pragmas -Wno-unused-parameter \
			-Wredundant-decls -Wdisabled-optimization \
			-Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder \
			-Woverloaded-virtual -Wsign-promo -Wsynth \
			-Wcast-qual -Wno-unused-parameter # -Wold-style-cast -Wshadow -Wundef

BASE		=	$(OSTYPE).$(ARCH).$(COMP).$(OPT).$(LINK)
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
# External Libraries
#-----------------------------------------------------------------------------

ifeq ($(ZLIB_LDFLAGS),)
ZLIB		=	false
endif
ifeq ($(ZLIB),true)
FLAGS		+=	-DWITH_ZLIB $(ZLIB_FLAGS)
LDFLAGS		+=	$(ZLIB_LDFLAGS)
endif

ifeq ($(READLINE_LDFLAGS),)
READLINE	=	false
endif
ifeq ($(READLINE),true)
FLAGS		+=	-DWITH_READLINE $(READLINE_FLAGS)
LDFLAGS		+=	$(READLINE_LDFLAGS)
endif

ifeq ($(ZIMPL),true)
FLAGS		+=	-DWITH_ZIMPL -Ilib/scip/lib/zimplinc
LDFLAGS		+=	-lzimpl.$(OSTYPE).$(ARCH).$(COMP) -lgmp -lz
endif


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
# OBJSCIP Library
#-----------------------------------------------------------------------------

OBJSCIPLIBNAME	=	objscip
OBJSCIPLIB	=	$(OBJSCIPLIBNAME).$(BASE)
OBJSCIPLIBFILE	=	$(LIBDIR)/scip/lib/lib$(OBJSCIPLIB).a


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

.PHONY: all
all:            $(MAINFILE)

.PHONY: lint
lint:		$(MAINSRC)
		$(LINT) lint/$(MAINNAME).lnt -os\(lint.out\) \
		$(FLAGS) -UNDEBUG $^

.PHONY: doc
doc:		
		cd doc; $(DOXY) $(MAINNAME).dxy

.PHONY: test
test:		
		cd check; \
		/bin/sh ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(MAIN).$(MACHINENAME) $(TIME) $(NODES) $(MEM);

$(OBJDIR):	
		-mkdir -p $(OBJDIR)

$(LIBDIR):
		-mkdir -p $(LIBDIR)

$(BINDIR):
		-mkdir -p $(BINDIR)

.PHONY: clean
clean:
		-rm -rf $(OBJDIR)/* $(MAINFILE)

.PHONY: depend
depend:
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-z\_]\{1,\}\)\.o|$$\(OBJDIR\)/\1.o|g'\'' \
		>$(MAINDEP)'

-include	$(MAINDEP)

ifeq ($(LINKER),C)
$(MAINFILE):	$(OBJDIR) $(BINDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(MAINXXX)
		@echo "-> linking $@"
ifeq ($(VERBOSE), true)
		$(CC) $(CFLAGS) $(MAINXXX) \
		-L$(LIBDIR)/scip/lib -l$(SCIPLIB) -l$(LPILIB) -l$(LPSLIB) -L$(LIBDIR) \
		$(LDFLAGS) -o $@
else
		@$(CC) $(CFLAGS) $(MAINXXX) \
		-L$(LIBDIR)/scip/lib -l$(SCIPLIB) -l$(LPILIB) -l$(LPSLIB) -L$(LIBDIR) \
		$(LDFLAGS) -o $@
endif
endif
ifeq ($(LINKER),CPP)
$(MAINFILE):	$(OBJDIR) $(BINDIR) $(SCIPLIBFILE) $(OBJSCIPLIBFILE) $(LPILIBFILE) $(MAINXXX)
		@echo "-> linking $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(CXXFLAGS) $(MAINXXX) \
		-L$(LIBDIR)/scip/lib -l$(SCIPLIB) -l$(OBJSCIPLIB) -l$(LPILIB) -l$(LPSLIB) -L$(LIBDIR) \
		$(LDFLAGS) -o $@
else
		@$(CXX) $(CXXFLAGS) $(MAINXXX) \
		-L$(LIBDIR)/scip/lib -l$(SCIPLIB) -l$(OBJSCIPLIB) -l$(LPILIB) -l$(LPSLIB) -L$(LIBDIR) \
		$(LDFLAGS) -o $@
endif
endif

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
ifeq ($(VERBOSE), true)
		$(CC) $(FLAGS) $(CFLAGS) -c $< -o $@
else
		@$(CC) $(FLAGS) $(CFLAGS) -c $< -o $@
endif

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(FLAGS) $(CXXFLAGS) -c $< -o $@
else
		@$(CXX) $(FLAGS) $(CXXFLAGS) -c $< -o $@
endif

# --- EOF ---------------------------------------------------------------------
