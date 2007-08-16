#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2006 Tobias Achterberg                              *
#*                                                                           *
#*                  2002-2006 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic Licence.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: Makefile,v 1.13 2007/08/16 14:08:55 bzfpfend Exp $

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
FEASTOL		=	default
TEST		=	miplib
SETTINGS        =       default

VERBOSE		=	false
OPT		=	opt
LPS		=	cpx
COMP		=	gnu
LIBEXT		=	a

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
CPLEX		=	cplex

FLAGS		=	-I$(SRCDIR) -I$(LIBDIR)/scip/src
OFLAGS		=
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

BASE		=	$(OSTYPE).$(ARCH).$(COMP).$(OPT)
OBJDIR		=	obj/O.$(BASE)
SRCDIR		=	src
BINDIR		=	bin
LIBDIR		=	lib
EXEEXTENSION	=

#-----------------------------------------------------------------------------
include make/make.$(BASE)
#-----------------------------------------------------------------------------

FLAGS		+=	$(USRFLAGS)
OFLAGS		+=	$(USROFLAGS)
CFLAGS		+=	$(USRCFLAGS)
CXXFLAGS	+=	$(USRCXXFLAGS)
LDFLAGS		+=	$(USRLDFLAGS)
ARFLAGS		+=	$(USRARFLAGS)
DFLAGS		+=	$(USRDFLAGS)


#-----------------------------------------------------------------------------
# Memory Management
#-----------------------------------------------------------------------------

#FLAGS		+=	-DNOSAFEMEM
#FLAGS		+=	-DNOBLOCKMEM


#-----------------------------------------------------------------------------
# LP Solver Interface
#-----------------------------------------------------------------------------

LPILIBNAME	=	lpi$(LPS)

ifeq ($(LPS),cpx)
LPSLDFLAGS	=	-lcplex.$(OSTYPE).$(ARCH).$(COMP)
endif

ifeq ($(LPS),spx)
LPSLDFLAGS	=	-lsoplex.$(OSTYPE).$(ARCH).$(COMP)
endif

ifeq ($(LPS),spxdbg)
LPSLDFLAGS	=	-lsoplexdbg.$(OSTYPE).$(ARCH).$(COMP)
endif

ifeq ($(LPS),spx121)
LPSLDFLAGS	=	-lsoplex121.$(OSTYPE).$(ARCH).$(COMP)
endif

ifeq ($(LPS),clp)
LPSLDFLAGS	=	-L$(LIBDIR)/clpinc/../lib -Wl,-rpath,$(LIBDIR)/clpinc/../lib -lClp -lCoin
endif

ifeq ($(LPS),clpdbg)
LPSLDFLAGS	=	-lclpdbg.$(OSTYPE).$(ARCH).$(COMP) -lcoindbg.$(OSTYPE).$(ARCH).$(COMP)
endif

LPILIB		=	$(LPILIBNAME).$(BASE)
LPILIBFILE	=	$(LIBDIR)/scip/lib/lib$(LPILIB).$(LIBEXT)


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
FLAGS		+=	-DWITH_ZIMPL -Ilib/zimplinc
LDFLAGS		+=	-lzimpl.$(OSTYPE).$(ARCH).$(COMP) -lgmp -lz
endif


#-----------------------------------------------------------------------------
# SCIP Library
#-----------------------------------------------------------------------------

SCIPLIBNAME	=	scip
SCIPLIB		=	$(SCIPLIBNAME).$(BASE)
SCIPLIBFILE	=	$(LIBDIR)/scip/lib/lib$(SCIPLIB).$(LIBEXT)


#-----------------------------------------------------------------------------
# OBJSCIP Library
#-----------------------------------------------------------------------------

OBJSCIPLIBNAME	=	objscip
OBJSCIPLIB	=	$(OBJSCIPLIBNAME).$(BASE)
OBJSCIPLIBFILE	=	$(LIBDIR)/scip/lib/lib$(OBJSCIPLIB).$(LIBEXT)


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	scipmip
MAINOBJ		=	cppmain.o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINDEP		=	$(SRCDIR)/depend.cppmain.$(OPT)

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))


#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

.PHONY: all
all:            $(LIBDIR) $(LIBDIR)/scip $(MAINFILE)

.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) lint/$(MAINNAME).lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

.PHONY: doc
doc:		
		cd doc; $(DOXY) $(MAINNAME).dxy

$(OBJDIR):	
		@-mkdir -p $(OBJDIR)

$(LIBDIR):
		@-mkdir -p $(LIBDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:
		-rm -rf $(OBJDIR)/* $(MAINFILE)

.PHONY: depend
depend:
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-z_/]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(MAINDEP)'

-include	$(MAINDEP)

$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(MAINOBJFILES) \
		-L$(LIBDIR)/scip/lib -l$(SCIPLIB) -l$(OBJSCIPLIB) -l$(LPILIB) $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) -o $@
else
		@$(CXX) $(MAINOBJFILES) \
		-L$(LIBDIR)/scip/lib -l$(SCIPLIB) -l$(OBJSCIPLIB) -l$(LPILIB) $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) -o $@
endif

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
ifeq ($(VERBOSE), true)
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< -o $@
else
		@$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< -o $@
endif

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
ifeq ($(VERBOSE), true)
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< -o $@
else
		@$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< -o $@
endif

$(LIBDIR)/scip:
		@$(SHELL) -ec 'if [ ! -e $@ ] ; \
                        then \
                                echo "** missing soft-link \"$@\"" ; \
                                echo -n "** enter full path of SCIP installation: " ; \
                                read TARGET ; \
                                if [ "$$TARGET" != "" ] ; \
                                then \
                                        echo "-> creating softlink \"$@\" -> \"$$TARGET\"" ; \
                                        ln -s $$TARGET $@ ; \
                                else \
                                        echo "-> skipped creation of softlink \"$@\"." ; \
                                fi ; \
                                echo ; \
                        fi'


# --- EOF ---------------------------------------------------------------------
