#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2007 Tobias Achterberg                              *
#*                                                                           *
#*                  2002-2007 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic Licence.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: Makefile,v 1.230.2.1 2009/06/10 17:47:12 bzfwolte Exp $

#@file    Makefile
#@brief   SCIP Makefile
#@author  Thorsten Koch
#@author  Tobias Achterberg


#-----------------------------------------------------------------------------
# detect host architecture
#-----------------------------------------------------------------------------
include make/make.detecthost



#-----------------------------------------------------------------------------
# default settings
#-----------------------------------------------------------------------------

VERSION		:=	1.00.3

TIME     	=  	3600
NODES           =       2100000000
MEM		=	1536
DISPFREQ	=	10000
FEASTOL		=	default
TEST		=	shortmiplib
SETTINGS        =       default
CONTINUE	=	false
LOCK		=	false
EXACTSOLVE	=	false
EXACTGETTESTSET	=	false

VERBOSE		=	false
OPT		=	opt
COMP		=	gnu
LPS		=	spx
LIBEXT		=	a
SHAREDLIBEXT	=	so
LINKER  	=	C
SOFTLINKS	=

MAKESOFTLINKS	=	true
READLINE	=	true
ZLIB		=	true
ZIMPL		=	true

CC		=	gcc
CC_c		=	-c # the trailing space is important
CC_o		=	-o # the trailing space is important
CXX		=	g++
CXX_c		=	-c # the trailing space is important
CXX_o		=	-o # the trailing space is important
LINKCC		=	gcc
LINKCC_L	=	-L
LINKCC_l	=	-l
LINKCC_o	=	-o # the trailing space is important
LINKCXX		=	g++
LINKCXX_L	=	-L
LINKCXX_l	=	-l
LINKCXX_o	=	-o # the trailing space is important
LINKLIBSUFFIX	=
DCC		=	gcc
DCXX		=	g++
AR		=	ar
AR_o		=
RANLIB		=	ranlib
LINT		=	flexelint
DOXY		=	doxygen
CPLEX		=	cplex
CBC		=	cbc

FLAGS		=	-I$(SRCDIR) -DWITH_SCIPDEF
OFLAGS		=
CFLAGS		=	
CXXFLAGS	=	
LDFLAGS		=	$(LINKCC_l)m$(LINKLIBSUFFIX)
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
BINOBJDIR	=	$(OBJDIR)/bin
LIBOBJDIR	=	$(OBJDIR)/lib
LIBOBJSUBDIRS	=	scip objscip blockmemshell tclique
SRCDIR		=	src
BINDIR		=	bin
LIBDIR		=	lib
EXEEXTENSION	=
ALLSRC		=

LINKSMARKERFILE	=	$(LIBDIR)/linkscreated.$(LPS).$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX).$(ZIMPL)
LASTSETTINGS	=	$(OBJDIR)/make.lastsettings

#-----------------------------------------------------------------------------
include make/make.$(BASE)
-include make/local/make.$(HOSTNAME)
-include make/local/make.$(HOSTNAME).$(COMP)
-include make/local/make.$(HOSTNAME).$(COMP).$(OPT)
#-----------------------------------------------------------------------------

FLAGS		+=	$(USRFLAGS)
OFLAGS		+=	$(USROFLAGS)
CFLAGS		+=	$(USRCFLAGS)
CXXFLAGS	+=	$(USRCXXFLAGS)
LDFLAGS		+=	$(USRLDFLAGS)
ARFLAGS		+=	$(USRARFLAGS)
DFLAGS		+=	$(USRDFLAGS)


#-----------------------------------------------------------------------------
# Exact MIP Solving
#-----------------------------------------------------------------------------
ifeq ($(EXACTSOLVE),true)
FLAGS		+=	-DROUNDING_FE
FLAGS		+=	-DEXACTSOLVE_OUT		# general info 
LDFLAGS		+=	$(LINKCC_l)gmp$(LINKLIBSUFFIX)
endif

ifeq ($(EXACTGETTESTSET),true)
FLAGS		+=	-DUNBNDVARSINFO			# info with nr of unbounded vars in presolved problem; for test set
endif


#-----------------------------------------------------------------------------
# Memory Management
#-----------------------------------------------------------------------------

#FLAGS		+=	-DNOSAFEMEM
#FLAGS		+=	-DNOBLOCKMEM


#-----------------------------------------------------------------------------
# LP Solver Interface
#-----------------------------------------------------------------------------

LPILIBSHORTNAME	=	lpi$(LPS)
LPILIBNAME	=	$(LPILIBSHORTNAME)-$(VERSION)
LPILIBOBJ	=
LPSOPTIONS	=

LPSOPTIONS	+=	cpx
ifeq ($(LPS),cpx)
FLAGS		+=	-I$(LIBDIR)/cpxinc
LPSLDFLAGS	=	$(LINKCC_l)cplex.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCC_l)pthread$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_cpx.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/cpxinc
SOFTLINKS	+=	$(LIBDIR)/libcplex.$(OSTYPE).$(ARCH).$(COMP).$(LIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libcplex.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
endif

LPSOPTIONS	+=	cpx903
ifeq ($(LPS),cpx903)
FLAGS		+=	-I$(LIBDIR)/cpx903inc
LPSLDFLAGS	=	$(LINKCC_l)cplex903.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCC_l)pthread$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_cpx903.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/cpx903inc
SOFTLINKS	+=	$(LIBDIR)/libcplex903.$(OSTYPE).$(ARCH).$(COMP).$(LIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libcplex903.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
endif

LPSOPTIONS	+=	xprs
ifeq ($(LPS),xprs)
FLAGS		+=	-I$(LIBDIR)/xprsinc
LPSLDFLAGS	=	$(LINKCC_l)xpress.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_xprs.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/xprsinc
SOFTLINKS	+=	$(LIBDIR)/libxpress.$(OSTYPE).$(ARCH).$(COMP).$(LIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libxpress.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
endif

LPSOPTIONS	+=	msk
ifeq ($(LPS),msk)
FLAGS		+=	-I$(LIBDIR)/mskinc
LPSLDFLAGS	=	$(LINKCC_l)mosek.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_msk.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/mskinc
SOFTLINKS	+=	$(LIBDIR)/libmosek.$(OSTYPE).$(ARCH).$(COMP).$(LIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libmosek.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
endif

LPSOPTIONS	+=	spx
ifeq ($(LPS),spx)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/spxinc 
LPSLDFLAGS	=	$(LINKCXX_l)soplex.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_spx.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/scip/lpi_spx.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/spxinc
SOFTLINKS	+=	$(LIBDIR)/libsoplex.$(OSTYPE).$(ARCH).$(COMP).$(LIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libsoplex.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
endif

LPSOPTIONS	+=	spxdbg
ifeq ($(LPS),spxdbg)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/spxinc 
LPSLDFLAGS	=	$(LINKCXX_l)soplexdbg.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_spx.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/scip/lpi_spx.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/spxinc
SOFTLINKS	+=	$(LIBDIR)/libsoplexdbg.$(OSTYPE).$(ARCH).$(COMP).$(LIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libsoplexdbg.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
endif

LPSOPTIONS	+=	spx121
ifeq ($(LPS),spx121)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/spx121inc 
LPSLDFLAGS	=	$(LINKCXX_l)soplex121.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_spx121.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/scip/lpi_spx121.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/spx121inc
SOFTLINKS	+=	$(LIBDIR)/libsoplex121.$(OSTYPE).$(ARCH).$(COMP).$(LIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libsoplex121.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
endif

LPSOPTIONS	+=	clp
ifeq ($(LPS),clp)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/clpinc
LPSLDFLAGS	=	$(LINKCXX_L)$(LIBDIR) -Wl,-rpath,$(LIBDIR) $(LINKCXX_l)clp.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCXX_l)coinutils.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) \
			$(LINKCXX_l)gfortran$(LINKLIBSUFFIX) $(LINKCXX_l)lapack$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_clp.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/scip/lpi_clp.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/clpinc
SOFTLINKS	+=	$(LIBDIR)/libclp.$(OSTYPE).$(ARCH).$(COMP).$(LIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libclp.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libcoinutils.$(OSTYPE).$(ARCH).$(COMP).$(LIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libcoinutils.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
endif

LPSOPTIONS	+=	clpdbg
ifeq ($(LPS),clpdbg)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/clpinc
LPSLDFLAGS	=	$(LINKCXX_L)$(LIBDIR) -Wl,-rpath,$(LIBDIR) $(LINKCXX_l)clpdbg.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCXX_l)coinutilsdbg.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) \
			$(LINKCXX_l)gfortran$(LINKLIBSUFFIX) $(LINKCXX_l)lapack$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_clp.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/scip/lpi_clp.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/clpinc
SOFTLINKS	+=	$(LIBDIR)/libclpdbg.$(OSTYPE).$(ARCH).$(COMP).$(LIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libclpdbg.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libcoinutilsdbg.$(OSTYPE).$(ARCH).$(COMP).$(LIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libcoinutilsdbg.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
endif

LPSOPTIONS	+=	none
ifeq ($(LPS),none)
LPILIBOBJ	=	scip/lpi_none.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
endif

LPILIB		=	$(LPILIBNAME).$(BASE)
LPILIBFILE	=	$(LIBDIR)/lib$(LPILIB).$(LIBEXT)
LPILIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(LPILIBOBJ))
LPILIBDEP	=	$(SRCDIR)/depend.lpilib.$(LPS).$(OPT)
LPILIBLINK	=	$(LIBDIR)/lib$(LPILIBSHORTNAME).$(BASE).$(LIBEXT)
ALLSRC		+=	$(LPILIBSRC)


#-----------------------------------------------------------------------------
# External Libraries
#-----------------------------------------------------------------------------

ZLIBDEP		:=	$(SRCDIR)/depend.zlib
ZLIBSRC		:=	$(shell cat $(ZLIBDEP))
ifeq ($(ZLIB_LDFLAGS),)
ZLIB		=	false
endif
ifeq ($(ZLIB),true)
FLAGS		+=	-DWITH_ZLIB $(ZLIB_FLAGS)
LDFLAGS		+=	$(ZLIB_LDFLAGS)
endif

READLINEDEP	:=	$(SRCDIR)/depend.readline
READLINESRC	:=	$(shell cat $(READLINEDEP))
ifeq ($(READLINE_LDFLAGS),)
READLINE	=	false
endif
ifeq ($(READLINE),true)
FLAGS		+=	-DWITH_READLINE $(READLINE_FLAGS)
LDFLAGS		+=	$(READLINE_LDFLAGS)
endif

ZIMPLDEP	:=	$(SRCDIR)/depend.zimpl
ZIMPLSRC	:=	$(shell cat $(ZIMPLDEP))
ifeq ($(ZIMPL),true)
ifeq ($(ZLIB),false)
$(error ZIMPL requires the ZLIB to be linked. Use either ZIMPL=false or ZLIB=true)
endif
FLAGS		+=	-DWITH_ZIMPL -I$(LIBDIR)/zimplinc $(ZIMPL_FLAGS)
LDFLAGS		+=	$(LINKCC_l)zimpl.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCC_l)gmp$(LINKLIBSUFFIX) $(ZIMPL_LDFLAGS)
DIRECTORIES	+=	$(LIBDIR)/zimplinc
SOFTLINKS	+=	$(LIBDIR)/zimplinc/zimpl
SOFTLINKS	+=	$(LIBDIR)/libzimpl.$(OSTYPE).$(ARCH).$(COMP).$(LIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libzimpl.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
endif


#-----------------------------------------------------------------------------
# SCIP Library
#-----------------------------------------------------------------------------

SCIPLIBSHORTNAME=	scip
SCIPLIBNAME	=	$(SCIPLIBSHORTNAME)-$(VERSION)
SCIPLIBOBJ	=	scip/branch.o \
			scip/buffer.o \
			scip/clock.o \
			scip/conflict.o \
			scip/cons.o \
			scip/cutpool.o \
			scip/debug.o \
			scip/dialog.o \
			scip/disp.o \
			scip/event.o \
			scip/fileio.o \
			scip/heur.o \
			scip/history.o \
			scip/implics.o \
			scip/interrupt.o \
			scip/intervalarith.o \
			scip/lp.o \
			scip/mem.o \
			scip/misc.o \
			scip/nodesel.o \
			scip/paramset.o \
			scip/presol.o \
			scip/pricestore.o \
			scip/pricer.o \
			scip/primal.o \
			scip/prob.o \
			scip/prop.o \
			scip/reader.o \
			scip/relax.o \
			scip/retcode.o \
			scip/scip.o \
			scip/scipdefplugins.o \
			scip/scipshell.o \
			scip/sepa.o \
			scip/sepastore.o \
			scip/set.o \
			scip/sol.o \
			scip/solve.o \
			scip/stat.o \
			scip/tree.o \
			scip/var.o \
			scip/vbc.o \
			scip/branch_allfullstrong.o \
			scip/branch_fullstrong.o \
			scip/branch_inference.o \
			scip/branch_mostinf.o \
			scip/branch_leastinf.o \
			scip/branch_pscost.o \
			scip/branch_random.o \
			scip/branch_relpscost.o \
			scip/cons_and.o \
			scip/cons_binpack.o \
			scip/cons_bounddisjunction.o \
			scip/cons_conjunction.o \
			scip/cons_eqknapsack.o \
			scip/cons_integral.o \
			scip/cons_invarknapsack.o \
			scip/cons_knapsack.o \
			scip/cons_linear.o \
			scip/cons_logicor.o \
			scip/cons_or.o \
			scip/cons_setppc.o \
			scip/cons_sos1.o \
			scip/cons_sos2.o \
			scip/cons_varbound.o \
			scip/cons_xor.o \
			scip/dialog_default.o \
			scip/disp_default.o \
			scip/heur_actconsdiving.o \
			scip/heur_coefdiving.o \
			scip/heur_crossover.o \
			scip/heur_feaspump.o \
			scip/heur_fixandinfer.o \
			scip/heur_fracdiving.o \
			scip/heur_guideddiving.o \
			scip/heur_intdiving.o \
			scip/heur_intshifting.o \
			scip/heur_linesearchdiving.o \
			scip/heur_localbranching.o \
			scip/heur_mutation.o \
			scip/heur_objpscostdiving.o \
			scip/heur_octane.o \
			scip/heur_oneopt.o \
			scip/heur_pscostdiving.o \
			scip/heur_rens.o \
			scip/heur_rins.o \
			scip/heur_rootsoldiving.o \
			scip/heur_rounding.o \
			scip/heur_shifting.o \
			scip/heur_simplerounding.o \
			scip/heur_veclendiving.o \
			scip/nodesel_bfs.o \
			scip/nodesel_dfs.o \
			scip/nodesel_estimate.o \
			scip/nodesel_hybridestim.o \
			scip/nodesel_restartdfs.o \
			scip/presol_dualfix.o \
			scip/presol_implics.o \
			scip/presol_inttobinary.o \
			scip/presol_probing.o \
			scip/presol_trivial.o \
			scip/prop_pseudoobj.o \
			scip/prop_rootredcost.o \
			scip/reader_cip.o \
			scip/reader_cnf.o \
			scip/reader_fix.o \
			scip/reader_lp.o \
			scip/reader_mps.o \
			scip/reader_rlp.o \
			scip/reader_sol.o \
			scip/reader_zpl.o \
			scip/sepa_clique.o \
			scip/sepa_cmir.o \
			scip/sepa_flowcover.o \
			scip/sepa_gomory.o \
			scip/sepa_impliedbounds.o \
			scip/sepa_intobj.o \
			scip/sepa_mcf.o \
			scip/sepa_redcost.o \
			scip/sepa_strongcg.o \
			tclique/tclique_branch.o \
			tclique/tclique_coloring.o \
			tclique/tclique_graph.o

SCIPLIB		=	$(SCIPLIBNAME).$(BASE)
SCIPLIBFILE	=	$(LIBDIR)/lib$(SCIPLIB).$(LIBEXT)
SCIPLIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(SCIPLIBOBJ))
SCIPLIBSRC	=	$(addprefix $(SRCDIR)/,$(SCIPLIBOBJ:.o=.c))
SCIPLIBDEP	=	$(SRCDIR)/depend.sciplib.$(OPT)
SCIPLIBLINK	=	$(LIBDIR)/lib$(SCIPLIBSHORTNAME).$(BASE).$(LIBEXT)
ALLSRC		+=	$(SCIPLIBSRC)


#-----------------------------------------------------------------------------
# Objective SCIP Library
#-----------------------------------------------------------------------------

OBJSCIPLIBSHORTNAME=	objscip
OBJSCIPLIBNAME	=	$(OBJSCIPLIBSHORTNAME)-$(VERSION)
OBJSCIPLIBOBJ	=	objscip/objbranchrule.o \
			objscip/objconshdlr.o \
			objscip/objdialog.o \
			objscip/objdisp.o \
			objscip/objeventhdlr.o \
			objscip/objheur.o \
			objscip/objmessagehdlr.o \
			objscip/objnodesel.o \
			objscip/objpresol.o \
			objscip/objpricer.o \
			objscip/objprobdata.o \
			objscip/objprop.o \
			objscip/objreader.o \
			objscip/objrelax.o \
			objscip/objsepa.o \
			objscip/objvardata.o

OBJSCIPLIB	=	$(OBJSCIPLIBNAME).$(BASE)
OBJSCIPLIBFILE	=	$(LIBDIR)/lib$(OBJSCIPLIB).$(LIBEXT)
OBJSCIPLIBOBJFILES=	$(addprefix $(LIBOBJDIR)/,$(OBJSCIPLIBOBJ))
OBJSCIPLIBSRC	=	$(addprefix $(SRCDIR)/,$(OBJSCIPLIBOBJ:.o=.cpp))
OBJSCIPLIBDEP	=	$(SRCDIR)/depend.objsciplib.$(OPT)
OBJSCIPLIBLINK	=	$(LIBDIR)/lib$(OBJSCIPLIBSHORTNAME).$(BASE).$(LIBEXT)
ALLSRC		+=	$(OBJSCIPLIBSRC)


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINSHORTNAME	=	scip
MAINNAME	=	$(MAINSHORTNAME)-$(VERSION)

ifeq ($(LINKER),C)
MAINOBJ		=	cmain.o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.c))
MAINDEP		=	$(SRCDIR)/depend.cmain.$(OPT)
endif
ifeq ($(LINKER),CPP)
MAINOBJ		=	cppmain.o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINDEP		=	$(SRCDIR)/depend.cppmain.$(OPT)
endif

MAINFILE	=	$(BINDIR)/$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINOBJFILES	=	$(addprefix $(BINOBJDIR)/,$(MAINOBJ))
MAINLINK	=	$(BINDIR)/$(MAINSHORTNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINSHORTLINK	=	$(BINDIR)/$(MAINSHORTNAME)$(EXEEXTENSION)
ALLSRC		+=	$(MAINSRC)



#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(SCIPLIBFILE) $(OBJSCIPLIBFILE) $(LPILIBFILE) $(LPILIBLINK) $(SCIPLIBLINK) $(OBJSCIPLIBLINK) $(MAINLINK) $(MAINSHORTLINK) \
		$(LPILIBOBJFILES) $(SCIPLIBOBJFILES) $(OBJSCIPLIBOBJFILES) $(MAINOBJFILES)
endif

.PHONY: all
all:            checklpsdefine $(LINKSMARKERFILE) $(SCIPLIBFILE) $(OBJSCIPLIBFILE) $(LPILIBFILE) $(MAINFILE) $(LPILIBLINK) $(SCIPLIBLINK) $(OBJSCIPLIBLINK) $(MAINLINK) $(MAINSHORTLINK)

.PHONY: lint
lint:		$(SCIPLIBSRC) $(OBJSCIPLIBSRC) $(LPILIBSRC) $(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) lint/$(MAINSHORTNAME).lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

.PHONY: doc
doc:		
		cd doc; $(DOXY) $(MAINSHORTNAME).dxy

.PHONY: check
check:		test

.PHONY: test
test:		
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION);

.PHONY: testpre
testpre:		
		cd check; \
		$(SHELL) ./checkpre.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(FEASTOL) $(DISPFREQ) $(CONTINUE);

.PHONY: testcplex
testcplex:		
		cd check; \
		$(SHELL) ./check_cplex.sh $(TEST) $(CPLEX) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(FEASTOL) 0.0 $(CONTINUE);

.PHONY: testcbc
testcbc:		
		cd check; \
		$(SHELL) ./check_cbc.sh $(TEST) $(CBC) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(FEASTOL) 0.0 $(CONTINUE);

$(LPILIBLINK):	$(LPILIBFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(LPILIBFILE)) $(notdir $@)

$(SCIPLIBLINK):	$(SCIPLIBFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(SCIPLIBFILE)) $(notdir $@)

$(OBJSCIPLIBLINK):	$(OBJSCIPLIBFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(OBJSCIPLIBFILE)) $(notdir $@)

$(MAINLINK) $(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):	
		@-mkdir -p $(OBJDIR)

$(BINOBJDIR):	$(OBJDIR)
		@-mkdir -p $(BINOBJDIR)

$(LIBOBJDIR):	$(OBJDIR)
		@-mkdir -p $(LIBOBJDIR)

$(LIBOBJSUBDIRS):	$(LIBOBJDIR)
		@-mkdir -p $(LIBOBJDIR)/$@

$(LIBDIR):
		@-mkdir -p $(LIBDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:
		-rm -rf $(OBJDIR)/* $(SCIPLIBFILE) $(OBJSCIPLIBFILE) $(LPILIBFILE) $(MAINFILE) \
		$(LPILIBLINK) $(SCIPLIBLINK) $(OBJSCIPLIBLINK) $(MAINLINK) $(MAINSHORTLINK)

.PHONY: lpidepend
lpidepend:
ifeq ($(LINKER),C)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(LPILIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(LPILIBDEP)'
endif
ifeq ($(LINKER),CPP)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(LPILIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(LPILIBDEP)'
endif

.PHONY: maindepend
maindepend:
ifeq ($(LINKER),C)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(BINOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(MAINDEP)'
endif
ifeq ($(LINKER),CPP)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(BINOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(MAINDEP)'
endif

.PHONY: depend
depend:		lpidepend maindepend
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(SCIPLIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(SCIPLIBDEP)'
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(OBJSCIPLIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(OBJSCIPLIBDEP)'
		@echo `grep -l "WITH_ZLIB" $(ALLSRC)` >$(ZLIBDEP)
		@echo `grep -l "WITH_READLINE" $(ALLSRC)` >$(READLINEDEP)
		@echo `grep -l "WITH_ZIMPL" $(ALLSRC)` >$(ZIMPLDEP)

-include	$(MAINDEP)
-include	$(SCIPLIBDEP)
-include	$(OBJSCIPLIBDEP)
-include 	$(LPILIBDEP)

$(MAINFILE):	$(BINDIR) $(BINOBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
ifeq ($(LINKER),C)
		$(LINKCC) $(MAINOBJFILES) \
		$(LINKCC_L)$(LIBDIR) $(LINKCC_l)$(SCIPLIB)$(LINKLIBSUFFIX) $(LINKCC_l)$(LPILIB)$(LINKLIBSUFFIX) \
		$(OFLAGS) $(LPSLDFLAGS) $(LDFLAGS) $(LINKCC_o)$@
endif
ifeq ($(LINKER),CPP)
		$(LINKCXX) $(MAINOBJFILES) \
		$(LINKCXX_L)$(LIBDIR) $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) \
		$(OFLAGS) $(LPSLDFLAGS) $(LDFLAGS) $(LINKCXX_o)$@
endif

$(SCIPLIBFILE):	checklpsdefine $(LIBOBJSUBDIRS) $(LIBDIR) touchexternal $(SCIPLIBOBJFILES) 
		@echo "-> generating library $@"
		-rm -f $@
		$(AR) $(ARFLAGS) $(AR_o)$@ $(SCIPLIBOBJFILES) 
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

$(OBJSCIPLIBFILE):	$(LIBOBJSUBDIRS) $(LIBDIR) $(OBJSCIPLIBOBJFILES) 
		@echo "-> generating library $@"
		-rm -f $@
		$(AR) $(ARFLAGS) $(AR_o)$@ $(OBJSCIPLIBOBJFILES) 
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

$(LPILIBFILE):	$(LIBOBJSUBDIRS) $(LIBDIR) $(LPILIBOBJFILES)
		@echo "-> generating library $@"
		-rm -f $@
		$(AR) $(ARFLAGS) $(AR_o)$@ $(LPILIBOBJFILES)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

$(BINOBJDIR)/%.o:	$(SRCDIR)/%.c $(BINOBJDIR)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) $(CC_c)$< $(CC_o)$@

$(BINOBJDIR)/%.o:	$(SRCDIR)/%.cpp $(BINOBJDIR)
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) $(CXX_c)$< $(CXX_o)$@

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.c $(LIBOBJDIR)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CFLAGS) $(CC_c)$< $(CC_o)$@

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.cpp $(LIBOBJDIR)
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CXXFLAGS) $(CXX_c)$< $(CXX_o)$@


-include $(LASTSETTINGS)

.PHONY: touchexternal
touchexternal:	$(ZLIBDEP) $(READLINEDEP) $(ZIMPLDEP)
ifneq ($(ZLIB),$(LAST_ZLIB))
		@-touch $(ZLIBSRC)
endif
ifneq ($(READLINE),$(LAST_READLINE))
		@-touch $(READLINESRC)
endif
ifneq ($(ZIMPL),$(LAST_ZIMPL))
		@-touch $(ZIMPLSRC)
endif
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_ZLIB=$(ZLIB)" >> $(LASTSETTINGS)
		@echo "LAST_READLINE=$(READLINE)" >> $(LASTSETTINGS)
		@echo "LAST_ZIMPL=$(ZIMPL)" >> $(LASTSETTINGS)

$(LINKSMARKERFILE):
		@$(MAKE) links

.PHONY: links
links:		echosoftlinks $(LIBDIR) $(DIRECTORIES) $(SOFTLINKS)
		@rm -f $(LINKSMARKERFILE)
		@echo "this is only a marker" > $(LINKSMARKERFILE)

.PHONY: echosoftlinks
echosoftlinks:
		@echo
		@echo "** creating softlinks: LPS=$(LPS) OSTYPE=$(OSTYPE) ARCH=$(ARCH) COMP=$(COMP) SUFFIX=$(LINKLIBSUFFIX) ZIMPL=$(ZIMPL)"
		@echo

$(DIRECTORIES):
		@echo "** creating directory \"$@\""
		@-mkdir -p $@

.PHONY: $(SOFTLINKS)
$(SOFTLINKS):
ifeq ($(MAKESOFTLINKS), true)
		@$(SHELL) -ec 'if test ! -e $@ ; \
			then \
				DIRNAME=`dirname $@` ; \
				BASENAMEA=`basename $@ .$(LIBEXT)` ; \
				BASENAMESO=`basename $@ .$(SHAREDLIBEXT)` ; \
				echo "** missing soft-link \"$@\"" ; \
				if [ -e $$DIRNAME/$$BASENAMEA.$(SHAREDLIBEXT) ] ; \
				then \
					echo "** this soft-link is not necessarily needed since \"$$DIRNAME/$$BASENAMEA.$(SHAREDLIBEXT)\" already exists - press return to skip" ; \
				fi ; \
				if [ -e $$DIRNAME/$$BASENAMESO.$(LIBEXT) ] ; \
				then \
					echo "** this soft-link is not necessarily needed since \"$$DIRNAME/$$BASENAMESO.$(LIBEXT)\" already exists - press return to skip" ; \
				fi ; \
				echo -n "** enter soft-link target file or directory for \"$@\" (return if not needed): " ; \
				read TARGET ; \
				if [ "$$TARGET" != "" ] ; \
				then \
					echo "-> creating softlink \"$@\" -> \"$$TARGET\"" ; \
					rm -f $@ ; \
					ln -s $$TARGET $@ ; \
				else \
					echo "-> skipped creation of softlink \"$@\". Call \"make links\" if needed later." ; \
				fi ; \
				echo ; \
			fi'
endif

.PHONY: checklpsdefine
checklpsdefine:
ifeq ($(LPILIBOBJ),)
		$(error invalid LP solver selected: LPS=$(LPS). Possible options are: $(LPSOPTIONS))
endif


# --- EOF ---------------------------------------------------------------------
# DO NOT DELETE
