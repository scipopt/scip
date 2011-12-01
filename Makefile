#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic Licence.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#@file    Makefile
#@brief   SCIP Makefile
#@author  Thorsten Koch
#@author  Tobias Achterberg
#@author  Marc Pfetsch
#@author  Kati Wolter

#-----------------------------------------------------------------------------
# detect host architecture
#-----------------------------------------------------------------------------
include make/make.detecthost



#-----------------------------------------------------------------------------
# default settings
#-----------------------------------------------------------------------------

VERSION		:=	1.2.0.8

TIME     	=  	3600
NODES           =       2100000000
MEM		=	1536
THREADS         =       1
DISPFREQ	=	10000
FEASTOL		=	default
TEST		=	shortmiplib
SETTINGS        =       default
CONTINUE	=	false
LOCK		=	false
REDUCEDSOLVE	=	true
EXACTSOLVE	=	true
BRANCHPLGS	=	true

VERBOSE		=	false
OPT		=	opt
COMP		=	gnu
LPS		=	spx
LPSEX		=	none
STATICLIBEXT	=	a
SHAREDLIBEXT	=	so
LIBEXT		=	$(STATICLIBEXT)
LINKER  	=	C
SOFTLINKS	=

MAKESOFTLINKS	=	true
READLINE	=	true
ZLIB		=	true
GMP             =       auto
ZIMPL		=	true
IPOPT		=	false
LPSOPT		=	opt
LPSEXOPT	=	opt
ZIMPLOPT	=	opt
IPOPTOPT	=	opt

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
LIBBUILD	=	$(AR)
LIBBUILD_o	=	$(AR_o)
LIBBUILDFLAGS	=       $(ARFLAGS)
LINT		=	flexelint
DOXY		=	doxygen
CPLEX		=	cplex
CBC		=	cbc
CBCPARALLEL	=	cbc-parallel
MOSEK           =       mosek
GUROBI          =       gurobi.sh
GLPK            =       glpsol
SYMPHONY        =       symphony
BLIS            =       blis


SHELL		= 	bash
READ		=	read -e


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
			-Wmissing-declarations -Wmissing-prototypes # -Wdeclaration-after-statement

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

# define to be able to locate library files
SCIPDIR		=	$(realpath .)

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
ifeq ($(REDUCEDSOLVE),true)
ifeq ($(ZIMPL),false)
$(error REDUCEDSOLVE requires ZIMPL. Use either REDUCEDSOLVE=false or ZIMPL=true)
endif
FLAGS		+=	-DREDUCEDSOLVE			# flag for switching between full and reduced version of SCIP
endif

ifeq ($(EXACTSOLVE),true)
ifeq ($(REDUCEDSOLVE),false)
$(error EXACTSOLVE requires REDUCEDSOLVE. Use either EXACTSOLVE=false or REDUCEDSOLVE=true)
endif
ifeq ($(GMP),false)
$(error EXACTSOLVE requires GMP. Use either EXACTSOLVE=false or GMP=true)
endif
FLAGS		+=	-DEXACTSOLVE			# flag for switching between exact and inexact mode
LDFLAGS		+=	$(LINKCC_l)mpfr$(LINKLIBSUFFIX)
LIBOBJSUBDIRS	+=	rectlu 
endif

ifeq ($(BRANCHPLGS),true)
ifeq ($(REDUCEDSOLVE),false)
$(error BRANCHPLGS requires REDUCEDSOLVE. Use either BRANCHPLGS=false or REDUCEDSOLVE=true)
endif
FLAGS		+=	-DBRANCHPLGS			# flag for including branching plugins or supporting only first fractional variable branching
endif



#-----------------------------------------------------------------------------
# Memory Management
#-----------------------------------------------------------------------------

#FLAGS		+=	-DBMS_NOSAFEMEM
#FLAGS		+=	-DBMS_NOBLOCKMEM


#-----------------------------------------------------------------------------
# LP Solver Interface
#-----------------------------------------------------------------------------

LPILIBSHORTNAME	=	lpi$(LPS)
LPILIBNAME	=	$(LPILIBSHORTNAME)-$(VERSION)
LPILIBOBJ	=
LPSOPTIONS	=
LPIINSTMSG	=


LPSOPTIONS	+=	cpx
ifeq ($(LPS),cpx)
FLAGS		+=	-I$(LIBDIR)/cpxinc
LPSLDFLAGS	=	$(LINKCC_l)cplex.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCC_l)pthread$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_cpx.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/cpxinc
SOFTLINKS	+=	$(LIBDIR)/libcplex.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libcplex.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
LPIINSTMSG	=	"  -> \"cpxinc\" is the path to the CPLEX \"include\" directory, e.g., \"<CPLEX-path>/include/ilcplex\".\n"
LPIINSTMSG	+=	" -> \"libcplex.*\" is the path to the CPLEX library, e.g., \"<CPLEX-path>/lib/x86_rhel4.0_3.4/static_pic/libcplex.a\""
endif

LPSOPTIONS	+=	xprs
ifeq ($(LPS),xprs)
FLAGS		+=	-I$(LIBDIR)/xprsinc
LPSLDFLAGS	=	$(LINKCC_l)xpress.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_xprs.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/xprsinc
SOFTLINKS	+=	$(LIBDIR)/libxpress.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libxpress.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
LPIINSTMSG	=	"  -> \"xprsinc\" is the path to the XPRESS \"include\" directory, e.g., \"<XPRESS-path>/include\".\n"
LPIINSTMSG	+=	" -> \"libpress.*\" is the path to the XPRESS library, e.g., \"<XPRESS-path>/lib/libxpress.a\""
endif

LPSOPTIONS	+=	msk
ifeq ($(LPS),msk)
FLAGS		+=	-I$(LIBDIR)/mskinc
LPSLDFLAGS	=	$(LINKCC_l)mosek.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) \
			$(LINKCXX_l)iomp5.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)$(LINKLIBSUFFIX) $(LINKCC_l)pthread$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_msk.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/mskinc
SOFTLINKS	+=	$(LIBDIR)/libmosek.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libmosek.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
LPIINSTMSG	=	"  -> \"mskinc\" is the path to the Mosek \"include\" directory, e.g., \"<Mosek-path>/include\".\n"
LPIINSTMSG	+=	" -> \"libmosek.*\" is the path to the Mosek library, e.g., \"<Mosek-path>/lib/libmosek.a\".\n"
LPIINSTMSG	+=	" -> \"libiomp5.*\" is the path to the libiomp5, e.g., \"<Mosek-path>/lib/libiomp5.a\""
endif

LPSOPTIONS	+=	spx
ifeq ($(LPS),spx)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/spxinc 
LPSLDFLAGS	=	$(LINKCXX_l)soplex.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_spx.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/scip/lpi_spx.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/spxinc
SOFTLINKS	+=	$(LIBDIR)/libsoplex.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libsoplex.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(SHAREDLIBEXT)
LPIINSTMSG	=	"  -> \"spxinc\" is the path to the SoPlex \"src\" directory, e.g., \"../../soplex/src\".\n"
LPIINSTMSG	+=	" -> \"libsoplex.*\" is the path to the SoPlex library, e.g., \"../../soplex/lib/libsoplex.linux.x86.gnu.opt.a\""
endif

LPSOPTIONS	+=	spx132
ifeq ($(LPS),spx132)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/spx132inc 
LPSLDFLAGS	=	$(LINKCXX_l)soplex132.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_spx132.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/scip/lpi_spx132.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/spx132inc
SOFTLINKS	+=	$(LIBDIR)/libsoplex132.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libsoplex132.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(SHAREDLIBEXT)
LPIINSTMSG	=	"  -> \"spxinc\" is the path to the SoPlex 1.3.2 \"src\" directory, e.g., \"../../soplex-132/src\".\n"
LPIINSTMSG	+=	" -> \"libsoplex.*\" is the path to the SoPlex library, e.g., \"../../soplex/lib/libsoplex-1.3.2.linux.x86.gnu.opt.a\""
endif

LPSOPTIONS	+=	clp
ifeq ($(LPS),clp)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/clpinc
LPSLDFLAGS	=	$(LINKCXX_L)$(LIBDIR) $(LINKCXX_l)clp.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)$(LINKLIBSUFFIX) \
			$(LINKCXX_l)coinutils.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)$(LINKLIBSUFFIX) \
			$(LINKCXX_l)bz2$(LINKLIBSUFFIX) $(LINKCXX_l)lapack$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_clp.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/scip/lpi_clp.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/clpinc
SOFTLINKS	+=	$(LIBDIR)/libclp.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libclp.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(SHAREDLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libcoinutils.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libcoinutils.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(SHAREDLIBEXT)
LPIINSTMSG	=	"  -> \"clpinc\" is the path to the Clp \"include\" directory, e.g., \"<Clp-path>/include/coin\".\n"
LPIINSTMSG	+=	" -> \"libclp.*\" is the path to the Clp library, e.g., \"<Clp-path>/lib/libclp.a\"\n"
LPIINSTMSG	+=	" -> \"libcoinutils.*\" is the path to the COIN utilities library, e.g., \"<Clp-path>/lib/libcoinutils.a\""
endif

LPSOPTIONS	+=	qso
ifeq ($(LPS),qso)
FLAGS         	+=      -I$(LIBDIR)/qsinc
LPSLDFLAGS    	=       $(LINKCC_l)qsopt.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCC_l)pthread$(LINKLIBSUFFIX)
LPILIBOBJ     	= 	scip/lpi_qso.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC     	=       $(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS     	+=      $(LIBDIR)/qsinc
SOFTLINKS     	+=      $(LIBDIR)/libqsopt.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
LPIINSTMSG	=	"  -> \"qsinc\" is the path to the QSopt \"include\" directory, e.g., \"<QSopt-path>\".\n"
LPIINSTMSG	+=	" -> \"libqsopt.*\" is the path to the QSopt library, e.g., \"<QSopt-path>/libqsopt.a\""
endif

LPSOPTIONS	+=	grb
ifeq ($(LPS),grb)
FLAGS		+=	-I$(LIBDIR)/grbinc
LPSLDFLAGS	=	$(LINKCC_l)gurobi.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCC_l)pthread$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_grb.o blockmemshell/memory.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/grbinc
SOFTLINKS	+=	$(LIBDIR)/libgurobi.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libgurobi.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
LPIINSTMSG	=	"  -> \"grbinc\" is the path to the Gurobi \"include\" directory, e.g., \"<Gurobi-path>/include\".\n"
LPIINSTMSG	+=	" -> \"libgurobi.*\" is the path to the Gurobi library, e.g., \"<Gurobi-path>/lib/libgurobi.so\""
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
# Exact LP Solver Interface
#-----------------------------------------------------------------------------

LPIEXLIBSHORTNAME =	lpiex$(LPSEX)
LPIEXLIBNAME	=	$(LPIEXLIBSHORTNAME)-$(VERSION)
LPIEXLIBOBJ	=
LPIEXLIBSRC	=
LPIEXUSED	=
LPSEXOPTIONS	=


LPSEXOPTIONS	+=	qsoex
ifeq ($(LPSEX),qsoex)
FLAGS         	+=      -I$(LIBDIR)/qsexinc
FLAGS         	+=      -I$(LIBDIR)/EGlibinc
LPSEXLDFLAGS   	=       $(LINKCC_l)qsoptex.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCC_l)pthread$(LINKLIBSUFFIX) \
			$(LINKCC_l)eglib.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
LPIEXLIBOBJ    	= 	scip/lpiex_qsoex.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPIEXLIBSRC    	=       $(addprefix $(SRCDIR)/,$(LPIEXLIBOBJ:.o=.c))
SOFTLINKS     	+=      $(LIBDIR)/qsexinc
SOFTLINKS     	+=      $(LIBDIR)/libqsoptex.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
SOFTLINKS     	+=      $(LIBDIR)/EGlibinc
SOFTLINKS     	+=      $(LIBDIR)/libeglib.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
LPIINSTMSG	+=	"\n  -> \"qsexinc\" is the path to the QSopt_ex \"include\" directory, e.g., \"<QSopt_ex-path>\".\n"
LPIINSTMSG	+=	" -> \"libqsoptex.*\" is the path to the QSopt_ex library, e.g., \"<QSopt_ex-path>/lib/qsoptex.a\".\n"
LPIINSTMSG	+=	" -> \"EGlibinc\" is the path to the EGlib \"include\" directory, e.g., \"<QSopt_ex-path/EGlib>\".\n"
LPIINSTMSG	+=	" -> \"libeglib.*\" is the path to the QSopt_ex library, e.g., \"<QSopt_ex-path>/EGlib/lib/EGlib.a\".\n"
LPIEXUSED	=	qsoex
endif

LPSEXOPTIONS	+=	none
ifeq ($(LPSEX),none)
LPSEXLDFLAGS   	=
LPIEXUSED	=	none
endif

LPIEXLIB	=	$(LPIEXLIBNAME).$(BASE)
LPIEXLIBFILE	=	$(LIBDIR)/lib$(LPIEXLIB).$(LIBEXT)
LPIEXLIBOBJFILES =	$(addprefix $(LIBOBJDIR)/,$(LPIEXLIBOBJ))
LPIEXLIBDEP	=	$(SRCDIR)/depend.lpiexlib.$(LPSEX).$(OPT)
LPIEXLIBLINK	=	$(LIBDIR)/lib$(LPIEXLIBSHORTNAME).$(BASE).$(LIBEXT)
ALLSRC		+=	$(LPIEXLIBSRC)

#-----------------------------------------------------------------------------
# External Libraries
#-----------------------------------------------------------------------------

ZLIBDEP		:=	$(SRCDIR)/depend.zlib
ZLIBSRC		:=	$(shell cat $(ZLIBDEP))
ifeq ($(ZLIB_LDFLAGS),)
ZLIB		=	false
endif

GMPDEP		:=	$(SRCDIR)/depend.gmp
GMPSRC		:=	$(shell cat $(GMPDEP))
ifeq ($(GMP),auto)
GMP		=	$(ZIMPL)
endif
ifeq ($(GMP_LDFLAGS),)
GMP		=	false
endif

READLINEDEP	:=	$(SRCDIR)/depend.readline
READLINESRC	:=	$(shell cat $(READLINEDEP))
ifeq ($(READLINE_LDFLAGS),)
READLINE	=	false
endif

ZIMPLDEP	:=	$(SRCDIR)/depend.zimpl
ZIMPLSRC	:=	$(shell cat $(ZIMPLDEP))
ifeq ($(ZIMPL),true)
ifeq ($(ZLIB),false)
$(error ZIMPL requires the ZLIB to be linked. Use either ZIMPL=false or ZLIB=true)
endif
ifeq ($(GMP),false)
$(error ZIMPL requires the GMP to be linked. Use either ZIMPL=false or GMP=auto)
endif
FLAGS		+=	-DWITH_ZIMPL -I$(LIBDIR)/zimplinc $(ZIMPL_FLAGS)
LDFLAGS		+=	$(LINKCC_l)zimpl.$(OSTYPE).$(ARCH).$(COMP).$(ZIMPLOPT)$(LINKLIBSUFFIX) $(ZIMPL_LDFLAGS)
DIRECTORIES	+=	$(LIBDIR)/zimplinc
SOFTLINKS	+=	$(LIBDIR)/zimplinc/zimpl
SOFTLINKS	+=	$(LIBDIR)/libzimpl.$(OSTYPE).$(ARCH).$(COMP).$(ZIMPLOPT).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libzimpl.$(OSTYPE).$(ARCH).$(COMP).$(ZIMPLOPT).$(SHAREDLIBEXT)
LPIINSTMSG	+=	"\n  -> \"zimplinc\" is a directory containing the path to the ZIMPL \"src\" directory, e.g., \"../../zimpl/src\".\n"
LPIINSTMSG	+=	" -> \"libzimpl.*\" is the path to the ZIMPL library, e.g., \"../../zimpl/lib/libzimpl.linux.x86.gnu.opt.a\""
endif

IPOPTDEP	:=	$(SRCDIR)/depend.ipopt
IPOPTSRC	:=	$(shell cat $(IPOPTDEP))
ifeq ($(IPOPT),true)
LINKER		=	CPP
FLAGS		+=	-DWITH_IPOPT -I$(LIBDIR)/ipoptinc $(IPOPT_FLAGS)
LDFLAGS		+=	$(LINKCXX_l)ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)$(LINKLIBSUFFIX) $(IPOPT_LDFLAGS)
ifeq ($(LIBEXT),$(STATICLIBEXT))
LDFLAGS		+=	`cat $(LIBDIR)/ipopt_addlibs.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT).$(STATICLIBEXT)`
else
LDFLAGS		+=	`cat $(LIBDIR)/ipopt_addlibs.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT).$(SHAREDLIBEXT)`
endif
SOFTLINKS	+=	$(LIBDIR)/ipoptinc
SOFTLINKS	+=	$(LIBDIR)/libipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT).$(SHAREDLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/ipopt_addlibs.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/ipopt_addlibs.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT).$(SHAREDLIBEXT)
LPIINSTMSG	+=	"\n  -> \"ipoptinc\" is a directory containing the ipopt header files, i.e., \"ipoptinc/IpIpoptApplication.hpp\" should exist.\n"
LPIINSTMSG	+=	" -> \"libipopt.*\" is the Ipopt library, e.g., \"/client/lib/ipopt-3.7.0/libipopt.a\"\n"
LPIINSTMSG	+=	" -> \"ipopt_addlibs.*\" is the Ipopt addlibs_cpp file, e.g., \"/client/lib/ipopt-3.7.0/ipopt_addlibs_cpp.txt\"\n"
endif

ifeq ($(READLINE),true)
FLAGS		+=	-DWITH_READLINE $(READLINE_FLAGS)
LDFLAGS		+=	$(READLINE_LDFLAGS)
endif

ifeq ($(GMP),true)
FLAGS		+=	-DWITH_GMP $(GMP_FLAGS)
LDFLAGS		+=	$(GMP_LDFLAGS)
endif

ifeq ($(ZLIB),true)
FLAGS		+=	-DWITH_ZLIB $(ZLIB_FLAGS)
LDFLAGS		+=	$(ZLIB_LDFLAGS)
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
			scip/nlpi.o \
			scip/nlpi_oracle.o \
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
			scip/cons_countsols.o \
			scip/cons_eqknapsack.o \
			scip/cons_indicator.o \
			scip/cons_integral.o \
			scip/cons_invarknapsack.o \
			scip/cons_knapsack.o \
			scip/cons_linear.o \
			scip/cons_logicor.o \
			scip/cons_or.o \
			scip/cons_orbitope.o \
			scip/cons_quadratic.o \
			scip/cons_setppc.o \
			scip/cons_soc.o \
			scip/cons_sos1.o \
			scip/cons_sos2.o \
			scip/cons_varbound.o \
			scip/cons_xor.o \
			scip/dialog_default.o \
			scip/disp_default.o \
			scip/heur_actconsdiving.o \
			scip/heur_coefdiving.o \
			scip/heur_crossover.o \
			scip/heur_dins.o \
			scip/heur_feaspump.o \
			scip/heur_fixandinfer.o \
			scip/heur_fracdiving.o \
			scip/heur_guideddiving.o \
			scip/heur_intdiving.o \
			scip/heur_intshifting.o \
			scip/heur_linesearchdiving.o \
			scip/heur_localbranching.o \
			scip/heur_mutation.o \
			scip/heur_nlp.o \
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
			scip/heur_trivial.o \
			scip/heur_trysol.o \
			scip/heur_twoopt.o \
			scip/heur_undercover.o \
			scip/heur_veclendiving.o \
			scip/heur_zirounding.o \
			scip/nodesel_bfs.o \
			scip/nodesel_dfs.o \
			scip/nodesel_estimate.o \
			scip/nodesel_hybridestim.o \
			scip/nodesel_restartdfs.o \
			scip/presol_boundshift.o \
			scip/presol_dualfix.o \
			scip/presol_implics.o \
			scip/presol_inttobinary.o \
			scip/presol_probing.o \
			scip/presol_trivial.o \
			scip/prop_pseudoobj.o \
			scip/prop_rootredcost.o \
			scip/prop_vbounds.o \
			scip/reader_ccg.o \
			scip/reader_cip.o \
			scip/reader_cnf.o \
			scip/reader_fix.o \
			scip/reader_fzn.o \
			scip/reader_gms.o \
			scip/reader_lp.o \
			scip/reader_mps.o \
			scip/reader_opb.o \
			scip/reader_ppm.o \
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
			scip/sepa_rapidlearning.o \
			scip/sepa_redcost.o \
			scip/sepa_strongcg.o \
			scip/sepa_zerohalf.o \
			tclique/tclique_branch.o \
			tclique/tclique_coloring.o \
			tclique/tclique_graph.o

ifeq ($(EXACTSOLVE),true)
SCIPLIBOBJ	+=	scip/cons_exactlp.o \
			scip/primalex.o \
			scip/solex.o \
			rectlu/rectlu_factor.o
endif

SCIPLIB		=	$(SCIPLIBNAME).$(BASE)
SCIPLIBFILE	=	$(LIBDIR)/lib$(SCIPLIB).$(LIBEXT)
SCIPLIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(SCIPLIBOBJ))
SCIPLIBSRC	=	$(addprefix $(SRCDIR)/,$(SCIPLIBOBJ:.o=.c))
SCIPLIBDEP	=	$(SRCDIR)/depend.sciplib.$(OPT)
SCIPLIBLINK	=	$(LIBDIR)/lib$(SCIPLIBSHORTNAME).$(BASE).$(LIBEXT)

ifeq ($(IPOPT),true)
SCIPLIBOBJ += scip/nlpi_ipopt.o
SCIPLIBOBJFILES += $(LIBOBJDIR)/scip/nlpi_ipopt.o
SCIPLIBSRC += $(SRCDIR)/scip/nlpi_ipopt.cpp
endif

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



LINKSMARKERFILE	=	$(LIBDIR)/linkscreated.$(LPS)-$(LPSOPT).$(LPSEX)-$(LPSEXOPT).$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX).$(ZIMPL)-$(ZIMPLOPT).$(IPOPT)-$(IPOPTOPT)
LASTSETTINGS	=	$(OBJDIR)/make.lastsettings

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(SCIPLIBFILE) $(OBJSCIPLIBFILE) $(LPILIBFILE) $(LPIEXLIBFILE) $(LPILIBLINK) $(LPIEXLIBLINK) $(SCIPLIBLINK) $(OBJSCIPLIBLINK) $(MAINLINK) $(MAINSHORTLINK) \
		$(LPILIBOBJFILES) $(LPIEXLIBOBJFILES) $(SCIPLIBOBJFILES) $(OBJSCIPLIBOBJFILES) $(MAINOBJFILES)
endif

.PHONY: all
all:            checklpsdefine checklpsexdefine $(LINKSMARKERFILE) $(SCIPLIBFILE) $(OBJSCIPLIBFILE) $(LPILIBFILE) $(LPIEXLIBFILE) $(MAINFILE) $(LPILIBLINK) $(LPIEXLIBLINK) $(SCIPLIBLINK) $(OBJSCIPLIBLINK) $(MAINLINK) $(MAINSHORTLINK)

.PHONY: lint
lint:		$(SCIPLIBSRC) $(OBJSCIPLIBSRC) $(LPILIBSRC) $(LPIEXLIBSRC) $(MAINSRC)
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
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS);

.PHONY: testcount
testcount:		
		cd check; \
		$(SHELL) ./checkcount.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS);

.PHONY: testcplex
testcplex:		
		cd check; \
		$(SHELL) ./check_cplex.sh $(TEST) $(CPLEX) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) 0.0 $(CONTINUE);

.PHONY: testmosek
testmosek:		
		cd check; \
		$(SHELL) ./check_mosek.sh $(TEST) $(MOSEK) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) 0.0 $(DISPFREQ) $(CONTINUE);

.PHONY: testcbc
testcbc:		
		cd check; \
		$(SHELL) ./check_cbc.sh $(TEST) $(CBC) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) 0.0 $(CONTINUE);

.PHONY: testcbcparallel
testcbcparallel:                
		cd check; \
                $(SHELL) ./check_cbc.sh $(TEST) $(CBCPARALLEL) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) 0.0 $(CONTINUE);

.PHONY: testcluster
testcluster:		
		cd check; \
		$(SHELL) ./check_cluster.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS);

.PHONY: testclustercbc
testclustercbc:		
		cd check; \
		$(SHELL) ./check_cluster_cbc.sh $(TEST) $(CBC) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(FEASTOL) 0.0 $(CONTINUE) $(LOCK);

.PHONY: testgurobi
testgurobi:		
		cd check; \
		$(SHELL) ./check_gurobi.sh $(TEST) $(GUROBI) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) 0.0 $(CONTINUE);

.PHONY: testglpk
testglpk:		
		cd check; \
		$(SHELL) ./check_glpk.sh $(TEST) $(GLPK) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) 0.0 $(CONTINUE);

.PHONY: testsymphony
testsymphony:		
		cd check; \
		$(SHELL) ./check_symphony.sh $(TEST) $(SYMPHONY) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) 0.0 $(CONTINUE);

.PHONY: testblis
testblis:		
		cd check; \
		$(SHELL) ./check_blis.sh $(TEST) $(BLIS) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) 0.0 $(CONTINUE);

.PHONY: tags
tags:
		rm -f TAGS; ctags -e -R -h ".c.cpp.h" --exclude=".*" src/; 

$(LPILIBLINK):	$(LPILIBFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(LPILIBFILE)) $(notdir $@)

$(LPIEXLIBLINK): $(LPIEXLIBFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(LPIEXLIBFILE)) $(notdir $@)

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
ifneq ($(LIBOBJDIR),)
		@cd $(LIBOBJDIR) && rm -f */*.o && rmdir $(LIBOBJSUBDIRS)
		@-rmdir $(LIBOBJDIR)
endif
ifneq ($(BINOBJDIR),)
		@-rm -f $(BINOBJDIR)/*.o && rmdir $(BINOBJDIR)
endif
ifneq ($(OBJDIR),)
		@-rm -f $(LASTSETTINGS)
		@-rmdir $(OBJDIR)
endif
		@echo "-> remove objective files"
		@-rm -f $(SCIPLIBFILE) $(OBJSCIPLIBFILE) $(LPILIBFILE) $(LPIEXLIBFILE) $(MAINFILE) \
		$(LPILIBLINK) $(LPIEXLIBLINK) $(SCIPLIBLINK) $(OBJSCIPLIBLINK) $(MAINLINK) $(MAINSHORTLINK)
		@echo "-> remove libraries"
		@echo "-> remove binary"

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

.PHONY: lpiexdepend
lpiexdepend:
ifneq ($(LPSEX),none)
ifeq ($(LINKER),C)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(LPIEXLIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(LPIEXLIBDEP)'
endif
ifeq ($(LINKER),CPP)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(LPIEXLIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(LPIEXLIBDEP)'
endif
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
depend:		lpidepend lpiexdepend maindepend
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(SCIPLIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(SCIPLIBDEP)'
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(OBJSCIPLIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(OBJSCIPLIBDEP)'
		@echo `grep -l "WITH_ZLIB" $(ALLSRC)` >$(ZLIBDEP)
		@echo `grep -l "WITH_GMP" $(ALLSRC)` >$(GMPDEP)
		@echo `grep -l "WITH_READLINE" $(ALLSRC)` >$(READLINEDEP)
		@echo `grep -l "WITH_ZIMPL" $(ALLSRC)` >$(ZIMPLDEP)
		@echo `grep -l "WITH_IPOPT" $(ALLSRC)` >$(IPOPTDEP)

-include	$(MAINDEP)
-include	$(SCIPLIBDEP)
-include	$(OBJSCIPLIBDEP)
-include 	$(LPILIBDEP)
-include 	$(LPIEXLIBDEP)

$(MAINFILE):	$(BINDIR) $(BINOBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(LPIEXLIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
ifeq ($(LINKER),C)
		$(LINKCC) $(MAINOBJFILES) \
		$(LINKCC_L)$(LIBDIR) $(LINKCC_l)$(SCIPLIB)$(LINKLIBSUFFIX) $(LINKCC_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCC_l)$(LPIEXLIB)$(LINKLIBSUFFIX) \
		$(OFLAGS) $(LPSLDFLAGS) $(LPSEXLDFLAGS) $(LDFLAGS) $(LINKCC_o)$@
endif
ifeq ($(LINKER),CPP)
		$(LINKCXX) $(MAINOBJFILES) \
		$(LINKCXX_L)$(LIBDIR) $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPIEXLIB)$(LINKLIBSUFFIX) \
		$(OFLAGS) $(LPSLDFLAGS) $(LPSEXLDFLAGS) $(LDFLAGS) $(LINKCXX_o)$@
endif

$(SCIPLIBFILE):	checklpsdefine checklpsexdefine $(LIBOBJSUBDIRS) $(LIBDIR) touchexternal $(SCIPLIBOBJFILES) 
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(SCIPLIBOBJFILES) 
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

$(OBJSCIPLIBFILE):	$(LIBOBJSUBDIRS) $(LIBDIR) $(OBJSCIPLIBOBJFILES) 
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(OBJSCIPLIBOBJFILES) 
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

$(LPILIBFILE):	$(LIBOBJSUBDIRS) $(LIBDIR) $(LPILIBOBJFILES)
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(LPILIBOBJFILES)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

$(LPIEXLIBFILE): $(LIBOBJSUBDIRS) $(LIBDIR) $(LPIEXLIBOBJFILES)
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(LPIEXLIBOBJFILES)
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
touchexternal:	$(ZLIBDEP) $(GMPDEP) $(READLINEDEP) $(ZIMPLDEP) $(IPOPTDEP)
ifneq ($(ZLIB),$(LAST_ZLIB))
		@-touch $(ZLIBSRC)
endif
ifneq ($(GMP),$(LAST_GMP))
		@-touch $(GMPSRC)
endif
ifneq ($(READLINE),$(LAST_READLINE))
		@-touch $(READLINESRC)
endif
ifneq ($(ZIMPL),$(LAST_ZIMPL))
		@-touch $(ZIMPLSRC)
endif
ifneq ($(IPOPT),$(LAST_IPOPT))
		@-touch $(IPOPTSRC)
endif
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_ZLIB=$(ZLIB)" >> $(LASTSETTINGS)
		@echo "LAST_GMP=$(GMP)" >> $(LASTSETTINGS)
		@echo "LAST_READLINE=$(READLINE)" >> $(LASTSETTINGS)
		@echo "LAST_ZIMPL=$(ZIMPL)" >> $(LASTSETTINGS)
		@echo "LAST_IPOPT=$(IPOPT)" >> $(LASTSETTINGS)

$(LINKSMARKERFILE):
		@$(MAKE) links

.PHONY: links
links:		echosoftlinks $(LIBDIR) $(DIRECTORIES) $(SOFTLINKS)
		@rm -f $(LINKSMARKERFILE)
		@echo "this is only a marker" > $(LINKSMARKERFILE)

.PHONY: echosoftlinks
echosoftlinks:
		@echo
		@echo "- Current settings: LPS=$(LPS) LPSEX=$(LPSEX) OSTYPE=$(OSTYPE) ARCH=$(ARCH) COMP=$(COMP) SUFFIX=$(LINKLIBSUFFIX) ZIMPL=$(ZIMPL) ZIMPLOPT=$(ZIMPLOPT) IPOPT=$(IPOPT) IPOPTOPT=$(IPOPTOPT)"
		@echo
		@echo "* SCIP needs some softlinks to external programs, in particular, LP-solvers."
		@echo "* Please insert the paths to the corresponding directories/libraries below."
		@echo "* The links will be installed in the 'lib' directory."
		@echo "* For more information and if you experience problems see the INSTALL file."
		@echo
		@echo -e $(LPIINSTMSG)

$(DIRECTORIES):
		@echo
		@echo "- creating directory \"$@\""
		@-mkdir -p $@



.PHONY: $(SOFTLINKS)
$(SOFTLINKS):
ifeq ($(MAKESOFTLINKS), true)
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
				cd $(LIBDIR) ; \
				eval $(READ) TARGET ; \
				cd .. ; \
				if test "$$TARGET" != "" ; \
				then \
					echo "-> creating softlink \"$@\" -> \"$$TARGET\"" ; \
					rm -f $@ ; \
					ln -s $$TARGET $@ ; \
				else \
					echo "* skipped creation of softlink \"$@\". Call \"make links\" if needed later." ; \
				fi ; \
				echo ; \
			fi'
endif

.PHONY: checklpsdefine
checklpsdefine:
ifeq ($(LPILIBOBJ),)
		$(error invalid LP solver selected: LPS=$(LPS). Possible options are: $(LPSOPTIONS))
endif

.PHONY: checklpsexdefine
checklpsexdefine:
ifeq ($(LPIEXUSED),)
		$(error invalid exact LP solver selected: LPSEX=$(LPSEX). Possible options are: $(LPSEXOPTIONS))
endif


# --- EOF ---------------------------------------------------------------------
# DO NOT DELETE
