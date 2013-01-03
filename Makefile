#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            *
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

#-----------------------------------------------------------------------------
# detect host architecture
#-----------------------------------------------------------------------------
include make/make.detecthost



#-----------------------------------------------------------------------------
# default settings
#-----------------------------------------------------------------------------

VERSION		:=	3.0.1

TIME     	=  	3600
NODES           =       2100000000
MEM		=	6144
THREADS         =       1
DISPFREQ	=	10000
FEASTOL		=	default
TEST		=	short
SETTINGS        =       default
CONTINUE	=	false
LOCK		=	false
VALGRIND	=	false

VERBOSE		=	false
OPT		=	opt
COMP		=	gnu
LPS		=	spx
STATICLIBEXT	=	a
SHAREDLIBEXT	=	so
LIBEXT		=	$(STATICLIBEXT)
LINKER  	=	C
SOFTLINKS	=

INSTALLDIR	=	

#will this be compiled for parascip, necessary for dbg-builds and cppad to make it threadsafe
PARASCIP	=	false

SHARED		=	false
MAKESOFTLINKS	=	true
READLINE	=	true
ZLIB		=	true
GMP		=	auto
ZIMPL		=	true
IPOPT		=	false
EXPRINT		=	cppad
LPSCHECK	=	false
LPSOPT		=	opt
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
LINKRPATH	=	-Wl,-rpath,
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
GAMS            =       gams


SHELL		= 	bash
READ		=	read -e
LN_s		= 	ln -s

FLAGS		=	-I$(SRCDIR) -DWITH_SCIPDEF
OFLAGS		=
CFLAGS		=	
CXXFLAGS	=	
LDFLAGS		=	$(LINKCC_l)m$(LINKLIBSUFFIX)
ARFLAGS		=	cr
DFLAGS		=	-MM

GCCWARN		=	-pedantic -Wno-long-long -Wall -W -Wpointer-arith -Wcast-align -Wwrite-strings -Wshadow \
			-Wno-unknown-pragmas -Wno-unused-parameter \
			-Wredundant-decls -Wdisabled-optimization \
			-Wsign-compare -Wstrict-prototypes \
			-Wmissing-declarations -Wmissing-prototypes # -Wdeclaration-after-statement

GXXWARN		=	-pedantic -Wno-long-long -Wall -W -Wpointer-arith -Wcast-align -Wwrite-strings -Wshadow \
			-Wno-unknown-pragmas -Wno-unused-parameter \
			-Wredundant-decls -Wdisabled-optimization \
			-Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder \
			-Woverloaded-virtual -Wsign-promo -Wsynth \
			-Wcast-qual -Wno-unused-parameter # -Wold-style-cast -Wshadow -Wundef

BASE		=	$(OSTYPE).$(ARCH).$(COMP).$(OPT)
OBJDIR		=	obj/O.$(BASE)
BINOBJDIR	=	$(OBJDIR)/bin
LIBOBJDIR	=	$(OBJDIR)/lib
LIBOBJSUBDIRS	=       scip objscip blockmemshell tclique nlpi xml dijkstra
SRCDIR		=	src
LIBDIR		=	lib
BINDIR		=	bin
INCLUDEDIR	=	include
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
# Memory Management
#-----------------------------------------------------------------------------

#FLAGS		+=	-DBMS_NOSAFEMEM
#FLAGS		+=	-DBMS_NOBLOCKMEM

#-----------------------------------------------------------------------------
# SHARED Libaries
#-----------------------------------------------------------------------------

ifeq ($(SHARED),true)
FLAGS		+=	-fPIC
LIBEXT		=	$(SHAREDLIBEXT)
LIBBUILD	=	$(LINKCC)
LIBBUILDFLAGS	+=      -shared
LIBBUILD_o	= 	-o # the trailing space is important
ARFLAGS		=
RANLIB		=
endif

#-----------------------------------------------------------------------------
# PARASCIP
#-----------------------------------------------------------------------------

ifeq ($(PARASCIP),false)
FLAGS		+=	-DNPARASCIP
else
LDFLAGS         +=      -lpthread
endif

#-----------------------------------------------------------------------------
# LP Solver Interface
#-----------------------------------------------------------------------------

LPILIBSHORTNAME	=	lpi$(LPS)
LPILIBNAME	=	$(LPILIBSHORTNAME)-$(VERSION)
LPILIBOBJ	=
LPSOPTIONS	=
LPIINSTMSG	=

LPSCHECKDEP	:=	$(SRCDIR)/depend.lpscheck
LPSCHECKSRC	:=	$(shell cat $(LPSCHECKDEP))

LPSOPTIONS	+=	cpx
ifeq ($(LPS),cpx)
FLAGS		+=	-I$(LIBDIR)/cpxinc
LPSLDFLAGS	+=	$(LINKCC_l)cplex.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCC_l)pthread$(LINKLIBSUFFIX)
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
LPSLDFLAGS	+=	$(LINKCC_l)xpress.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)  $(LINKCC_l)pthread$(LINKLIBSUFFIX)  $(LINKCC_l)dl$(LINKLIBSUFFIX)
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
LPSLDFLAGS	+=	$(LINKCC_l)mosek.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) \
			$(LINKCXX_l)iomp5.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCC_l)pthread$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_msk.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/mskinc
SOFTLINKS	+=	$(LIBDIR)/libmosek.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libmosek.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libiomp5.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libiomp5.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
LPIINSTMSG	=	"  -> \"mskinc\" is the path to the Mosek \"include\" directory, e.g., \"<Mosek-path>/include\".\n"
LPIINSTMSG	+=	" -> \"libmosek.*\" is the path to the Mosek library, e.g., \"<Mosek-path>/lib/libmosek.a\".\n"
LPIINSTMSG	+=	" -> \"libiomp5.*\" is the path to the libiomp5, e.g., \"<Mosek-path>/lib/libiomp5.a\""
endif

LPSOPTIONS	+=	spx
ifeq ($(LPS),spx)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/spxinc
LPSLDFLAGS	+=	$(LINKCXX_l)soplex.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)$(LINKLIBSUFFIX)
LPILIBOBJ	=	scip/lpi_spx.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/scip/lpi_spx.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/spxinc
SOFTLINKS	+=	$(LIBDIR)/libsoplex.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libsoplex.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(SHAREDLIBEXT)
LPIINSTMSG	=	"  -> \"spxinc\" is the path to the SoPlex \"src\" directory, e.g., \"../../soplex/src\".\n"
LPIINSTMSG	+=	" -> \"libsoplex.*\" is the path to the SoPlex library, e.g., \"../../soplex/lib/libsoplex.linux.x86.gnu.opt.a\""
ifeq ($(LPSCHECK),true)
FLAGS		+=	-DWITH_LPSCHECK -I$(LIBDIR)/cpxinc
LPSLDFLAGS	+=	$(LINKCC_l)cplex.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCC_l)pthread$(LINKLIBSUFFIX)
SOFTLINKS	+=	$(LIBDIR)/cpxinc
SOFTLINKS	+=	$(LIBDIR)/libcplex.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/libcplex.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
LPIINSTMSG	+=	"  -> \"cpxinc\" is the path to the CPLEX \"include\" directory, e.g., \"<CPLEX-path>/include/ilcplex\".\n"
LPIINSTMSG	+=	" -> \"libcplex.*\" is the path to the CPLEX library, e.g., \"<CPLEX-path>/lib/x86_rhel4.0_3.4/static_pic/libcplex.a\""
endif
endif

LPSOPTIONS	+=	spx132
ifeq ($(LPS),spx132)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/spx132inc
LPSLDFLAGS	+=	$(LINKCXX_l)soplex132.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)$(LINKLIBSUFFIX)
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
CLPDIR		= 	$(LIBDIR)/clp.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)
FLAGS		+=	-I$(CLPDIR)/include/coin
# for newer Clp versions all linker flags are in share/coin/doc/Clp/clp_addlibs.txt
LPSLDFLAGS	+=	$(shell test -e $(CLPDIR)/share/coin/doc/Clp/clp_addlibs.txt && cat $(CLPDIR)/share/coin/doc/Clp/clp_addlibs.txt)
# if we could not find clp_addlibs file, try to guess linker flags
ifeq ($(LPSLDFLAGS),)
LPSLDFLAGS	+=	$(LINKCXX_L)$(CLPDIR)/lib $(LINKCXX_l)Clp$(LINKLIBSUFFIX) \
			$(LINKCXX_l)CoinUtils$(LINKLIBSUFFIX) \
			$(LINKCXX_l)bz2$(LINKLIBSUFFIX) $(LINKCXX_l)lapack$(LINKLIBSUFFIX)
endif
# ensure that also shared libraries are found while running the binary
ifneq ($(LINKRPATH),)
CLPFULLPATH	:=	$(realpath $(CURDIR)/$(CLPDIR))
LPSLDFLAGS	+=	$(LINKRPATH)$(CLPFULLPATH)/lib
endif
LPILIBOBJ	=	scip/lpi_clp.o scip/bitencode.o blockmemshell/memory.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/scip/lpi_clp.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(CLPDIR)
LPIINSTMSG	=	"  -> \"clp.*\" is a directory containing the Clp installation, i.e., \"clp.*/include/coin/ClpModel.hpp\" should exist.\n"
endif

LPSOPTIONS	+=	qso
ifeq ($(LPS),qso)
FLAGS         	+=      -I$(LIBDIR)/qsinc
LPSLDFLAGS    	+=       $(LINKCC_l)qsopt.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCC_l)pthread$(LINKLIBSUFFIX)
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
LPSLDFLAGS	+=	$(LINKCC_l)gurobi.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX) $(LINKCC_l)pthread$(LINKLIBSUFFIX)
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
LPILIBOBJ	=	scip/lpi_none.o blockmemshell/memory.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
endif

LPILIB		=	$(LPILIBNAME).$(BASE)
LPILIBFILE	=	$(LIBDIR)/lib$(LPILIB).$(LIBEXT)
LPILIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(LPILIBOBJ))
LPILIBDEP	=	$(SRCDIR)/depend.lpilib.$(LPS).$(OPT)
LPILIBLINK	=	$(LIBDIR)/lib$(LPILIBSHORTNAME).$(BASE).$(LIBEXT)
LPILIBSHORTLINK = 	$(LIBDIR)/lib$(LPILIBSHORTNAME).$(LIBEXT)
ALLSRC		+=	$(LPILIBSRC)


#-----------------------------------------------------------------------------
# NLP Solver Interfaces and expression interpreter 
#-----------------------------------------------------------------------------

NLPILIBCOBJ	= nlpi/nlpi.o \
		  nlpi/nlpioracle.o \
		  nlpi/expr.o

NLPILIBCXXOBJ	= nlpi/intervalarithext.o

NLPILIBSCIPOBJ	= blockmemshell/memory.o \
		  scip/misc.o \
		  scip/intervalarith.o \
		  scip/interrupt.o \
		  scip/message.o

ifeq ($(EXPRINT),none)
NLPILIBCOBJ += 	nlpi/exprinterpret_none.o
endif
ifeq ($(EXPRINT),cppad)
NLPILIBCXXOBJ += nlpi/exprinterpret_cppad.o
NLPILIBSHORTNAMECPPAD = .cppad
endif

ifeq ($(IPOPT),true)
NLPILIBSHORTNAME = $(NLPILIBSHORTNAME).ipopt
NLPILIBSHORTNAMEIPOPT = .ipopt
NLPILIBCXXOBJ	+= nlpi/nlpi_ipopt.o
else
NLPILIBCOBJ	+= nlpi/nlpi_ipopt_dummy.o
endif

NLPILIBSHORTNAME = nlpi$(NLPILIBSHORTNAMECPPAD)$(NLPILIBSHORTNAMEIPOPT)
NLPILIBNAME	=	$(NLPILIBSHORTNAME)-$(VERSION)
NLPILIB		=	$(NLPILIBNAME).$(BASE)
NLPILIBFILE	=	$(LIBDIR)/lib$(NLPILIB).$(LIBEXT)
NLPILIBOBJFILES =	$(addprefix $(LIBOBJDIR)/,$(NLPILIBCOBJ)) $(addprefix $(LIBOBJDIR)/,$(NLPILIBCXXOBJ))
NLPILIBSCIPOBJFILES =	$(addprefix $(LIBOBJDIR)/,$(NLPILIBSCIPOBJ))
NLPILIBSRC	=	$(addprefix $(SRCDIR)/,$(NLPILIBCOBJ:.o=.c)) $(addprefix $(SRCDIR)/,$(NLPILIBCXXOBJ:.o=.cpp))
NLPILIBDEP	=	$(SRCDIR)/depend.nlpilib$(NLPILIBSHORTNAMECPPAD)$(NLPILIBSHORTNAMEIPOPT).$(OPT)
NLPILIBLINK	=	$(LIBDIR)/lib$(NLPILIBSHORTNAME).$(BASE).$(LIBEXT)
NLPILIBSHORTLINK	=	$(LIBDIR)/lib$(NLPILIBSHORTNAME).$(LIBEXT)
ALLSRC		+=	$(NLPILIBSRC)

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
LPIINSTMSG	+=	"\n  -> \"zimplinc\" is a directory containing the path to the ZIMPL \"src\" directory, e.g., \"../../../zimpl/src\".\n"
LPIINSTMSG	+=	" -> \"libzimpl.*\" is the path to the ZIMPL library, e.g., \"../../zimpl/lib/libzimpl.linux.x86.gnu.opt.a\""
endif

ifeq ($(IPOPT),true)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)/include/coin $(IPOPT_FLAGS)
# for Ipopt >= 3.9.0, all linker flags are in share/coin/doc/Ipopt/ipopt_addlibs_cpp.txt
# for Ipopt < 3.9.0, we need to link against libipopt from the lib directory, plus the additional flags given in share/doc/coin/Ipopt/ipopt_addlibs_cpp.txt
# finally, if no ipopt_addlibs_cpp.txt is present in the ipopt installation but a user put one into SCIP's libdir, take this one
LDFLAGS		+=	$(shell test -e $(LIBDIR)/ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)/share/coin/doc/Ipopt/ipopt_addlibs_cpp.txt && \
  cat $(LIBDIR)/ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)/share/coin/doc/Ipopt/ipopt_addlibs_cpp.txt)
LDFLAGS		+=	$(shell test -e $(LIBDIR)/ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)/share/doc/coin/Ipopt/ipopt_addlibs_cpp.txt && \
  (echo $(LINKCXX_L)$(LIBDIR)/ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)/lib $(LINKCXX_l)ipopt$(LINKLIBSUFFIX); cat $(LIBDIR)/ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)/share/doc/coin/Ipopt/ipopt_addlibs_cpp.txt))
LDFLAGS		+=	$(shell test -e $(LIBDIR)/ipopt_addlibs_cpp.txt && \
  (echo $(LINKCXX_L)$(LIBDIR)/ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)/lib $(LINKCXX_l)ipopt; cat $(LIBDIR)/ipopt_addlibs_cpp.txt))
# ensure that also shared libraries are found while running the binary
# for Ipopt 3.9.x, the libraries are installed in the lib/coin and lib/coin/ThirdParty subdirectories
# for Ipopt != 3.9.x, they are installed into the lib subdirectory
ifneq ($(LINKRPATH),)
IPOPTFULLPATH	:=	$(realpath $(CURDIR)/$(LIBDIR)/ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT))
LDFLAGS		+=	$(LINKRPATH)$(IPOPTFULLPATH)/lib
LDFLAGS		+=	$(shell test -e $(IPOPTFULLPATH)/lib/coin && echo $(LINKRPATH)$(IPOPTFULLPATH)/lib/coin)
LDFLAGS		+=	$(shell test -e $(IPOPTFULLPATH)/lib/coin/ThirdParty && echo $(LINKRPATH)$(IPOPTFULLPATH)/lib/coin/ThirdParty)
endif
SOFTLINKS	+=	$(LIBDIR)/ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)
LPIINSTMSG	+=	"\n  -> \"ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)\" is a directory containing the ipopt installation, i.e., \"ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)/include/coin/IpIpoptApplication.hpp\", \"ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)/lib/libipopt*\", ... should exist.\n"
endif

ifeq ($(EXPRINT),cppad)
LINKER		=	CPP
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
SCIPPLUGINLIBOBJ=       scip/branch_allfullstrong.o \
			scip/branch_fullstrong.o \
			scip/branch_inference.o \
			scip/branch_leastinf.o \
			scip/branch_mostinf.o \
			scip/branch_pscost.o \
			scip/branch_random.o \
			scip/branch_relpscost.o \
			scip/cons_abspower.o \
			scip/cons_and.o \
			scip/cons_bivariate.o \
			scip/cons_bounddisjunction.o \
			scip/cons_conjunction.o \
			scip/cons_countsols.o \
			scip/cons_cumulative.o \
			scip/cons_disjunction.o \
			scip/cons_indicator.o \
			scip/cons_integral.o \
			scip/cons_knapsack.o \
			scip/cons_linear.o \
			scip/cons_linking.o \
			scip/cons_logicor.o \
			scip/cons_nonlinear.o \
			scip/cons_or.o \
			scip/cons_orbitope.o \
			scip/cons_pseudoboolean.o \
			scip/cons_quadratic.o \
			scip/cons_setppc.o \
			scip/cons_soc.o \
			scip/cons_sos1.o \
			scip/cons_sos2.o \
			scip/cons_superindicator.o \
			scip/cons_varbound.o \
			scip/cons_xor.o \
			scip/dialog_default.o \
			scip/disp_default.o \
			scip/heur_actconsdiving.o \
			scip/heur_clique.o \
			scip/heur_coefdiving.o \
			scip/heur_crossover.o \
			scip/heur_dins.o \
			scip/heur_feaspump.o \
			scip/heur_fixandinfer.o \
			scip/heur_fracdiving.o \
			scip/heur_guideddiving.o \
			scip/heur_zeroobj.o \
			scip/heur_intdiving.o \
			scip/heur_intshifting.o \
			scip/heur_linesearchdiving.o \
			scip/heur_localbranching.o \
			scip/heur_mutation.o \
			scip/heur_nlpdiving.o \
			scip/heur_objpscostdiving.o \
			scip/heur_octane.o \
			scip/heur_oneopt.o \
			scip/heur_pscostdiving.o \
			scip/heur_rens.o \
			scip/heur_rins.o \
			scip/heur_rootsoldiving.o \
			scip/heur_rounding.o \
			scip/heur_shiftandpropagate.o \
			scip/heur_shifting.o \
			scip/heur_simplerounding.o \
			scip/heur_subnlp.o \
			scip/heur_trivial.o \
			scip/heur_trysol.o \
			scip/heur_twoopt.o \
			scip/heur_undercover.o \
			scip/heur_vbounds.o \
			scip/heur_veclendiving.o \
			scip/heur_zirounding.o \
			scip/message_default.o \
			scip/nodesel_bfs.o \
			scip/nodesel_dfs.o \
			scip/nodesel_estimate.o \
			scip/nodesel_hybridestim.o \
			scip/nodesel_restartdfs.o \
			scip/presol_boundshift.o \
			scip/presol_components.o \
			scip/presol_convertinttobin.o \
			scip/presol_domcol.o\
			scip/presol_dualfix.o \
			scip/presol_gateextraction.o \
			scip/presol_implics.o \
			scip/presol_inttobinary.o \
			scip/presol_trivial.o \
			scip/prop_genvbounds.o \
			scip/prop_obbt.o \
			scip/prop_probing.o \
			scip/prop_pseudoobj.o \
			scip/prop_redcost.o \
			scip/prop_rootredcost.o \
			scip/prop_vbounds.o \
			scip/reader_bnd.o \
			scip/reader_ccg.o \
			scip/reader_cip.o \
			scip/reader_cnf.o \
			scip/reader_fix.o \
			scip/reader_fzn.o \
			scip/reader_gms.o \
			scip/reader_lp.o \
			scip/reader_mps.o \
			scip/reader_opb.o \
			scip/reader_osil.o \
			scip/reader_pip.o \
			scip/reader_ppm.o \
			scip/reader_rlp.o \
			scip/reader_sol.o \
			scip/reader_wbo.o \
			scip/reader_zpl.o \
			scip/sepa_cgmip.o \
			scip/sepa_clique.o \
			scip/sepa_closecuts.o \
			scip/sepa_cmir.o \
			scip/sepa_flowcover.o \
			scip/sepa_gomory.o \
			scip/sepa_impliedbounds.o \
			scip/sepa_intobj.o \
			scip/sepa_mcf.o \
			scip/sepa_oddcycle.o \
			scip/sepa_rapidlearning.o \
			scip/sepa_strongcg.o \
			scip/sepa_zerohalf.o

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
			scip/nlp.o \
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
			scip/scipgithash.o \
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
			tclique/tclique_branch.o \
			tclique/tclique_coloring.o \
			tclique/tclique_graph.o \
			dijkstra/dijkstra.o \
			xml/xmlparse.o

SCIPLIB		=	$(SCIPLIBNAME).$(BASE)
SCIPLIBFILE	=	$(LIBDIR)/lib$(SCIPLIB).$(LIBEXT)
SCIPLIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(SCIPPLUGINLIBOBJ))
SCIPLIBOBJFILES	+=	$(addprefix $(LIBOBJDIR)/,$(SCIPLIBOBJ))
SCIPLIBSRC	=	$(addprefix $(SRCDIR)/,$(SCIPPLUGINLIBOBJ:.o=.c))
SCIPLIBSRC	+=	$(addprefix $(SRCDIR)/,$(SCIPLIBOBJ:.o=.c))
SCIPPLUGININCSRC=	$(addprefix $(SRCDIR)/,$(SCIPPLUGINLIBOBJ:.o=.h))
SCIPLIBDEP	=	$(SRCDIR)/depend.sciplib.$(OPT)
SCIPLIBLINK	=	$(LIBDIR)/lib$(SCIPLIBSHORTNAME).$(BASE).$(LIBEXT)
SCIPLIBSHORTLINK = 	$(LIBDIR)/lib$(SCIPLIBSHORTNAME).$(LIBEXT)

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
OBJSCIPINCSRC	=	$(addprefix $(SRCDIR)/,$(OBJSCIPLIBOBJ:.o=.h))
OBJSCIPLIBDEP	=	$(SRCDIR)/depend.objsciplib.$(OPT)
OBJSCIPLIBLINK	=	$(LIBDIR)/lib$(OBJSCIPLIBSHORTNAME).$(BASE).$(LIBEXT)
OBJSCIPLIBSHORTLINK=	$(LIBDIR)/lib$(OBJSCIPLIBSHORTNAME).$(LIBEXT)
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

ifneq ($(LINKRPATH),)
LDFLAGS		+=	$(LINKRPATH)$(SCIPDIR)/$(LIBDIR)
endif


LINKSMARKERFILE	=	$(LIBDIR)/linkscreated.$(LPS)-$(LPSOPT).$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX).$(ZIMPL)-$(ZIMPLOPT).$(IPOPT)-$(IPOPTOPT)
LASTSETTINGS	=	$(OBJDIR)/make.lastsettings

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(SCIPLIBFILE) $(OBJSCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) \
		$(LPILIBLINK) $(LPILIBSHORTLINK) $(SCIPLIBLINK) $(SCIPLIBSHORTLINK) \
		$(OBJSCIPLIBLINK) $(OBJSCIPLIBSHORTLINK) $(NLPILIBLINK) $(NLPILIBSHORTLINK) \
		$(MAINLINK) $(MAINSHORTLINK) \
		$(LPILIBOBJFILES) $(NLPILIBOBJFILES) $(SCIPLIBOBJFILES) $(OBJSCIPLIBOBJFILES) $(MAINOBJFILES)
endif

all: 		githash libs $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)

libs: 		$(LINKSMARKERFILE) $(SCIPLIBFILE) $(OBJSCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(LPILIBLINK) $(LPILIBSHORTLINK) $(NLPILIBLINK) $(NLPILIBSHORTLINK) $(SCIPLIBLINK) $(SCIPLIBSHORTLINK) $(OBJSCIPLIBLINK) $(OBJSCIPLIBSHORTLINK) 

.PHONY: lint
lint:		$(SCIPLIBSRC) $(OBJSCIPLIBSRC) $(LPILIBSRC) $(NLPILIBSRC) $(MAINSRC)
		-rm -f lint.out
ifeq ($(FILES),)
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) lint/$(MAINSHORTNAME).lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'
else
		$(SHELL) -ec  'for i in $(FILES); \
			do \
			echo $$i; \
			$(LINT) lint/$(MAINSHORTNAME).lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'
endif

.PHONY: doc
doc: 		
		cd doc; $(DOXY) $(MAINSHORTNAME).dxy ; $(DOXY) $(MAINSHORTNAME)devel.dxy


.PHONY: check
check:		test

.PHONY: test
test:
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(VALGRIND);

.PHONY: testcount
testcount:		
		cd check; \
		$(SHELL) ./check_count.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS);

.PHONY: testcplex
testcplex:		
		cd check; \
		$(SHELL) ./check_cplex.sh $(TEST) $(CPLEX) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE);

.PHONY: testmosek
testmosek:		
		cd check; \
		$(SHELL) ./check_mosek.sh $(TEST) $(MOSEK) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE);

.PHONY: testcbc
testcbc:		
		cd check; \
		$(SHELL) ./check_cbc.sh $(TEST) $(CBC) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(CONTINUE);

.PHONY: testcbcparallel
testcbcparallel:                
		cd check; \
                $(SHELL) ./check_cbc.sh $(TEST) $(CBCPARALLEL) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(CONTINUE);

.PHONY: testgurobi
testgurobi:		
		cd check; \
		$(SHELL) ./check_gurobi.sh $(TEST) $(GUROBI) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE);

.PHONY: testglpk
testglpk:		
		cd check; \
		$(SHELL) ./check_glpk.sh $(TEST) $(GLPK) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE);

.PHONY: testsymphony
testsymphony:		
		cd check; \
		$(SHELL) ./check_symphony.sh $(TEST) $(SYMPHONY) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE);

.PHONY: testblis
testblis:		
		cd check; \
		$(SHELL) ./check_blis.sh $(TEST) $(BLIS) $(SETTINGS) $(OSTYPE).$(ARCH).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE);

.PHONY: tags
tags:
		rm -f TAGS; ctags -e -R -h ".c.cpp.h" --exclude=".*" src/; 

# include target to detect the current git hash 
-include make/local/make.detectgithash

# this empty target is needed for the SCIP release versions
githash::	# do not remove the double-colon

# include local targets 
-include make/local/make.targets

# include install/uninstall targets
-include make/make.install

# the testgams target need to come after make/local/make.targets has been included (if any), because the latter may assign a value to CLIENTTMPDIR
.PHONY: testgams
testgams:
ifeq ($(CLIENTTMPDIR),)
		CLIENTTMPDIR=/tmp
endif
		cd check; \
		$(SHELL) ./check_gamscluster.sh $(TEST) $(GAMS) "$(GAMSSOLVER)" $(SETTINGS) $(OSTYPE).$(ARCH) $(TIME) $(NODES) "$(GAP)" $(THREADS) $(CONTINUE) "$(CONVERTSCIP)" local dummy dummy $(CLIENTTMPDIR) 1 true;

$(LPILIBLINK):	$(LPILIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(LPILIBFILE)) $(notdir $@)

$(LPILIBSHORTLINK):	$(LPILIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(LPILIBFILE)) $(notdir $@)

$(NLPILIBLINK):	$(NLPILIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(NLPILIBFILE)) $(notdir $@)

$(NLPILIBSHORTLINK):	$(NLPILIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(NLPILIBFILE)) $(notdir $@)

$(SCIPLIBLINK):	$(SCIPLIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(SCIPLIBFILE)) $(notdir $@)

$(SCIPLIBSHORTLINK):	$(SCIPLIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(SCIPLIBFILE)) $(notdir $@)

$(OBJSCIPLIBLINK):	$(OBJSCIPLIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(OBJSCIPLIBFILE)) $(notdir $@)

$(OBJSCIPLIBSHORTLINK):	$(OBJSCIPLIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(OBJSCIPLIBFILE)) $(notdir $@)

$(MAINLINK) $(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(MAINFILE)) $(notdir $@)

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
clean:          cleanlib cleanbin $(LIBOBJSUBDIRS) $(LIBOBJDIR) $(BINOBJDIR) $(OBJDIR)
ifneq ($(LIBOBJDIR),)
		@-(cd $(LIBOBJDIR) && rm -f */*.o && rmdir $(LIBOBJSUBDIRS));
		@-rmdir $(LIBOBJDIR)
endif
ifneq ($(BINOBJDIR),)
		@-rm -f $(BINOBJDIR)/*.o && rmdir $(BINOBJDIR)
endif
ifneq ($(OBJDIR),)
		@-rm -f $(LASTSETTINGS)
		@-rmdir $(OBJDIR)
endif

.PHONY: cleanlib
cleanlib:       $(LIBDIR)
		@echo "-> remove library $(SCIPLIBFILE)"
		@-rm -f $(SCIPLIBFILE) $(SCIPLIBLINK) $(SCIPLIBSHORTLINK)
		@echo "-> remove library $(OBJSCIPLIBFILE)"
		@-rm -f $(OBJSCIPLIBFILE) $(OBJSCIPLIBLINK) $(OBJSCIPLIBSHORTLINK)
		@echo "-> remove library $(LPILIBFILE)"
		@-rm -f $(LPILIBFILE) $(LPILIBLINK) $(LPILIBSHORTLINK)
		@echo "-> remove library $(NLPILIBFILE)"
		@-rm -f $(NLPILIBFILE) $(NLPILIBLINK) $(NLPILIBSHORTLINK)

.PHONY: cleanbin
cleanbin:       $(BINDIR) 
		@echo "-> remove binary $(MAINFILE)"
		@-rm -f $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)

.PHONY: lpidepend
lpidepend:
ifeq ($(LINKER),C)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(LPILIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		| sed '\''s|$(LIBDIR)/cpxinc/cpxconst.h||g'\'' \
		>$(LPILIBDEP)'
endif
ifeq ($(LINKER),CPP)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(LPILIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		| sed '\''s|$(LIBDIR)/clp[^ ]*||g'\'' \
		| sed '\''s|$(LIBDIR)/cpxinc/cpxconst.h||g'\'' \
		| sed '\''s|$(LIBDIR)/spxinc[^ ]*||g'\'' \
		>$(LPILIBDEP)'
endif
		@#we explicitely add all lpi's here, since the content of depend.lpscheck should be independent of the currently selected LPI, but contain all LPI's that use the WITH_LPSCHECK define
		@echo `grep -l "WITH_LPSCHECK" $(SCIPLIBSRC) $(OBJSCIPLIBSRC) $(MAINSRC) $(NLPILIBSRC) src/scip/lpi*.{c,cpp}` >$(LPSCHECKDEP)

.PHONY: nlpidepend
nlpidepend:
ifeq ($(LINKER),C)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(NLPILIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(NLPILIBDEP)'
endif
ifeq ($(LINKER),CPP)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(NLPILIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		| sed '\''s|$(LIBDIR)/ipopt[^ ]*||g'\'' \
		>$(NLPILIBDEP)'
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

.PHONY: scipdepend
scipdepend:
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(SCIPLIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(SCIPLIBDEP)'
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(OBJSCIPLIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(OBJSCIPLIBDEP)'
		@echo `grep -l "WITH_ZLIB" $(ALLSRC)` >$(ZLIBDEP)
		@echo `grep -l "WITH_GMP" $(ALLSRC)` >$(GMPDEP)
		@echo `grep -l "WITH_READLINE" $(ALLSRC)` >$(READLINEDEP)
		@echo `grep -l "WITH_ZIMPL" $(ALLSRC)` >$(ZIMPLDEP)

depend:		scipdepend lpidepend nlpidepend maindepend

-include	$(MAINDEP)
-include	$(SCIPLIBDEP)
-include	$(OBJSCIPLIBDEP)
-include	$(LPILIBDEP)
-include	$(NLPILIBDEP)

$(MAINFILE):	$(BINDIR) $(BINOBJDIR) $(LIBOBJSUBDIRS) $(SCIPLIBOBJFILES) $(LPILIBOBJFILES) $(NLPILIBOBJFILES) $(MAINOBJFILES)
		@echo "-> linking $@"
ifeq ($(LINKER),C)
		$(LINKCC) $(MAINOBJFILES) \
		$(LINKCC_L)$(LIBDIR) $(SCIPLIBOBJFILES) $(LPILIBOBJFILES) $(NLPILIBOBJFILES) \
		$(OFLAGS) $(LPSLDFLAGS) $(LDFLAGS) $(LINKCC_o)$@
endif
ifeq ($(LINKER),CPP)
		$(LINKCXX) $(MAINOBJFILES) \
		$(LINKCXX_L)$(LIBDIR) $(SCIPLIBOBJFILES) $(LPILIBOBJFILES) $(NLPILIBOBJFILES) \
		$(OFLAGS) $(LPSLDFLAGS) $(LDFLAGS) $(LINKCXX_o)$@
endif

$(SCIPLIBFILE):	checklpsdefine $(LIBDIR) $(LIBOBJSUBDIRS) touchexternal $(SCIPLIBOBJFILES)
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

$(NLPILIBFILE):	$(LIBOBJSUBDIRS) $(LIBDIR) $(NLPILIBOBJFILES) $(NLPILIBSCIPOBJFILES) 
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(NLPILIBOBJFILES) $(NLPILIBSCIPOBJFILES) 
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
touchexternal:	$(ZLIBDEP) $(GMPDEP) $(READLINEDEP) $(ZIMPLDEP) $(LPSCHECKDEP)
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
ifneq ($(LPSCHECK),$(LAST_LPSCHECK))
		@-touch $(LPSCHECKSRC)
endif
ifneq ($(SHARED),$(LAST_SHARED))
		@-touch $(ALLSRC)
endif
ifneq ($(USRCFLAGS),$(LAST_USRCFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRCXXFLAGS),$(LAST_USRCXXFLAGS))
		@-touch $(ALLSRC)
endif
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_ZLIB=$(ZLIB)" >> $(LASTSETTINGS)
		@echo "LAST_GMP=$(GMP)" >> $(LASTSETTINGS)
		@echo "LAST_READLINE=$(READLINE)" >> $(LASTSETTINGS)
		@echo "LAST_ZIMPL=$(ZIMPL)" >> $(LASTSETTINGS)
		@echo "LAST_LPSCHECK=$(LPSCHECK)" >> $(LASTSETTINGS)
		@echo "LAST_SHARED=$(SHARED)" >> $(LASTSETTINGS)
		@echo "LAST_USRCFLAGS=$(USRCFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRCXXFLAGS=$(USRCXXFLAGS)" >> $(LASTSETTINGS)

$(LINKSMARKERFILE):
		@$(MAKE) links

.PHONY: links
links:		echosoftlinks $(LIBDIR) $(DIRECTORIES) $(SOFTLINKS)
		@rm -f $(LINKSMARKERFILE)
		@echo "this is only a marker" > $(LINKSMARKERFILE)

.PHONY: echosoftlinks
echosoftlinks:
		@echo
		@echo "- Current settings: LPS=$(LPS) OSTYPE=$(OSTYPE) ARCH=$(ARCH) COMP=$(COMP) SUFFIX=$(LINKLIBSUFFIX) ZIMPL=$(ZIMPL) ZIMPLOPT=$(ZIMPLOPT) IPOPT=$(IPOPT) IPOPTOPT=$(IPOPTOPT) EXPRINT=$(EXPRINT)"
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
				cd $$DIRNAME ; \
				eval $(READ) TARGET ; \
				cd $(SCIPDIR) ; \
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
endif

.PHONY: checklpsdefine
checklpsdefine:
ifeq ($(LPILIBOBJ),)
		$(error invalid LP solver selected: LPS=$(LPS). Possible options are: $(LPSOPTIONS))
endif

# --- EOF ---------------------------------------------------------------------
# DO NOT DELETE
