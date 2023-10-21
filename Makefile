#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#
#*                                                                           *#
#*                  This file is part of the program and library             *#
#*         SCIP --- Solving Constraint Integer Programs                      *#
#*                                                                           *#
#*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      *#
#*                                                                           *#
#*  Licensed under the Apache License, Version 2.0 (the "License");          *#
#*  you may not use this file except in compliance with the License.         *#
#*  You may obtain a copy of the License at                                  *#
#*                                                                           *#
#*      http://www.apache.org/licenses/LICENSE-2.0                           *#
#*                                                                           *#
#*  Unless required by applicable law or agreed to in writing, software      *#
#*  distributed under the License is distributed on an "AS IS" BASIS,        *#
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *#
#*  See the License for the specific language governing permissions and      *#
#*  limitations under the License.                                           *#
#*                                                                           *#
#*  You should have received a copy of the Apache-2.0 license                *#
#*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         *#
#*                                                                           *#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *#

#@file    Makefile
#@brief   SCIP Makefile
#@author  Thorsten Koch
#@author  Tobias Achterberg
#@author  Marc Pfetsch
#@author  Timo Berthold

#-----------------------------------------------------------------------------
# paths variables
#-----------------------------------------------------------------------------

SCIPDIR		=	./

INSTALLDIR	=

#-----------------------------------------------------------------------------
# include make.project file
#-----------------------------------------------------------------------------

# do not use other open source projects; needs to be set before including make.project
ifeq ($(OPENSOURCE),false)
	override EXPRINT	=	none
	override GMP		=	false
	override READLINE	=	false
	override ZLIB		=	false
	override ZIMPL		=	false
	override IPOPT		=	false
	override AMPL		=	false
endif

# mark that this is a SCIP internal makefile
SCIPINTERNAL	=	true

# use THREADSAFE=true if compiled with TPI not equal to none
ifeq ($(TPI),omp)
	override THREADSAFE      =       true
endif
ifeq ($(TPI),tny)
        override THREADSAFE      =       true
endif


# load default settings and detect host architecture
include $(SCIPDIR)/make/make.project

#-----------------------------------------------------------------------------
# default settings
#-----------------------------------------------------------------------------

VERSION		=	$(SCIP_VERSION)
SCIPGITHASH	=
SOFTLINKS	=
MAKESOFTLINKS	=	true
TOUCHLINKS	=	false

#-----------------------------------------------------------------------------
# define build flags
#-----------------------------------------------------------------------------
BUILDFLAGS =	" ARCH=$(ARCH)\\n\
		COMP=$(COMP)\\n\
		DEBUGSOL=$(DEBUGSOL)\\n\
		EXPRINT=$(EXPRINT)\\n\
		GMP=$(GMP)\\n\
		IPOPT=$(IPOPT)\\n\
		IPOPTOPT=$(IPOPTOPT)\\n\
		LPS=$(LPS)\\n\
		LPSCHECK=$(LPSCHECK)\\n\
		LPSOPT=$(LPSOPT)\\n\
		NOBLKBUFMEM=$(NOBLKBUFMEM)\\n\
		NOBLKMEM=$(NOBLKMEM)\\n\
		NOBUFMEM=$(NOBUFMEM)\\n\
		OPT=$(OPT)\\n\
		OSTYPE=$(OSTYPE)\\n\
		THREADSAFE=$(THREADSAFE)\\n\
		PAPILO=$(PAPILO)\\n\
		READLINE=$(READLINE)\\n\
		SANITIZE=$(SANITIZE)\\n\
		SHARED=$(SHARED)\\n\
		SYM=$(SYM)\\n\
		USRARFLAGS=$(USRARFLAGS)\\n\
		USRCFLAGS=$(USRCFLAGS)\\n\
		USRCXXFLAGS=$(USRCXXFLAGS)\\n\
		USRFLAGS=$(USRFLAGS)\\n\
		USRLDFLAGS=$(USRLDFLAGS)\\n\
		USROFLAGS=$(USROFLAGS)\\n\
		VERSION=$(VERSION)\\n\
		WORHP=$(WORHP)\\n\
		WORHPOPT=$(WORHPOPT)\\n\
		ZIMPL=$(ZIMPL)\\n\
		ZIMPLOPT=$(ZIMPLOPT)\\n\
		AMPL=$(AMPL)\\n\
		ZLIB=$(ZLIB)"

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
FLAGS		+=	-I$(LIBDIR)/include/cpxinc
LPILIBOBJ	=	lpi/lpi_cpx.o scip/bitencode.o blockmemshell/memory.o scip/rbtree.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/include/cpxinc
ifeq ($(SHARED),true)
SOFTLINKS	+=	$(LIBDIR)/shared/libcplex.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
else
SOFTLINKS	+=	$(LIBDIR)/static/libcplex.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
endif
LPIINSTMSG	=	"  -> \"cpxinc\" is the path to the CPLEX \"include\" directory, e.g., \"<CPLEX-path>/include/ilcplex\".\n"
LPIINSTMSG	+=	" -> \"libcplex.*.a\" is the path to the CPLEX library, e.g., \"<CPLEX-path>/lib/x86-64_linux/static_pic/libcplex.a\"\n"
LPIINSTMSG	+=	" -> \"libcplex.*.so\" is the path to the CPLEX library, e.g., \"<CPLEX-path>/bin/x86-64_linux/libcplex1263.so\""
endif

# XPRESS is only available as shared library
LPSOPTIONS	+=	xprs
ifeq ($(LPS),xprs)
FLAGS		+=	-I$(LIBDIR)/include/xprsinc
LPILIBOBJ	=	lpi/lpi_xprs.o scip/bitencode.o blockmemshell/memory.o scip/rbtree.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/include/xprsinc
SOFTLINKS	+=	$(LIBDIR)/shared/libxprs.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
LPIINSTMSG	=	"  -> \"xprsinc\" is the path to the XPRESS \"include\" directory, e.g., \"<XPRESS-path>/include\".\n"
LPIINSTMSG	+=	" -> \"libxprs.*\" is the path to the XPRESS library, e.g., \"<XPRESS-path>/lib/libxprs.so\""
endif

# mosek only supports shared libraries
LPSOPTIONS	+=	msk
ifeq ($(LPS),msk)
FLAGS		+=	-I$(LIBDIR)/include/mskinc
LPILIBOBJ	=	lpi/lpi_msk.o scip/bitencode.o blockmemshell/memory.o scip/rbtree.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/include/mskinc
SOFTLINKS	+=	$(LIBDIR)/shared/libmosek.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/shared/libiomp5.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)    # for Mosek < 9
SOFTLINKS	+=	$(LIBDIR)/shared/libcilkrts.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)  # for Mosek < 10
LPIINSTMSG	=	"  -> \"mskinc\" is the path to the Mosek \"include\" directory, e.g., \"<Mosek-path>/tools/platform/linux64x86/h\".\n"
LPIINSTMSG	+=	" -> \"libmosek.*\" is the path to the Mosek library, e.g., \"<Mosek-path>/tools/platform/linux64x86/bin/libmosek64.so\".\n"
LPIINSTMSG	+=	" -> \"libiomp5.*\" is the path to the Intel OpenMP library, e.g., \"<Mosek-path>/tools/platform/linux64x86/bin/libiomp5.so\" (required for Mosek < 9.0 only).\n"
LPIINSTMSG	+=	" -> \"libcilkrts.*\" is the path to the cilk library, e.g., \"<Mosek-path>/tools/platform/linux64x86/bin/libcilkrts.so.5\" (required for Mosek < 10.0 only).\n"
endif

LPSOPTIONS	+=	spx1
ifeq ($(LPS),spx1)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/include/spxinc
LPILIBOBJ	=	lpi/lpi_spx1.o scip/bitencode.o blockmemshell/memory.o scip/rbtree.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/lpi/lpi_spx1.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/rbtree.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/include/spxinc
ifeq ($(SHARED),true)
SOFTLINKS	+=	$(LIBDIR)/shared/libsoplex.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(SHAREDLIBEXT)
else
SOFTLINKS	+=	$(LIBDIR)/static/libsoplex.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(STATICLIBEXT)
endif
LPIINSTMSG	=	"  -> \"spxinc\" is the path to the SoPlex \"src\" directory, e.g., \"<SoPlex-path>/src\".\n"
LPIINSTMSG	+=	" -> \"libsoplex.*\" is the path to the SoPlex library, e.g., \"<SoPlex-path>/lib/libsoplex.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(STATICLIBEXT)\""
ifeq ($(LPSCHECK),true)
FLAGS		+=	-DSCIP_WITH_LPSCHECK -I$(LIBDIR)/include/cpxinc
SOFTLINKS	+=	$(LIBDIR)/include/cpxinc
ifeq ($(SHARED),true)
SOFTLINKS	+=	$(LIBDIR)/shared/libcplex.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
else
SOFTLINKS	+=	$(LIBDIR)/static/libcplex.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
endif
LPIINSTMSG	+=	"  -> \"cpxinc\" is the path to the CPLEX \"include\" directory, e.g., \"<CPLEX-path>/include/ilcplex\".\n"
LPIINSTMSG	+=	" -> \"libcplex.*.a\" is the path to the CPLEX library, e.g., \"<CPLEX-path>/lib/x86-64_linux/static_pic/libcplex.a\"\n"
LPIINSTMSG	+=	" -> \"libcplex.*.so\" is the path to the CPLEX library, e.g., \"<CPLEX-path>/bin/x86-64_linux/libcplex1263.so\""
endif
endif

LPSOPTIONS	+=	spx ( = spx2)
ifeq ($(LPS),spx2)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/include/spxinc
LPILIBOBJ	=	lpi/lpi_spx2.o scip/bitencode.o blockmemshell/memory.o scip/rbtree.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/lpi/lpi_spx2.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/rbtree.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/include/spxinc
ifeq ($(SHARED),true)
SOFTLINKS	+=	$(LIBDIR)/shared/libsoplex.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(SHAREDLIBEXT)
else
SOFTLINKS	+=	$(LIBDIR)/static/libsoplex.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT).$(STATICLIBEXT)
endif
LPIINSTMSG	=	"  -> \"spxinc\" is the path to the SoPlex \"src\" directory, e.g., \"<SoPlex-path>/src\".\n"
LPIINSTMSG	+=	" -> \"libsoplex.*\" is the path to the SoPlex library, e.g., \"<SoPlex-path>/lib/libsoplex.linux.x86.gnu.opt.a\""
ifeq ($(LPSCHECK),true)
FLAGS		+=	-DSCIP_WITH_LPSCHECK -I$(LIBDIR)/include/cpxinc
SOFTLINKS	+=	$(LIBDIR)/include/cpxinc
ifeq ($(SHARED),true)
SOFTLINKS	+=	$(LIBDIR)/shared/libcplex.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
else
SOFTLINKS	+=	$(LIBDIR)/static/libcplex.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
endif
LPIINSTMSG	+=	"  -> \"cpxinc\" is the path to the CPLEX \"include\" directory, e.g., \"<CPLEX-path>/include/ilcplex\".\n"
LPIINSTMSG	+=	" -> \"libcplex.*.a\" is the path to the CPLEX library, e.g., \"<CPLEX-path>/lib/x86_rhel4.0_3.4/static_pic/libcplex.a\"\n"
LPIINSTMSG	+=	" -> \"libcplex.*.so\" is the path to the CPLEX library, e.g., \"<CPLEX-path>/bin/x86-64_linux/libcplex1263.so\""
endif
endif

ifeq ($(DEBUGSOL),true)
FLAGS		+=	-DWITH_DEBUG_SOLUTION
endif

LPSOPTIONS	+=	clp
ifeq ($(LPS),clp)
LINKER		=	CPP
FLAGS		+=	-I$(LIBDIR)/$(LIBTYPE)/clp.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)/include/coin -I$(LIBDIR)/$(LIBTYPE)/clp.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)/include/coin-or
LPILIBOBJ	=	lpi/lpi_clp.o scip/bitencode.o blockmemshell/memory.o scip/rbtree.o scip/message.o
LPILIBSRC	=	$(SRCDIR)/lpi/lpi_clp.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/blockmemshell/memory.c $(SRCDIR)/scip/rbtree.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/$(LIBTYPE)/clp.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)
LPIINSTMSG	=	"  -> \"clp.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)\" is the path to the Clp installation directory, i.e., \"<Clp-path>/include/coin/ClpModel.hpp\" should exist.\n"
endif

LPSOPTIONS	+=	qso
ifeq ($(LPS),qso)
FLAGS         	+=      -I$(LIBDIR)/include/qsinc
LPILIBOBJ     	= 	lpi/lpi_qso.o scip/bitencode.o blockmemshell/memory.o scip/rbtree.o scip/message.o
LPILIBSRC     	=       $(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS     	+=      $(LIBDIR)/include/qsinc
SOFTLINKS     	+=      $(LIBDIR)/static/libqsopt.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
LPIINSTMSG	=	"  -> \"qsinc\" is the path to the QSopt \"include\" directory, e.g., \"<QSopt-path>\".\n"
LPIINSTMSG	+=	" -> \"libqsopt.*\" is the path to the QSopt library, e.g., \"<QSopt-path>/libqsopt.a\""
endif

# Gurobi only supports shared libraries
LPSOPTIONS	+=	grb
ifeq ($(LPS),grb)
FLAGS		+=	-I$(LIBDIR)/include/grbinc
LPILIBOBJ	=	lpi/lpi_grb.o blockmemshell/memory.o scip/rbtree.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
SOFTLINKS	+=	$(LIBDIR)/include/grbinc
SOFTLINKS	+=	$(LIBDIR)/shared/libgurobi.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
LPIINSTMSG	=	"  -> \"grbinc\" is the path to the Gurobi \"include\" directory, e.g., \"<Gurobi-path>/include\".\n"
LPIINSTMSG	+=	" -> \"libgurobi.*\" is the path to the Gurobi library, e.g., \"<Gurobi-path>/lib/libgurobi.so\""
endif

# glop only supports shared libraries
LPSOPTIONS	+=	glop
ifeq ($(LPS),glop)
LPILIBOBJ	=	lpi/lpi_glop.o scip/bitencode.o scip/rbtree.o scip/message.o
LPILIBSRC  	=	$(SRCDIR)/lpi/lpi_glop.cpp $(SRCDIR)/scip/bitencode.c $(SRCDIR)/scip/rbtree.c $(SRCDIR)/scip/message.c
SOFTLINKS	+=	$(LIBDIR)/shared/ortools
LPIINSTMSG	=	"  -> \"ortools\" is the path to the OR-Tools directory.\n"
endif

LPSOPTIONS	+=	none
ifeq ($(LPS),none)
LPILIBOBJ	=	lpi/lpi_none.o blockmemshell/memory.o scip/rbtree.o scip/message.o
LPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(LPILIBOBJ:.o=.c))
endif

LPILIB		=	$(LPILIBNAME).$(BASE)
LPILIBFILE	=	$(LIBDIR)/$(LIBTYPE)/lib$(LPILIB).$(LIBEXT)
LPILIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(LPILIBOBJ))
LPILIBLINK	=	$(LIBDIR)/$(LIBTYPE)/lib$(LPILIBSHORTNAME).$(BASE).$(LIBEXT)
LPILIBSHORTLINK = 	$(LIBDIR)/$(LIBTYPE)/lib$(LPILIBSHORTNAME).$(LIBEXT)
ALLSRC		+=	$(LPILIBSRC)

ifeq ($(SHARED),true)
LPILIBEXTLIBS	=	$(LIBBUILD_L)$(LIBDIR)/$(LIBTYPE) $(LPSLDFLAGS)
ifneq ($(LINKRPATH),)
LPILIBEXTLIBS	+=	$(LINKRPATH)$(realpath $(LIBDIR)/$(LIBTYPE))
endif
endif

#-----------------------------------------------------------------------------
# Parallel Interface
#-----------------------------------------------------------------------------

TPILIBSHORTNAME	=	tpi$(TPI)
TPILIBNAME	=	$(TPILIBSHORTNAME)-$(VERSION)
TPILIBOBJ	=

TPIOPTIONS	+=	none
ifeq ($(TPI),none)
TPILIBOBJ	=	tpi/tpi_none.o
FLAGS		+=	-DTPI_NONE
endif

TPIOPTIONS	+=	omp
ifeq ($(TPI),omp)
TPILIBOBJ	=	tpi/tpi_openmp.o
FLAGS		+=	-DTPI_OMP
endif

TPIOPTIONS	+=	tny
ifeq ($(TPI),tny)
TPILIBOBJ	=	tpi/tpi_tnycthrd.o \
			tinycthread/tinycthread.o
FLAGS		+=	-DTPI_TNY
endif

TPILIBSRC  	=	$(addprefix $(SRCDIR)/,$(TPILIBOBJ:.o=.c))
TPILIB		=	$(TPILIBNAME).$(BASE)
TPILIBFILE	=	$(LIBDIR)/$(LIBTYPE)/lib$(TPILIB).$(LIBEXT)
TPILIBOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(TPILIBOBJ))
TPILIBLINK	=	$(LIBDIR)/$(LIBTYPE)/lib$(TPILIBSHORTNAME).$(BASE).$(LIBEXT)
TPILIBSHORTLINK = 	$(LIBDIR)/$(LIBTYPE)/lib$(TPILIBSHORTNAME).$(LIBEXT)
ALLSRC		+=	$(TPILIBSRC)

#-----------------------------------------------------------------------------
# Symmetry Interface
#-----------------------------------------------------------------------------

SYMOPTIONS	+=	none
ifeq ($(SYM),none)
SYMOBJ		=	symmetry/compute_symmetry_none.o
SYMOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(SYMOBJ))
SYMSRC  	=	$(addprefix $(SRCDIR)/,$(SYMOBJ:.o=.cpp))
ALLSRC		+=	$(SYMSRC)
endif

SYMOPTIONS	+=	bliss
ifeq ($(SYM),bliss)
SYMOBJ		=	symmetry/compute_symmetry_bliss.o
SYMOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(SYMOBJ))
SYMSRC		=       $(addprefix $(SRCDIR)/,$(SYMOBJ:.o=.cpp))
ifeq ($(BLISSEXTERNAL),false)
FLAGS		+=	-I$(SRCDIR)/bliss/src -I$(SRCDIR)/bliss/include
BLISSOBJ	=	bliss/src/abstractgraph.o
BLISSOBJ	+=	bliss/src/bliss.o
BLISSOBJ	+=	bliss/src/bliss_C.o
BLISSOBJ	+=	bliss/src/defs.o
BLISSOBJ	+=	bliss/src/digraph.o
BLISSOBJ	+=	bliss/src/graph.o
BLISSOBJ	+=	bliss/src/orbit.o
BLISSOBJ	+=	bliss/src/partition.o
BLISSOBJ	+=	bliss/src/uintseqhash.o
BLISSOBJ	+=	bliss/src/utils.o
SYMOBJFILES	+=	$(addprefix $(LIBOBJDIR)/,$(BLISSOBJ))
SYMSRC  	+=	$(addprefix $(SRCDIR)/,$(BLISSOBJ:.o=.cc))
else
FLAGS		+=	-I$(LIBDIR)/include/
endif
ALLSRC		+=	$(SYMSRC)
ifeq ($(BLISSEXTERNAL),true)
SOFTLINKS	+=	$(LIBDIR)/include/bliss
ifeq ($(SHARED),true)
SOFTLINKS	+=	$(LIBDIR)/shared/libbliss.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
else
SOFTLINKS	+=	$(LIBDIR)/static/libbliss.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
endif
endif
LPIINSTMSG	+=	"\n  -> \"blissinc\" is the path to the BLISS directory, e.g., \"<BLISS-path>\".\n"
LPIINSTMSG	+=	" -> \"libbliss.*.a\" is the path to the BLISS library, e.g., \"<BLISS-path>/libbliss.a\"\n"
LPIINSTMSG	+=	" -> \"libbliss.*.so\" is the path to the BLISS library, e.g., \"<BLISS-path>/libbliss.so\""
endif

SYMOPTIONS	+=	sbliss
ifeq ($(SYM),sbliss)
SYMOBJ		=	symmetry/compute_symmetry_sassy.o
SYMOBJFILES	=	$(addprefix $(LIBOBJDIR)/,$(SYMOBJ))
SYMSRC  	=	$(addprefix $(SRCDIR)/,$(SYMOBJ:.o=.cpp))
ifeq ($(BLISSEXTERNAL),false)
FLAGS		+=	-I$(SRCDIR)/bliss/src -I$(SRCDIR)/bliss/include
BLISSOBJ	=	bliss/src/abstractgraph.o
BLISSOBJ	+=	bliss/src/bliss.o
BLISSOBJ	+=	bliss/src/bliss_C.o
BLISSOBJ	+=	bliss/src/defs.o
BLISSOBJ	+=	bliss/src/digraph.o
BLISSOBJ	+=	bliss/src/graph.o
BLISSOBJ	+=	bliss/src/orbit.o
BLISSOBJ	+=	bliss/src/partition.o
BLISSOBJ	+=	bliss/src/uintseqhash.o
BLISSOBJ	+=	bliss/src/utils.o
SYMOBJFILES	+=	$(addprefix $(LIBOBJDIR)/,$(BLISSOBJ))
SYMSRC  	+=	$(addprefix $(SRCDIR)/,$(BLISSOBJ:.o=.cc))
else
FLAGS		+=	-I$(LIBDIR)/include/
endif
ALLSRC		+=	$(SYMSRC)
ifeq ($(BLISSEXTERNAL),true)
SOFTLINKS	+=	$(LIBDIR)/include/bliss
ifeq ($(SHARED),true)
SOFTLINKS	+=	$(LIBDIR)/shared/libbliss.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
else
SOFTLINKS	+=	$(LIBDIR)/static/libbliss.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
endif
endif
LPIINSTMSG	+=	"\n  -> \"blissinc\" is the path to the BLISS directory, e.g., \"<BLISS-path>\".\n"
LPIINSTMSG	+=	" -> \"libbliss.*.a\" is the path to the BLISS library, e.g., \"<BLISS-path>/libbliss.a\"\n"
LPIINSTMSG	+=	" -> \"libbliss.*.so\" is the path to the BLISS library, e.g., \"<BLISS-path>/libbliss.so\""
endif

SYMOPTIONS      +=      nauty
ifeq ($(SYM),nauty)
SYMOBJ          =       symmetry/compute_symmetry_nauty.o
SYMOBJFILES     =       $(addprefix $(LIBOBJDIR)/,$(SYMOBJ))
SYMSRC          =       $(addprefix $(SRCDIR)/,$(SYMOBJ:.o=.c))
ifeq ($(NAUTYEXTERNAL),false)
FLAGS           +=      -I$(SRCDIR)/nauty/src -I$(SRCDIR)/nauty/include
NAUTYOBJ        =       nauty/nauty.o
NAUTYOBJ        +=      nauty/nautil.o
NAUTYOBJ        +=      nauty/nausparse.o
NAUTYOBJ        +=      nauty/naugraph.o
NAUTYOBJ        +=      nauty/schreier.o
NAUTYOBJ        +=      nauty/naurng.o
SYMOBJFILES     +=      $(addprefix $(LIBOBJDIR)/,$(NAUTYOBJ))
SYMSRC          +=      $(addprefix $(SRCDIR)/,$(NAUTYOBJ:.o=.c))
else
FLAGS           +=      -I$(LIBDIR)/include/
endif
ALLSRC          +=      $(SYMSRC)
ifeq ($(NAUTYEXTERNAL),true)
SOFTLINKS       +=      $(LIBDIR)/include/nauty
SOFTLINKS       +=      $(LIBDIR)/static/libnauty.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
LPIINSTMSG      +=      "\n  -> \"nautyinc\" is the path to the Nauty directory, e.g., \"<Nauty-path>\".\n"
LPIINSTMSG      +=      " -> \"libnauty.*.a\" is the path to the Nauty library, e.g., \"<Nauty-path>/nauty.a\"\n"
endif
endif

#-----------------------------------------------------------------------------
# PaPILO Library
#-----------------------------------------------------------------------------

ifeq ($(PAPILO),true)
FLAGS        +=    -DSCIP_WITH_PAPILO -DPAPILO_NO_CMAKE_CONFIG -I$(LIBDIR)/include/tbb/include -I$(LIBDIR)/include/papilo/external -I$(LIBDIR)/include/papilo/src
SOFTLINKS    +=    $(LIBDIR)/include/papilo
LPIINSTMSG    +=    "\n  -> \"papilo\" is the path to the PaPILO directory\n"
SOFTLINKS    +=    $(LIBDIR)/include/boost
LPIINSTMSG    +=    "\n  -> \"boost\" is the path to the boost include folder\n"
SOFTLINKS    +=    $(LIBDIR)/include/tbb
LPIINSTMSG    +=    "\n  -> \"tbb\" is the path to the tbb include folder\n"
SOFTLINKS    +=    $(LIBDIR)/shared/libtbb.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
LPIINSTMSG    +=    " -> \"libtbb.*.so\" is the path to the oneTBB library, e.g., \"<oneTBB-path>/libtbb.so\""
endif

#-----------------------------------------------------------------------------
# External Libraries
#-----------------------------------------------------------------------------

ZLIBDEP		:=	$(SRCDIR)/depend.zlib
ZLIBSRC		:=	$(shell cat $(ZLIBDEP))

GMPDEP		:=	$(SRCDIR)/depend.gmp
GMPSRC		:=	$(shell cat $(GMPDEP))

READLINEDEP	:=	$(SRCDIR)/depend.readline
READLINESRC	:=	$(shell cat $(READLINEDEP))

ZIMPLDEP	:=	$(SRCDIR)/depend.zimpl
ZIMPLSRC	:=	$(shell cat $(ZIMPLDEP))

AMPLDEP	:=	$(SRCDIR)/depend.ampl
AMPLSRC	:=	$(shell cat $(AMPLDEP))

THREADSAFEDEP	:=	$(SRCDIR)/depend.threadsafe
THREADSAFESRC	:=	$(shell cat $(THREADSAFEDEP))

ifeq ($(ZIMPL),true)
ifeq ($(GMP),false)
$(error ZIMPL requires the GMP to be linked. Use either ZIMPL=false or GMP=true.)
endif
FLAGS		+=	-DSCIP_WITH_ZIMPL -I$(LIBDIR)/include/zimplinc $(ZIMPL_FLAGS)
DIRECTORIES	+=	$(LIBDIR)/include/zimplinc
SOFTLINKS	+=	$(LIBDIR)/include/zimplinc/zimpl
SOFTLINKS	+=	$(LIBDIR)/$(LIBTYPE)/libzimpl.$(OSTYPE).$(ARCH).$(COMP).$(ZIMPLOPT).$(STATICLIBEXT)
LPIINSTMSG	+=	"\n  -> \"zimplinc\" is a directory containing the path to the ZIMPL \"src\" directory, e.g., \"<ZIMPL-path>/src\".\n"
LPIINSTMSG	+=	" -> \"libzimpl.*\" is the path to the ZIMPL library, e.g., \"<ZIMPL-path>/lib/libzimpl.$(OSTYPE).$(ARCH).$(COMP).$(ZIMPLOPT).$(STATICLIBEXT)\""
endif

ifeq ($(GMP),true)
ifeq ($(COMP),msvc)
SOFTLINKS	+=	$(LIBDIR)/mpir.$(ARCH)
SOFTLINKS	+=	$(LIBDIR)/$(LIBTYPE)/libmpir.$(ARCH).$(OPT).lib
SOFTLINKS	+=	$(LIBDIR)/$(LIBTYPE)/libpcre.$(ARCH).$(OPT).lib
LPIINSTMSG	+=	"\n  -> \"mpir.$(ARCH)\" is a directory containing the mpir installation, i.e., \"mpir.$(ARCH)/gmp.h\" should exist.\n"
LPIINSTMSG	+=	" -> \"libmpir.*\" is the path to the MPIR library\n"
LPIINSTMSG	+=	" -> \"libpcre.*\" is the path to the PCRE library"
endif
endif

ifeq ($(IPOPT),true)
SOFTLINKS	+=	$(LIBDIR)/$(LIBTYPE)/ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)
LPIINSTMSG	+=	"\n  -> \"ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)\" is a directory containing the ipopt installation, i.e., \"ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)/include/coin/IpIpoptApplication.hpp\", \"ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)/lib/libipopt*\", ... should exist.\n"
endif

ifeq ($(FILTERSQP),true)
FLAGS		+=	-DFNAME_$(FORTRAN_NAMING_CONVENTION)
SOFTLINKS	+=	$(LIBDIR)/$(LIBTYPE)/libfiltersqp.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
#SOFTLINKS	+=	$(LIBDIR)/$(LIBTYPE)/libfiltersqp.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
SOFTLINKS	+=	$(LIBDIR)/$(LIBTYPE)/libbqpd.$(OSTYPE).$(ARCH).$(COMP).$(STATICLIBEXT)
#SOFTLINKS	+=	$(LIBDIR)/$(LIBTYPE)/libbqpd.$(OSTYPE).$(ARCH).$(COMP).$(SHAREDLIBEXT)
LPIINSTMSG	+=	" -> \"libfiltersqp.$(OSTYPE).$(ARCH).$(COMP).*\" is the path to the filterSQP library.\n"
LPIINSTMSG	+=	" -> \"libbqpd.$(OSTYPE).$(ARCH).$(COMP).*\" is the path to the BQPD library.\n"
endif

# WORHP provides only shared libraries
ifeq ($(WORHP),true)
SOFTLINKS	+=	$(LIBDIR)/$(LIBTYPE)/worhp.$(OSTYPE).$(ARCH).$(COMP).$(WORHPOPT)
LPIINSTMSG	+=	"\n  -> \"worhp.$(OSTYPE).$(ARCH).$(COMP).$(WORHPOPT)\" is a directory containing the WORHP installation, i.e., \"worhp.$(OSTYPE).$(ARCH).$(COMP).$(WORHPOPT)/include/worhp/worhp.h\" should exist.\n"
endif

ifeq ($(AMPL),true)
FLAGS		+=	-DSCIP_WITH_AMPL -I$(SRCDIR)/amplmp/include
LINKER		=	CPP
endif

ifeq ($(SHARED),true)
SCIPLIBEXTLIBS	=	$(LIBBUILD_L)$(LIBDIR)/$(LIBTYPE) $(IPOPTLIBS) $(FILTERSQPLIBS)
ifeq ($(IPOPT),true)
ifneq ($(LINKRPATH),)
SCIPLIBEXTLIBS	+=	 $(LINKRPATH)$(realpath $(LIBDIR)/$(LIBTYPE)/ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)/lib)
endif
endif
ifeq ($(WORHP),true)
ifneq ($(LINKRPATH),)
SCIPLIBEXTLIBS	+=	$(LINKRPATH)$(realpath $(LIBDIR)/$(LIBTYPE)/worhp.$(OSTYPE).$(ARCH).$(COMP).$(WORHPOPT)/lib)
endif
endif
ifeq ($(ZLIB),true)
SCIPLIBEXTLIBS	+=	$(ZLIB_LDFLAGS)
endif
ifeq ($(GMP),true)
SCIPLIBEXTLIBS	+=	$(GMP_LDFLAGS)
endif
ifeq ($(READLINE_LDFLAGS),true)
SCIPLIBEXTLIBS	+=	$(READLINE_LDFLAGS)
endif
SCIPLIBEXTLIBS	+=	$(ZIMPLLIB)
ifneq ($(LINKRPATH),)
SCIPLIBEXTLIBS	+=	$(LINKRPATH)$(realpath $(LIBDIR)/$(LIBTYPE))
endif
endif


#-----------------------------------------------------------------------------
# SCIP Library
#-----------------------------------------------------------------------------

SCIPLIBBASESHORTNAME=	scipbase
SCIPLIBBASENAME	=	$(SCIPLIBBASESHORTNAME)-$(VERSION)
SCIPPLUGINLIBOBJ=	scip/benders_default.o \
			scip/benderscut_feas.o \
			scip/benderscut_feasalt.o \
			scip/benderscut_int.o \
			scip/benderscut_nogood.o \
			scip/benderscut_opt.o \
			scip/branch_allfullstrong.o \
			scip/branch_cloud.o \
			scip/branch_distribution.o \
			scip/branch_fullstrong.o \
			scip/branch_inference.o \
			scip/branch_leastinf.o \
			scip/branch_lookahead.o \
			scip/branch_mostinf.o \
			scip/branch_multaggr.o \
			scip/branch_nodereopt.o \
			scip/branch_pscost.o \
			scip/branch_random.o \
			scip/branch_relpscost.o \
			scip/branch_vanillafullstrong.o \
			scip/compr_largestrepr.o \
			scip/compr_weakcompr.o \
			scip/concsolver_scip.o \
			scip/cons_abspower.o \
			scip/cons_and.o \
			scip/cons_benders.o \
			scip/cons_benderslp.o \
			scip/cons_bounddisjunction.o \
			scip/cons_cardinality.o \
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
			scip/cons_orbisack.o \
			scip/cons_orbitope.o \
			scip/cons_pseudoboolean.o \
			scip/cons_quadratic.o \
			scip/cons_setppc.o \
			scip/cons_soc.o \
			scip/cons_sos1.o \
			scip/cons_sos2.o \
			scip/cons_superindicator.o \
			scip/cons_symresack.o \
			scip/cons_varbound.o \
			scip/cons_xor.o \
			scip/cons_components.o \
			scip/cutsel_hybrid.o \
			scip/dialog_default.o \
			scip/event_softtimelimit.o \
			scip/disp_default.o \
			scip/event_solvingphase.o \
			scip/prop_sync.o \
			scip/event_globalbnd.o \
			scip/event_estim.o \
			scip/expr_abs.o \
			scip/expr_entropy.o \
			scip/expr_erf.o \
			scip/expr_exp.o \
			scip/expr_log.o \
			scip/expr_pow.o \
			scip/expr_product.o \
			scip/expr_sum.o \
			scip/expr_trig.o \
			scip/expr_value.o \
			scip/expr_var.o \
			scip/expr_varidx.o \
			scip/heur_sync.o \
			scip/heur_actconsdiving.o \
			scip/heur_adaptivediving.o \
			scip/heur_bound.o \
			scip/heur_clique.o \
			scip/heur_coefdiving.o \
			scip/heur_completesol.o \
			scip/heur_conflictdiving.o \
			scip/heur_crossover.o \
			scip/heur_dins.o \
			scip/heur_distributiondiving.o \
			scip/heur_dps.o \
			scip/heur_dualval.o \
			scip/heur_farkasdiving.o \
			scip/heur_feaspump.o \
			scip/heur_fixandinfer.o \
			scip/heur_fracdiving.o \
			scip/heur_gins.o \
			scip/heur_guideddiving.o \
			scip/heur_indicator.o \
			scip/heur_intdiving.o \
			scip/heur_intshifting.o \
			scip/heur_linesearchdiving.o \
			scip/heur_localbranching.o \
			scip/heur_lpface.o \
			scip/heur_alns.o \
			scip/heur_locks.o \
			scip/heur_multistart.o \
			scip/heur_mutation.o \
			scip/heur_mpec.o \
			scip/heur_nlpdiving.o \
			scip/heur_objpscostdiving.o \
			scip/heur_octane.o \
			scip/heur_ofins.o \
			scip/heur_oneopt.o \
			scip/heur_padm.o \
			scip/heur_proximity.o \
			scip/heur_pscostdiving.o \
			scip/heur_reoptsols.o \
			scip/heur_repair.o \
			scip/heur_randrounding.o \
			scip/heur_rens.o \
			scip/heur_rins.o \
			scip/heur_rootsoldiving.o \
			scip/heur_rounding.o \
			scip/heur_shiftandpropagate.o \
			scip/heur_shifting.o \
			scip/heur_simplerounding.o \
			scip/heur_subnlp.o \
			scip/heur_trivial.o \
			scip/heur_trivialnegation.o \
			scip/heur_trustregion.o \
			scip/heur_trysol.o \
			scip/heur_twoopt.o \
			scip/heur_undercover.o \
			scip/heur_vbounds.o \
			scip/heur_veclendiving.o \
			scip/heur_zeroobj.o \
			scip/heur_zirounding.o \
			scip/message_default.o \
			scip/nlhdlr_bilinear.o \
			scip/nlhdlr_convex.o \
			scip/nlhdlr_default.o \
			scip/nlhdlr_perspective.o \
			scip/nlhdlr_quadratic.o \
			scip/nlhdlr_quotient.o \
			scip/nlhdlr_soc.o \
			scip/nlpi_all.o \
			scip/nodesel_bfs.o \
			scip/nodesel_breadthfirst.o \
			scip/nodesel_dfs.o \
			scip/nodesel_estimate.o \
			scip/nodesel_hybridestim.o \
			scip/nodesel_restartdfs.o \
			scip/nodesel_uct.o \
			scip/presol_boundshift.o \
			scip/presol_convertinttobin.o \
			scip/presol_domcol.o\
			scip/presol_dualagg.o\
			scip/presol_dualcomp.o\
			scip/presol_dualinfer.o\
			scip/presol_gateextraction.o \
			scip/presol_implics.o \
			scip/presol_inttobinary.o \
			scip/presol_qpkktref.o \
			scip/presol_redvub.o \
			scip/presol_trivial.o \
			scip/presol_tworowbnd.o \
			scip/presol_sparsify.o \
			scip/presol_dualsparsify.o \
			scip/presol_stuffing.o \
			scip/prop_dualfix.o \
			scip/prop_genvbounds.o \
			scip/prop_nlobbt.o \
			scip/prop_obbt.o \
			scip/prop_probing.o \
			scip/prop_pseudoobj.o \
			scip/prop_redcost.o \
			scip/prop_rootredcost.o \
			scip/prop_symmetry.o \
			scip/prop_vbounds.o \
			scip/reader_bnd.o \
			scip/reader_ccg.o \
			scip/reader_cip.o \
			scip/reader_cnf.o \
			scip/reader_cor.o \
			scip/reader_dec.o \
			scip/reader_diff.o \
			scip/reader_fix.o \
			scip/reader_fzn.o \
			scip/reader_gms.o \
			scip/reader_lp.o \
			scip/reader_mps.o \
			scip/reader_mst.o \
			scip/reader_opb.o \
			scip/reader_osil.o \
			scip/reader_pip.o \
			scip/reader_pbm.o \
			scip/reader_ppm.o \
			scip/reader_rlp.o \
			scip/reader_smps.o \
			scip/reader_sol.o \
			scip/reader_sto.o \
			scip/reader_tim.o \
			scip/reader_wbo.o \
			scip/reader_zpl.o \
			scip/sepa_aggregation.o \
			scip/sepa_cgmip.o \
			scip/sepa_clique.o \
			scip/sepa_closecuts.o \
			scip/sepa_convexproj.o \
			scip/sepa_disjunctive.o \
			scip/sepa_eccuts.o \
			scip/sepa_gauge.o \
			scip/sepa_gomory.o \
			scip/sepa_impliedbounds.o \
			scip/sepa_interminor.o \
			scip/sepa_intobj.o \
			scip/sepa_mcf.o \
			scip/sepa_minor.o \
			scip/sepa_mixing.o \
			scip/sepa_oddcycle.o \
			scip/sepa_rapidlearning.o \
			scip/sepa_rlt.o \
			scip/sepa_zerohalf.o \
			scip/table_default.o


SCIPPLUGINLIBCPPOBJ =	scip/presol_milp.o

ifeq ($(EXPRINT),none)
SCIPPLUGINLIBOBJ 	+=	scip/exprinterpret_none.o
endif
ifeq ($(EXPRINT),cppad)
SCIPPLUGINLIBCPPOBJ 	+= 	scip/exprinterpret_cppad.o
endif

ifeq ($(IPOPT),true)
SCIPPLUGINLIBCPPOBJ	+= 	scip/nlpi_ipopt.o
else
SCIPPLUGINLIBOBJ	+= 	scip/nlpi_ipopt_dummy.o
endif

ifeq ($(FILTERSQP),true)
SCIPPLUGINLIBOBJ	+= scip/nlpi_filtersqp.o
else
SCIPPLUGINLIBOBJ	+= scip/nlpi_filtersqp_dummy.o
endif

ifeq ($(WORHP),true)
SCIPPLUGINLIBOBJ	+= 	scip/nlpi_worhp.o
else
SCIPPLUGINLIBOBJ	+= 	scip/nlpi_worhp_dummy.o
endif

ifeq ($(AMPL),true)
SCIPPLUGINLIBCPPOBJ += scip/reader_nl.o
SCIPPLUGINLIBCPPOBJ += amplmp/src/format.o amplmp/src/expr-info.o amplmp/src/nl-reader.o amplmp/src/os.o amplmp/src/posix.o
endif

SCIPLIBOBJ	=	scip/boundstore.o \
			scip/branch.o \
			scip/bandit.o \
			scip/bandit_epsgreedy.o \
			scip/bandit_exp3.o \
			scip/bandit_ucb.o \
			scip/benders.o \
			scip/benderscut.o \
			scip/bendersdefcuts.o \
			scip/clock.o \
			scip/concsolver.o \
			scip/concurrent.o \
			scip/conflict.o \
			scip/conflictstore.o \
			scip/cons.o \
			scip/cutpool.o \
			scip/cuts.o \
			scip/cutsel.o \
			scip/debug.o \
			scip/dcmp.o \
			scip/dialog.o \
			scip/disp.o \
			scip/event.o \
			scip/expr.o \
			scip/exprcurv.o \
			scip/expriter.o \
			scip/fileio.o \
			scip/heur.o \
			scip/heuristics.o \
			scip/compr.o \
			scip/history.o \
			scip/implics.o \
			scip/interrupt.o \
			scip/intervalarith.o \
			scip/lp.o \
			scip/matrix.o \
			scip/mem.o \
			scip/misc.o \
			scip/misc_linear.o \
			scip/misc_rowprep.o \
			scip/nlhdlr.o \
			scip/nlp.o \
			scip/nlpi.o \
			scip/nlpioracle.o \
			scip/nodesel.o \
			scip/paramset.o \
			scip/presol.o \
			scip/presolve.o \
			scip/pricestore.o \
			scip/pricer.o \
			scip/primal.o \
			scip/prob.o \
			scip/prop.o \
			scip/reader.o \
			scip/relax.o \
			scip/reopt.o \
			scip/retcode.o \
			scip/scip_benders.o \
			scip/scip_branch.o \
			scip/scip_compr.o \
			scip/scip_concurrent.o \
			scip/scip_conflict.o \
			scip/scip_cons.o \
			scip/scip_copy.o \
			scip/scip_cut.o \
			scip/scip_cutsel.o \
			scip/scip_datastructures.o\
			scip/scip_debug.o \
			scip/scip_dcmp.o \
			scip/scip_dialog.o \
			scip/scip_disp.o \
			scip/scip_event.o \
			scip/scip_expr.o \
			scip/scip_general.o \
			scip/scip_heur.o \
			scip/scip_lp.o \
			scip/scip_mem.o \
			scip/scip_message.o \
			scip/scip_nlp.o \
			scip/scip_nlpi.o \
			scip/scip_nodesel.o \
			scip/scip_numerics.o \
			scip/scip_param.o \
			scip/scip_presol.o \
			scip/scip_pricer.o \
			scip/scip_prob.o \
			scip/scip_probing.o \
			scip/scip_prop.o \
			scip/scip_randnumgen.o \
			scip/scip_reader.o \
			scip/scip_relax.o \
			scip/scip_reopt.o \
			scip/scip_sepa.o \
			scip/scip_sol.o \
			scip/scip_solve.o \
			scip/scip_solvingstats.o \
			scip/scip_table.o \
			scip/scip_timing.o \
			scip/scip_tree.o \
			scip/scip_validation.o \
			scip/scip_var.o \
			scip/scip_bandit.o \
			scip/scipbuildflags.o \
			scip/scipcoreplugins.o \
			scip/scipdefplugins.o \
			scip/scipgithash.o \
			scip/scipshell.o \
			scip/sepa.o \
			scip/sepastore.o \
			scip/set.o \
			scip/sol.o \
			scip/solve.o \
			scip/stat.o \
			scip/symmetry.o \
			scip/syncstore.o \
			scip/table.o \
			scip/tree.o \
			scip/treemodel.o \
			scip/var.o \
			scip/visual.o \
			tclique/tclique_branch.o \
			tclique/tclique_coloring.o \
			tclique/tclique_graph.o \
			dijkstra/dijkstra.o \
			xml/xmlparse.o

SCIPLIBBASE	=	$(SCIPLIBBASENAME).$(BASE)
SCIPLIBBASEFILE	=	$(LIBDIR)/$(LIBTYPE)/lib$(SCIPLIBBASE).$(LIBEXT)
SCIPLIBBASEOBJFILES =	$(addprefix $(LIBOBJDIR)/,$(SCIPPLUGINLIBOBJ))
SCIPLIBBASEOBJFILES +=	$(addprefix $(LIBOBJDIR)/,$(SCIPPLUGINLIBCPPOBJ))
SCIPLIBBASEOBJFILES +=	$(addprefix $(LIBOBJDIR)/,$(SCIPLIBOBJ))
SCIPLIBBASESRC	=	$(addprefix $(SRCDIR)/,$(SCIPPLUGINLIBOBJ:.o=.c))
SCIPLIBBASESRC	+=	$(addprefix $(SRCDIR)/,$(SCIPPLUGINLIBCPPOBJ:.o=.cpp))
SCIPLIBBASESRC	+=	$(addprefix $(SRCDIR)/,$(SCIPLIBOBJ:.o=.c))
SCIPLIBBASELINK	 =	$(LIBDIR)/$(LIBTYPE)/lib$(SCIPLIBBASESHORTNAME).$(BASE).$(LIBEXT)
SCIPLIBBASESHORTLINK = 	$(LIBDIR)/$(LIBTYPE)/lib$(SCIPLIBBASESHORTNAME).$(LIBEXT)

# define library that contains everything
SCIPLIBSHORTNAME = scip
SCIPLIBNAME =		$(SCIPLIBSHORTNAME)-$(VERSION)
SCIPLIB	=		$(SCIPLIBNAME).$(BASE).$(LPS)
SCIPLIBFILE =		$(LIBDIR)/$(LIBTYPE)/lib$(SCIPLIB).$(LIBEXT)
SCIPLIBLINK =		$(LIBDIR)/$(LIBTYPE)/lib$(SCIPLIBSHORTNAME).$(BASE).$(LPS).$(LIBEXT)
SCIPLIBSHORTLINK =	$(LIBDIR)/$(LIBTYPE)/lib$(SCIPLIBSHORTNAME).$(LIBEXT)

# for backward compatibility
SCIPLIBSOLVERLINK =	$(LIBDIR)/$(LIBTYPE)/libscipsolver.$(BASE).$(LPS).$(LIBEXT)
SCIPLIBSOLVERSHORTLINK = $(LIBDIR)/$(LIBTYPE)/libscipsolver.$(LIBEXT)


ALLSRC		+=	$(SCIPLIBBASESRC)

SCIPGITHASHFILE	= 	$(SRCDIR)/scip/githash.c
SCIPBUILDFLAGSFILE = 	$(SRCDIR)/scip/buildflags.c

#-----------------------------------------------------------------------------
# Objective SCIP Library
#-----------------------------------------------------------------------------

OBJSCIPLIBSHORTNAME=	objscip
OBJSCIPLIBNAME	=	$(OBJSCIPLIBSHORTNAME)-$(VERSION)
OBJSCIPLIBOBJ	=	objscip/objbenders.o \
			objscip/objbenderscut.o \
			objscip/objbranchrule.o \
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
			objscip/objtable.o \
			objscip/objvardata.o

OBJSCIPLIB	=	$(OBJSCIPLIBNAME).$(BASE)
OBJSCIPLIBFILE	=	$(LIBDIR)/$(LIBTYPE)/lib$(OBJSCIPLIB).$(LIBEXT)
OBJSCIPLIBOBJFILES=	$(addprefix $(LIBOBJDIR)/,$(OBJSCIPLIBOBJ))
OBJSCIPLIBSRC	=	$(addprefix $(SRCDIR)/,$(OBJSCIPLIBOBJ:.o=.cpp))
OBJSCIPINCSRC	=	$(addprefix $(SRCDIR)/,$(OBJSCIPLIBOBJ:.o=.h))
OBJSCIPLIBLINK	=	$(LIBDIR)/$(LIBTYPE)/lib$(OBJSCIPLIBSHORTNAME).$(BASE).$(LIBEXT)
OBJSCIPLIBSHORTLINK=	$(LIBDIR)/$(LIBTYPE)/lib$(OBJSCIPLIBSHORTNAME).$(LIBEXT)
ALLSRC		+=	$(OBJSCIPLIBSRC)


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINSHORTNAME	=	scip
MAINNAME	=	$(MAINSHORTNAME)-$(VERSION)

MAINOBJ		=	main.o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.c))

MAINFILE	=	$(BINDIR)/$(MAINNAME).$(BASE).$(LPS).$(TPI)$(EXEEXTENSION)
MAINOBJFILES	=	$(addprefix $(BINOBJDIR)/,$(MAINOBJ))
MAINLINK	=	$(BINDIR)/$(MAINSHORTNAME).$(BASE).$(LPS).$(TPI)$(EXEEXTENSION)
MAINSHORTLINK	=	$(BINDIR)/$(MAINSHORTNAME)$(EXEEXTENSION)
ALLSRC		+=	$(MAINSRC)

ifeq ($(SHARED),true)
WINLIBFILENAME	=	lib$(MAINNAME).$(BASE).$(LPS).dll
else
WINLIBFILENAME	=	lib$(MAINNAME).$(BASE).$(LPS).lib
endif

LINKSMARKERFILE	=	$(LIBDIR)/$(LIBTYPE)/linkscreated.$(LPS)-$(LPSOPT).$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX).$(ZIMPL)-$(ZIMPLOPT).$(IPOPT)-$(IPOPTOPT).$(FILTERSQP).$(SYM).$(PAPILO)
LASTSETTINGS	=	$(OBJDIR)/make.lastsettings

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(SCIPLIBBASEFILE) $(OBJSCIPLIBFILE) $(LPILIBFILE) $(TPILIBFILE) \
		$(LPILIBLINK) $(LPILIBSHORTLINK) $(TPILIBLINK) $(TPILIBSHORTLINK) $(SCIPLIBBASELINK) $(SCIPLIBBASESHORTLINK) \
		$(OBJSCIPLIBLINK) $(OBJSCIPLIBSHORTLINK) \
		$(MAINLINK) $(MAINSHORTLINK) \
		$(LPILIBOBJFILES) $(TPILIBOBJFILES) $(SCIPLIBBASEOBJFILES) $(OBJSCIPLIBOBJFILES) $(MAINOBJFILES) $(SYMOBJFILES) \
		$(SCIPLIBFILE) $(SCIPLIBLINK) $(SCIPLIBSHORTLINK)
MAKE		+= -s
endif

.PHONY: all
all:		libs
		@$(MAKE) $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)

.PHONY: libs
libs:		libscipbase libobjscip liblpi libtpi
ifeq ($(SHARED),true)
		@$(MAKE) libscip
endif

.PHONY: preprocess
preprocess:     checkdefines
		@$(SHELL) -ec 'if test ! -e $(LINKSMARKERFILE) ; \
			then \
				echo "-> generating necessary links" ; \
				$(MAKE) -j1 $(LINKSMARKERFILE) ; \
			fi'
		@$(MAKE) touchexternal

.PHONY: lint
lint:		$(SCIPLIBBASESRC) $(OBJSCIPLIBSRC) $(LPILIBSRC) $(TPILIBSRC) $(MAINSRC) $(SYMSRC) githash
		-rm -f lint.out

		@$(SHELL) -ec 'if test -e lint/co-gcc.mak ; \
			then \
				echo "-> generating gcc-include-path lint-file" ; \
				cd lint; $(MAKE) -f co-gcc.mak ; \
			else \
				echo "-> lint Makefile not found"; \
			fi'
ifeq ($(FILES),)
		$(SHELL) -ec 'for i in $^; \
			do \
				echo $$i; \
				$(LINT) lint/main-gcc.lnt +os\(lint.out\) -u -zero \
				$(USRFLAGS) $(FLAGS) -I/usr/include -UNDEBUG -USCIP_WITH_READLINE -USCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'
else
		$(SHELL) -ec  'for i in $(FILES); \
			do \
				echo $$i; \
				$(LINT) lint/main-gcc.lnt +os\(lint.out\) -u -zero \
				$(USRFLAGS) $(FLAGS) -I/usr/include -UNDEBUG -USCIP_WITH_READLINE -USCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'
endif

.PHONY: pclint
pclint:		$(SCIPLIBBASESRC) $(OBJSCIPLIBSRC) $(LPILIBSRC) $(TPILIBSRC) $(MAINSRC) $(SYMSRC)
		-rm -f pclint.out

		@$(SHELL) -ec 'if ! test -e pclint/co-gcc.h ; \
			then \
				echo "-> running pclint configuration"; \
				CCPATH=`which $(CC)`; \
				echo "-> path to compiler: "$${CCPATH}; \
				cd pclint; \
				echo "-> running $(PCLINTCONFIG)"; \
				python $(PCLINTCONFIG) --compiler=$(CC) --compiler-bin=$${CCPATH} --config-output-lnt-file=co-gcc.lnt --config-output-header-file=co-gcc.h --generate-compiler-config ; \
			fi'
ifeq ($(FILES),)
		@$(SHELL) -ec 'echo "-> running pclint on $^..."; \
			$(PCLINT) pclint/main-gcc.lnt +os\(pclint.out\) -u -zero -max_threads=$(MAXJOBS) \
			$(USRFLAGS) $(FLAGS) -Ipclint -uNDEBUG -uSCIP_WITH_READLINE -uSCIP_ROUNDING_FE -D_BSD_SOURCE $^;'
else
		@$(SHELL) -ec  'echo "-> running pclint on files $(FILES) ..."; \
			$(PCLINT) pclint/main-gcc.lnt +os\(pclint.out\) -u -zero -max_threads=$(MAXJOBS) \
			$(USRFLAGS) $(FLAGS) -Ipclint -uNDEBUG -uSCIP_WITH_READLINE -uSCIP_ROUNDING_FE -D_BSD_SOURCE $(FILES);'
endif

.PHONY: splint
splint:		$(SCIPLIBBASESRC) $(OBJSCIPLIBSRC) $(LPILIBSRC) $(TPILIBSRC) $(MAINSRC) $(SYMSRC)
		-rm -f splint.out
ifeq ($(FILES),)
		$(SHELL) -c '$(SPLINT) -I$(SRCDIR) -I/usr/include/linux $(FLAGS) $(SPLINTFLAGS) $(filter %.c %.h,$^) >> splint.out;'
else
		$(SHELL) -c '$(SPLINT) -I$(SRCDIR) -I/usr/include/linux $(FLAGS) $(SPLINTFLAGS) $(filter %.c %.h,$(FILES)) >> splint.out;'
endif

.PHONY: doc
doc:
		cd doc; $(SHELL) builddoc.sh;

.PHONY: docpreview
docpreview:
# generates preview for a list of files
ifneq ($(FILES),)
		echo "generating doxygen preview for $(FILES)"
		cd doc; ( cat $(MAINSHORTNAME).dxy && echo 'FILE_PATTERNS = $(FILES)' ) | $(DOXY) -
else
		echo "please specify file(s) for which preview should be created"
endif
.PHONY: check
check:		test

.PHONY: test
test:
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(EXECUTABLE) $(SETTINGS) $(BINID) $(OUTPUTDIR) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) \
		$(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(DEBUGTOOL) $(CLIENTTMPDIR) $(REOPT) $(OPTCOMMAND) $(SETCUTOFF) $(MAXJOBS) $(VISUALIZE) $(PERMUTE) \
                $(SEEDS) $(GLBSEEDSHIFT) $(STARTPERM) $(PYTHON);

.PHONY: testcount
testcount:
		cd check; \
		$(SHELL) ./check_count.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(FEASTOL) \
		$(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(OUTPUTDIR);

.PHONY: tags
tags:
		rm -f TAGS; ctags -e -R -h ".c.cpp.h" --exclude=".*" src/; sed -i 's!\#undef .*!!g' TAGS

# include target to detect the current git hash
-include make/local/make.detectgithash

# this empty target is needed for the SCIP release versions
githash::      # do not remove the double-colon

# include local targets
-include make/local/make.targets

# include install/uninstall targets
-include make/make.install

# the testgams target need to come after make/local/make.targets has been included (if any), because the latter may assign a value to CLIENTTMPDIR
.PHONY: testgams
testgams:
		cd check; \
		$(SHELL) ./check_gamscluster.sh $(TEST) $(GAMS) "$(GAMSSOLVER)" $(SETTINGS) $(OSTYPE).$(ARCH) $(TIME) $(NODES) $(MEM) "$(GAP)" \
		$(THREADS) $(CONTINUE) "$(CONVERTSCIP)" local dummy dummy "$(CLIENTTMPDIR)" 1 true $(SETCUTOFF);

$(LPILIBLINK):	$(LPILIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(LPILIBFILE)) $(notdir $@)

# the short link targets should be phony such that they are always updated and point to the files with last make options, even if nothing needed to be rebuilt
.PHONY: $(LPILIBSHORTLINK)
$(LPILIBSHORTLINK):	$(LPILIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(LPILIBFILE)) $(notdir $@)

$(TPILIBLINK):	$(TPILIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(TPILIBFILE)) $(notdir $@)

# the short link targets should be phony such that they are always updated and point to the files with last make options, even if nothing needed to be rebuilt
.PHONY: $(TPILIBSHORTLINK)
$(TPILIBSHORTLINK):	$(TPILIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(TPILIBFILE)) $(notdir $@)

$(SCIPLIBBASELINK): 	$(SCIPLIBBASEFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(SCIPLIBBASEFILE)) $(notdir $@)

# the short link targets should be phony such that they are always updated and point to the files with last make options, even if nothing needed to be rebuilt
.PHONY: $(SCIPLIBBASESHORTLINK)
$(SCIPLIBBASESHORTLINK):	$(SCIPLIBBASEFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(SCIPLIBBASEFILE)) $(notdir $@)

$(OBJSCIPLIBLINK):	$(OBJSCIPLIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(OBJSCIPLIBFILE)) $(notdir $@)

# the short link targets should be phony such that they are always updated and point to the files with last make options, even if nothing needed to be rebuilt
.PHONY: $(OBJSCIPLIBSHORTLINK)
$(OBJSCIPLIBSHORTLINK):	$(OBJSCIPLIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(OBJSCIPLIBFILE)) $(notdir $@)

$(SCIPLIBLINK): $(SCIPLIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(SCIPLIBFILE)) $(notdir $@)

# the short link targets should be phony such that they are always updated and point to the files with last make options, even if nothing needed to be rebuilt
.PHONY: $(SCIPLIBSHORTLINK)
$(SCIPLIBSHORTLINK): $(SCIPLIBFILE)
		@rm -f $@ $(SCIPLIBSOLVERLINK) $(SCIPLIBSOLVERSHORTLINK)
		cd $(dir $@) && $(LN_s) $(notdir $(SCIPLIBFILE)) $(notdir $@)
		# for backward compatibility:
		cd $(dir $@) && $(LN_s) $(notdir $(SCIPLIBFILE)) $(notdir $(SCIPLIBSOLVERLINK))
		cd $(dir $@) && $(LN_s) $(notdir $(SCIPLIBFILE)) $(notdir $(SCIPLIBSOLVERSHORTLINK))

# the short link targets should be phony such that they are always updated and point to the files with last make options, even if nothing needed to be rebuilt
.PHONY: $(MAINSHORTLINK)
$(MAINLINK) $(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINOBJDIR):	| $(OBJDIR)
		@-mkdir -p $(BINOBJDIR)

$(LIBOBJDIR):	| $(OBJDIR)
		@-mkdir -p $(LIBOBJDIR)

$(LIBOBJSUBDIRS):	| $(LIBOBJDIR)
		@-mkdir -p $(LIBOBJSUBDIRS)

$(LIBDIR):
		@-mkdir -p $(LIBDIR)

$(LIBDIR)/static: $(LIBDIR)
		@-mkdir -p $(LIBDIR)/static

$(LIBDIR)/shared: $(LIBDIR)
		@-mkdir -p $(LIBDIR)/shared

$(LIBDIR)/include: $(LIBDIR)
		@-mkdir -p $(LIBDIR)/include

$(LIBDIR)/src: $(LIBDIR)
		@-mkdir -p $(LIBDIR)/src

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:          cleanlibs cleanbin | $(LIBOBJSUBDIRS) $(LIBOBJDIR) $(BINOBJDIR) $(OBJDIR)
ifneq ($(LIBOBJDIR),)
		@-(cd $(LIBOBJDIR) && rm -f */*.o */*.d */*/*.o */*/*.d)
		@-rmdir $(LIBOBJSUBDIRS)
		@-rmdir $(LIBOBJDIR)
endif
ifneq ($(BINOBJDIR),)
		@-rm -f $(BINOBJDIR)/*.o $(BINOBJDIR)/*.d && rmdir $(BINOBJDIR)
endif
ifneq ($(OBJDIR),)
		@-rm -f $(LASTSETTINGS)
		@-rmdir $(OBJDIR)
endif

.PHONY: cleanlibs
cleanlibs:      | $(LIBDIR)/$(LIBTYPE)
		@echo "-> remove library $(SCIPLIBBASEFILE)"
		@-rm -f $(SCIPLIBBASEFILE) $(SCIPLIBBASELINK) $(SCIPLIBBASESHORTLINK)
		@echo "-> remove library $(OBJSCIPLIBFILE)"
		@-rm -f $(OBJSCIPLIBFILE) $(OBJSCIPLIBLINK) $(OBJSCIPLIBSHORTLINK)
		@echo "-> remove library $(LPILIBFILE)"
		@-rm -f $(LPILIBFILE) $(LPILIBLINK) $(LPILIBSHORTLINK)
		@echo "-> remove library $(TPILIBFILE)"
		@-rm -f $(TPILIBFILE) $(TPILIBLINK) $(TPILIBSHORTLINK)
		@echo "-> remove library $(SCIPLIBFILE)"
		@-rm -f $(SCIPLIBFILE) $(SCIPLIBLINK) $(SCIPLIBSHORTLINK) $(SCIPLIBSOLVERLINK) $(SCIPLIBSOLVERSHORTLINK)

.PHONY: cleanbin
cleanbin:       | $(BINDIR)
		@echo "-> remove binary $(MAINFILE)"
		@-rm -f $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)

depend:
# We explicitely add all lpi's here, since the content of depend.lpscheck should be independent of the currently selected LPI,
# but contain all LPI's that use the SCIP_WITH_LPSCHECK define.
		@echo `grep -l "SCIP_WITH_LPSCHECK" $(SCIPLIBBASESRC) $(OBJSCIPLIBSRC) $(MAINSRC) src/lpi/lpi*.{c,cpp}` >$(LPSCHECKDEP)
		@echo `grep -l "SCIP_WITH_ZLIB" $(ALLSRC)` >$(ZLIBDEP)
		@echo `grep -l "SCIP_WITH_GMP" $(ALLSRC)` >$(GMPDEP)
		@echo `grep -l "SCIP_WITH_READLINE" $(ALLSRC)` >$(READLINEDEP)
		@echo `grep -l "SCIP_WITH_ZIMPL" $(ALLSRC)` >$(ZIMPLDEP)
		@echo `grep -l "SCIP_WITH_AMPL" $(ALLSRC)` >$(AMPLDEP)
		@echo `grep -l "SCIP_THREADSAFE" $(ALLSRC) src/lpi/lpi*.{c,cpp} src/scip/nlpi*.{c,cpp}` >$(THREADSAFEDEP)

# do not attempt to include .d files if there will definitely be any (empty DFLAGS), because it slows down the build on Windows considerably
ifneq ($(DFLAGS),)
-include $(MAINOBJFILES:.o=.d)
-include $(SCIPLIBBASEOBJFILES:.o=.d)
-include $(OBJSCIPOBJFILES:.o=.d)
-include $(LPILIBOBJFILES:.o=.d)
-include $(TPILIBOBJFILES:.o=.d)
else
ifeq ($(VERBOSE),true)
$(info No compilation dependencies. If changing header files, do a make clean before building.)
endif
endif

# make binary
$(MAINFILE):	$(MAINOBJFILES) $(SCIPLIBBASEFILE) $(OBJSCIPLIBFILE) $(LPILIBFILE) $(TPILIBFILE) | $(BINDIR) $(BINOBJDIR) $(LIBOBJSUBDIRS)
		@echo "-> linking $@"
ifeq ($(LINKER),C)
		$(LINKCC) $(MAINOBJFILES) $(LINKCCSCIPALL) $(LINKCC_o)$@ \
		|| ($(MAKE) errorhints && false)
endif
ifeq ($(LINKER),CPP)
		$(LINKCXX) $(MAINOBJFILES) $(LINKCCSCIPALL) $(LINKCXX_o)$@ \
		|| ($(MAKE) errorhints && false)
endif

.PHONY: libscipbase
libscipbase:	preprocess
		@$(MAKE) $(SCIPLIBBASEFILE) $(SCIPLIBBASELINK) $(SCIPLIBBASESHORTLINK)

$(SCIPLIBBASEFILE):	$(SCIPLIBBASEOBJFILES) $(SYMOBJFILES) | $(LIBDIR)/$(LIBTYPE) $(LIBOBJSUBDIRS)
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(SCIPLIBBASEOBJFILES) $(SYMOBJFILES) $(SCIPLIBEXTLIBS)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

.PHONY: libobjscip
libobjscip:	preprocess
		@$(MAKE) $(OBJSCIPLIBFILE) $(OBJSCIPLIBLINK) $(OBJSCIPLIBSHORTLINK)

$(OBJSCIPLIBFILE):	$(OBJSCIPLIBOBJFILES) | $(LIBOBJSUBDIRS) $(LIBDIR)/$(LIBTYPE)
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(OBJSCIPLIBOBJFILES)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

.PHONY: liblpi
liblpi:		preprocess
		@$(MAKE) $(LPILIBFILE) $(LPILIBLINK) $(LPILIBSHORTLINK)

$(LPILIBFILE):	$(LPILIBOBJFILES) | $(LIBOBJSUBDIRS) $(LIBDIR)/$(LIBTYPE)
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(LPILIBOBJFILES) $(LPILIBEXTLIBS)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

.PHONY: libtpi
libtpi:		preprocess
		@$(MAKE) $(TPILIBFILE) $(TPILIBLINK) $(TPILIBSHORTLINK)

$(TPILIBFILE):	$(TPILIBOBJFILES) | $(LIBOBJSUBDIRS) $(LIBDIR)/$(LIBTYPE)
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(TPILIBOBJFILES)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

ifneq ($(LINKRPATH),)
LIBSCIPRPATHARG = $(LINKRPATH)$(SCIPREALPATH)/$(LIBDIR)/shared
endif

.PHONY: libscip
libscip:	preprocess
		@$(MAKE) $(SCIPLIBFILE) $(SCIPLIBLINK) $(SCIPLIBSHORTLINK)

$(SCIPLIBFILE): $(SCIPLIBBASEOBJFILES) $(LPILIBOBJFILES) $(TPILIBOBJFILES) $(SYMOBJFILES) $(OBJSCIPLIBOBJFILES) | $(LIBDIR)/$(LIBTYPE) $(LIBOBJSUBDIRS)
		@echo "-> generating library $@"
		-rm -f $@
ifeq ($(SHARED),false)
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(SCIPLIBBASEOBJFILES) $(LPILIBOBJFILES) $(TPILIBOBJFILES) $(SYMOBJFILES) $(OBJSCIPLIBOBJFILES)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif
else
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ $(SCIPLIBBASEOBJFILES) $(LPILIBOBJFILES) $(LPILIBEXTLIBS) $(TPILIBOBJFILES) $(SYMOBJFILES) $(OBJSCIPLIBOBJFILES) $(SCIPLIBEXTLIBS) \
		$(LPSLDFLAGS) $(LDFLAGS) $(LIBSCIPRPATHARG)
endif

$(BINOBJDIR)/%.o:	$(SRCDIR)/%.c | $(BINOBJDIR)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) $(DFLAGS) $(TPICFLAGS) $(CC_c)$< $(CC_o)$@

$(BINOBJDIR)/%.o:	$(SRCDIR)/%.cpp | $(BINOBJDIR)
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) $(DFLAGS) $(TPICFLAGS) $(CXX_c)$< $(CXX_o)$@

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.c | $(LIBOBJDIR) $(LIBOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CFLAGS) $(DFLAGS) $(TPICFLAGS) $(CC_c)$< $(CC_o)$@

# add special target for papilo to avoid sanitizers (leading to false positives in TBB)
ifeq ($(SANITIZE),true)
$(LIBOBJDIR)/scip/presol_milp.o: $(SRCDIR)/scip/presol_milp.cpp | $(LIBOBJDIR) $(LIBOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CXXFLAGS) $(DFLAGS) $(TPICFLAGS) $(CXX_c)$< $(CXX_o)$@ -fno-sanitize=all
endif

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.cpp | $(LIBOBJDIR) $(LIBOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CXXFLAGS) $(DFLAGS) $(TPICFLAGS) $(CXX_c)$< $(CXX_o)$@

$(LIBOBJDIR)/%.o:	$(SRCDIR)/%.cc | $(LIBOBJDIR) $(LIBOBJSUBDIRS)
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(LIBOFLAGS) $(CXXFLAGS) $(DFLAGS) $(TPICFLAGS) $(CXX_c)$< $(CXX_o)$@

-include $(LASTSETTINGS)

.PHONY: windowslib
windowslib: $(SCIPLIBBASEOBJFILES) $(MAINOBJFILES) $(LPILIBOBJFILES) $(OBJSCIPLIBOBJFILES) $(TPILIBOBJFILES) $(SYMOBJFILES) | $(LIBOBJSUBDIRS) $(LIBDIR)/$(LIBTYPE)
		@echo "-> generating Windows library $@"
ifeq ($(SHARED),true)
		$(LINKCC) $(LIBBUILDFLAGS) $(LINKCC_L)$(LIBDIR)/$(LIBTYPE) $(LIBBUILD_o)$(LIBDIR)/$(LIBTYPE)/$(WINLIBFILENAME) \
			$(SCIPLIBBASEOBJFILES) $(OBJSCIPLIBOBJFILES) $(LPILIBOBJFILES) $(TPILIBOBJFILES) $(SYMOBJFILES) \
			$(LPSLDFLAGS) $(LDFLAGS)
else
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LINKCC_L)$(LIBDIR)/$(LIBTYPE) $(LIBBUILD_o)$(LIBDIR)/$(LIBTYPE)/$(WINLIBFILENAME) \
			$(SCIPLIBBASEOBJFILES) $(OBJSCIPLIBOBJFILES) $(LPILIBOBJFILES) $(TPILIBOBJFILES) $(SYMOBJFILES) \
			$(LPSLDFLAGS) $(LDFLAGS)
endif

.PHONY: touchexternal
touchexternal:	$(ZLIBDEP) $(GMPDEP) $(READLINEDEP) $(ZIMPLDEP) $(AMPLDEP) $(LPSCHECKDEP) $(THREADSAFEDEP) | $(LIBOBJDIR)
ifeq ($(TOUCHLINKS),true)
		@-touch $(ZLIBSRC)
		@-touch $(GMPSRC)
		@-touch $(READLINESRC)
		@-touch $(ZIMPLSRC)
		@-touch $(AMPLSRC)
		@-touch $(LPSCHECKSRC)
		@-touch $(LPILIBSRC)
endif
ifneq ($(SCIPGITHASH),$(LAST_SCIPGITHASH))
		@$(MAKE) githash
endif
		@$(SHELL) -ec 'if test ! -e $(SCIPGITHASHFILE) ; \
			then \
				echo "-> generating $(SCIPGITHASHFILE)" ; \
				$(MAKE) githash ; \
			fi'
ifneq ($(subst \\n,\n,$(BUILDFLAGS)),$(LAST_BUILDFLAGS))
		@echo "#define SCIP_BUILDFLAGS \"$(BUILDFLAGS)\"" > $(SCIPBUILDFLAGSFILE)
endif
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
ifneq ($(AMPL),$(LAST_AMPL))
		@-touch $(AMPLSRC)
endif
ifneq ($(SYM),$(LAST_SYM))
		@-touch $(SYMSRC)
endif
ifneq ($(LPSCHECK),$(LAST_LPSCHECK))
		@-touch $(LPSCHECKSRC)
endif
ifneq ($(THREADSAFE),$(LAST_THREADSAFE))
		@-touch $(THREADSAFESRC)
endif
ifneq ($(USRFLAGS),$(LAST_USRFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USROFLAGS),$(LAST_USROFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRCFLAGS),$(LAST_USRCFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRCXXFLAGS),$(LAST_USRCXXFLAGS))
		@-touch $(ALLSRC)
endif
ifneq ($(USRLDFLAGS),$(LAST_USRLDFLAGS))
		@-touch -c $(SCIPLIBBASEOBJFILES) $(LPILIBOBJFILES) $(TPILIBOBJFILES) $(TPILIBOBJFILES) $(MAINOBJFILES)
endif
ifneq ($(USRARFLAGS),$(LAST_USRARFLAGS))
		@-touch -c $(SCIPLIBBASEOBJFILES) $(OBJSCIPLIBOBJFILES) $(LPILIBOBJFILES) $(TPILIBOBJFILES)
endif
ifneq ($(NOBLKMEM),$(LAST_NOBLKMEM))
		@-touch -c $(ALLSRC)
endif
ifneq ($(NOBUFMEM),$(LAST_NOBUFMEM))
		@-touch -c $(ALLSRC)
endif
ifneq ($(NOBLKBUFMEM),$(LAST_NOBLKBUFMEM))
		@-touch -c $(ALLSRC)
endif
ifneq ($(SANITIZE),$(LAST_SANITIZE))
		@-touch -c $(ALLSRC)
endif
ifneq ($(TPI),$(LAST_TPI))
		@-touch -c $(ALLSRC)
endif
ifneq ($(DEBUGSOL),$(LAST_DEBUGSOL))
		@-touch -c $(ALLSRC)
endif
ifneq ($(PAPILO),$(LAST_PAPILO))
		@-touch -c $(ALLSRC)
endif
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_BUILDFLAGS=\"$(BUILDFLAGS)\"" >> $(LASTSETTINGS)
		@echo "LAST_SCIPGITHASH=$(SCIPGITHASH)" >> $(LASTSETTINGS)
		@echo "LAST_ZLIB=$(ZLIB)" >> $(LASTSETTINGS)
		@echo "LAST_GMP=$(GMP)" >> $(LASTSETTINGS)
		@echo "LAST_READLINE=$(READLINE)" >> $(LASTSETTINGS)
		@echo "LAST_ZIMPL=$(ZIMPL)" >> $(LASTSETTINGS)
		@echo "LAST_AMPL=$(AMPL)" >> $(LASTSETTINGS)
		@echo "LAST_IPOPT=$(IPOPT)" >> $(LASTSETTINGS)
		@echo "LAST_WORHP=$(WORHP)" >> $(LASTSETTINGS)
		@echo "LAST_SYM=$(SYM)" >> $(LASTSETTINGS)
		@echo "LAST_THREADSAFE=$(THREADSAFE)" >> $(LASTSETTINGS)
		@echo "LAST_LPSCHECK=$(LPSCHECK)" >> $(LASTSETTINGS)
		@echo "LAST_USRFLAGS=$(USRFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USROFLAGS=$(USROFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRCFLAGS=$(USRCFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRCXXFLAGS=$(USRCXXFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRLDFLAGS=$(USRLDFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_USRARFLAGS=$(USRARFLAGS)" >> $(LASTSETTINGS)
		@echo "LAST_NOBLKMEM=$(NOBLKMEM)" >> $(LASTSETTINGS)
		@echo "LAST_NOBUFMEM=$(NOBUFMEM)" >> $(LASTSETTINGS)
		@echo "LAST_NOBLKBUFMEM=$(NOBLKBUFMEM)" >> $(LASTSETTINGS)
		@echo "LAST_SANITIZE=$(SANITIZE)" >> $(LASTSETTINGS)
		@echo "LAST_TPI=$(TPI)" >> $(LASTSETTINGS)
		@echo "LAST_DEBUGSOL=$(DEBUGSOL)" >> $(LASTSETTINGS)
		@echo "LAST_PAPILO=$(PAPILO)" >> $(LASTSETTINGS)

$(LINKSMARKERFILE):
		@$(MAKE) links

.PHONY: links
links:		| $(LIBDIR)/static $(LIBDIR)/shared $(LIBDIR)/include $(DIRECTORIES) echosoftlinks $(SOFTLINKS)
		@rm -f $(LINKSMARKERFILE)
		@echo "this is only a marker" > $(LINKSMARKERFILE)

.PHONY: echosoftlinks
echosoftlinks:
		@echo
		@echo "- Current settings: LPS=$(LPS) OSTYPE=$(OSTYPE) ARCH=$(ARCH) COMP=$(COMP) SHARED=$(SHARED) SUFFIX=$(LINKLIBSUFFIX) ZIMPL=$(ZIMPL) ZIMPLOPT=$(ZIMPLOPT) AMPL=$(AMPL) IPOPT=$(IPOPT) IPOPTOPT=$(IPOPTOPT) WORHP=$(WORHP) WORHPOPT=$(WORHPOPT) FILTERSQP=$(FILTERSQP) EXPRINT=$(EXPRINT) SYM=$(SYM)"
		@echo
		@echo "* SCIP needs some softlinks to external programs, in particular, LP-solvers."
		@echo "* Please insert the paths to the corresponding directories/libraries below."
		@echo "* The links will be installed in the 'lib/include' and 'lib/$(LIBTYPE)' directories."
		@echo "* For more information and if you experience problems see the INSTALL file."
		@echo
		@echo -e $(LPIINSTMSG)

$(DIRECTORIES):
		@echo
		@echo "- creating directory \"$@\""
		@-mkdir -p $@


# Create softlinks for external libraries. The user can enter the
# filename and the link is created in the corresponding directories
# (lib/include, lib/shared, lib/static).
#
# We also test whether older links in lib/ still exist. In later
# releases this feature can be removed.
.PHONY: $(SOFTLINKS)
$(SOFTLINKS):
ifeq ($(MAKESOFTLINKS), true)
		@$(SHELL) -ec 'if test ! -e $@ ; \
			then \
				DIRNAME=`dirname $@` ; \
				echo ; \
				echo "> Enter soft-link target file or directory for \"$@\" (return if not needed): " ; \
				echo -n "> " ; \
				cd $$DIRNAME ; \
				eval $(READ) TARGET ; \
				cd $(SCIPREALPATH) ; \
				if test "$$TARGET" != "" ; \
				then \
					echo "-> creating softlink \"$@\" -> \"$$TARGET\"" ; \
					rm -f $@ ; \
					$(LN_s) $$TARGET $@ ; \
				else \
					echo "* skipped creation of softlink \"$@\". Call \"make links\" if needed later." ; \
				fi ; \
				FILENAME=$@ ; \
				FNAME=$${FILENAME#lib/static} ; \
				FNAME=$${FNAME#lib/include} ; \
				FNAME="lib"$$FNAME ; \
				if test -e $$FNAME ; \
				then \
					echo ; \
					echo "The link "$$FNAME" still exists. Consider removing it." ; \
				fi ; \
				echo ; \
			fi'
endif

.PHONY: checkdefines
checkdefines:
ifeq ($(LPILIBOBJ),)
		$(error invalid LP solver selected: LPS=$(LPS). Possible options are: $(LPSOPTIONS))
endif
ifneq ($(TPI),none)
ifneq ($(TPI),omp)
ifneq ($(TPI),tny)
		$(error invalid TPI flag selected: TPI=$(TPI). Possible options are: $(TPIOPTIONS))
endif
endif
endif
ifneq ($(GMP),true)
ifneq ($(GMP),false)
		$(error invalid GMP flag selected: GMP=$(GMP). Possible options are: true false)
endif
endif
ifneq ($(ZIMPL),true)
ifneq ($(ZIMPL),false)
		$(error invalid ZIMPL flag selected: ZIMPL=$(ZIMPL). Possible options are: true false auto)
endif
endif
ifneq ($(AMPL),true)
ifneq ($(AMPL),false)
		$(error invalid AMPL flag selected: AMPL=$(AMPL). Possible options are: true false)
endif
endif
ifneq ($(IPOPT),true)
ifneq ($(IPOPT),false)
		$(error invalid IPOPT flag selected: IPOPT=$(IPOPT). Possible options are: true false)
endif
endif
ifneq ($(FILTERSQP),true)
ifneq ($(FILTERSQP),false)
		$(error invalid FILTERSQP flag selected: FILTERSQP=$(FILTERSQP). Possible options are: true false)
endif
endif
ifneq ($(WORHP),true)
ifneq ($(WORHP),false)
		$(error invalid WORHP flag selected: WORHP=$(WORHP). Possible options are: true false)
endif
endif
ifneq ($(READLINE),true)
ifneq ($(READLINE),false)
		$(error invalid READLINE flag selected: READLINE=$(READLINE). Possible options are: true false)
endif
endif
ifneq ($(ZLIB),true)
ifneq ($(ZLIB),false)
		$(error invalid ZLIB flag selected: ZLIB=$(ZLIB). Possible options are: true false)
endif
endif
ifneq ($(THREADSAFE),true)
ifneq ($(THREADSAFE),false)
		$(error invalid THREADSAFE flag selected: THREADSAFE=$(THREADSAFE). Possible options are: true false)
endif
endif
ifeq ($(SHARED),true)
ifeq ($(COMP),msvc)
		$(error invalid flags selected: SHARED=$(SHARED) and COMP=$(COMP). Please use 'make dll' to generate a dynamic library with MSVC)
endif
endif
ifneq ($(SYM),bliss)
ifneq ($(SYM),sbliss)
ifneq ($(SYM),nauty)
ifneq ($(SYM),none)
		$(error invalid SYM flag selected: SYM=$(SYM). Possible options are: $(SYMOPTIONS))
endif
endif
endif
endif
ifneq ($(PAPILO),true)
ifneq ($(PAPILO),false)
		$(error invalid PAPILO flag selected: PAPILO=$(PAPILO). Possible options are: true false)
endif
endif

.PHONY: errorhints
errorhints:
ifeq ($(READLINE),true)
		@echo "build failed with READLINE=true: if readline is not available, try building with READLINE=false"
endif
ifeq ($(ZLIB),true)
		@echo "build failed with ZLIB=true: if ZLIB is not available, try building with ZLIB=false"
endif
ifeq ($(GMP),true)
		@echo "build failed with GMP=true: if GMP is not available, try building with GMP=false (note that this will deactivate Zimpl support)"
endif
ifeq ($(GMP),false)
ifeq ($(LPS),spx1)
		@echo "build failed with GMP=false and LPS=spx1: use GMP=true or make sure that SoPlex is also built without GMP support (make GMP=false)"
endif
ifeq ($(LPS),spx2)
		@echo "build failed with GMP=false and LPS=spx2: use GMP=true or make sure that SoPlex is also built without GMP support (make GMP=false)"
endif
endif

.PHONY: help
help:
		@echo "Use the SCIP makefile system."
		@echo
		@echo "  The main options for the SCIP makefile system are as follows:"
		@echo
		@echo "  General options:"
		@echo "  - OPT={dbg|opt}: Use debug or optimized (default) mode, respectively."
		@echo "  - LPS={clp|cpx|grb|glop|msk|qso|spx|xprs|none}: Determine LP-solver."
		@echo "      clp: COIN-OR Clp LP-solver"
		@echo "      cpx: CPLEX LP-solver"
		@echo "      glop: Glop LP-solver"
		@echo "      grb: Gurobi LP-solver"
		@echo "      msk: Mosek LP-solver"
		@echo "      qso: QSopt LP-solver"
		@echo "      spx: old SoPlex LP-solver (for versions < 2)"
		@echo "      spx2: new SoPlex LP-solver (default) (from version 2)"
		@echo "      xprs: XPress LP-solver"
		@echo "      none: no LP-solver"
		@echo "  - COMP={clang|gnu|intel}: Determine compiler."
		@echo "  - TPI={none|omp|tny}: Determine parallel interface for concurrent solve (default: none)."
		@echo "  - SHARED={true|false}: Build shared libraries or not (default)."
		@echo
		@echo "  More detailed options:"
		@echo "  - ZIMPL=<true|false>: Turn ZIMPL support on or off (default)."
		@echo "  - ZIMPLOPT=<dbg|opt>: Use debug or optimized (default) mode for ZIMPL."
		@echo "  - AMPL=<true|false>: Turn AMPL .nl support on (default) or off."
		@echo "  - LPSOPT=<dbg|opt>: Use debug or optimized (default) mode for LP-solver (SoPlex and Clp only)."
		@echo "  - READLINE=<true|false>: Turns support via the readline library on (default) or off."
		@echo "  - IPOPT=<true|false>: Turns support of IPOPT on or off (default)."
		@echo "  - EXPRINT=<cppad|none>: Use CppAD as expressions interpreter (default) or no expressions interpreter."
		@echo "  - SYM=<none|bliss|nauty|sbliss>: To choose type of symmetry handling."
		@echo "  - PARASCIP=<true|false>: Build for ParaSCIP (deprecated, use THREADSAFE)."
		@echo "  - THREADSAFE=<true|false>: Build thread safe."
		@echo "  - NOBLKMEM=<true|false>: Turn off block memory or on (default)."
		@echo "  - NOBUFMEM=<true|false>>: Turn off buffer memory or on (default)."
		@echo "  - NOBLKBUFMEM=<true|false>: Turn usage of internal memory functions off or on (default)."
		@echo "  - VERBOSE=<true|false>: Turn on verbose messages of makefile or off (default)."
		@echo
		@echo "  Main targets:"
		@echo "  - all (default): Build SCIP libraries and binary."
		@echo "  - libscip: Build standalone SCIP library."
		@echo "  - libscipbase: Build library for the main part of SCIP."
		@echo "  - libobjscip: Build library for the C++-interface of SCIP."
		@echo "  - liblpi: Build library for the LP interface in SCIP."
		@echo "  - libtpi: Build library for the parallel task interface in SCIP."
		@echo "  - links: Reconfigures the links in the \"lib\" directory."
		@echo "  - doc: Creates documentation in the \"doc\" directory."
		@echo "  - libs: Create all SCIP libraries."
		@echo "  - lint: Run lint on all SCIP files. (Needs flexelint.)"
		@echo "  - pclint: Run pclint on all SCIP files. (Needs pclint.)"
		@echo "  - splint: Run splint on all C SCIP files. (Needs splint.)"
		@echo "  - clean: Removes all object files."
		@echo "  - cleanlibs: Remove all SCIP libraries."
		@echo "  - depend: Updates dependencies files. This is only needed if you modify SCIP source."
		@echo "  - tags: Creates TAGS file that can be used in (x)emacs."
		@echo "  - check or test: Runs the check/test script, see the online documentation."

# --- EOF ---------------------------------------------------------------------
# DO NOT DELETE
