#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id$

VERBOSE		=	false
CXX		=	clang++
CXX_o		=	-o # the trailing space is important
LINKCXX		=	clang++
LINKCXX_L	=	-L
LINKCXX_l	=	-l
LINKCXX_o	=	-o # the trailing space is important
FLAGS		=	
OFLAGS		=	
CXXFLAGS	=	-I/usr/local/include
LDFLAGS		=	-L/usr/local/lib

DOXY		=	doxygen

#-----------------------------------------------------------------------------
# GMP Library
#-----------------------------------------------------------------------------
GMP_LDFLAGS = $(LINKCXX_l)gmp

#-----------------------------------------------------------------------------
# ZLIB Library
#-----------------------------------------------------------------------------
ZLIB_LDFLAGS = $(LINKCXX_l)z

#-----------------------------------------------------------------------------
# Local Setting (include local make file setting if these exist)
#-----------------------------------------------------------------------------
-include make.local

#-----------------------------------------------------------------------------
# USER Flags
#-----------------------------------------------------------------------------
FLAGS		+=	$(USRFLAGS)
OFLAGS		+=	$(USROFLAGS)
CXXFLAGS	+=	$(USRCXXFLAGS)
LDFLAGS		+=	$(USRLDFLAGS)

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------
MAINNAME	=	solchecker

MAINOBJ		=	gmputils.o \
			main.o \
			model.o \
			mpsinput.o 

BINDIR          =       bin
OBJDIR          =       obj
SRCDIR		= 	src

MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))

MAIN		=	$(MAINNAME)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))



#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES)
endif

.PHONY: all
all:            $(MAINFILE) $(MAINOBJFILES)

$(OBJDIR):	
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:          $(BINDIR) $(OBJDIR)
ifneq ($(OBJDIR),)
		@echo "-> remove objective files"
		@-rm -f $(OBJDIR)/*.o
		@-rmdir $(OBJDIR)
endif
		@echo "-> remove binary"
		@-rm -f $(MAINFILE)
		@-rmdir $(BINDIR)

.PHONY: doc
doc: 
		cd doc; $(DOXY) $(MAINNAME).dxy;

$(MAINFILE):	$(BINDIR) $(OBJDIR) $(MAINOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) $(LDFLAGS) $(ZLIB_LDFLAGS) $(GMP_LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
