#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic Licence.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#@file    Makefile
#@brief   Makefile for cumulative scheduling example
#@author  Stefan Heinz


#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

SCIPDIR         =       ../..
CPOPT		=	false

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	scip
MAINOBJ		=	cons_optcumulative.o \
			heur_listscheduling.o \
			heur_optcumulative.o \
			reader_cmin.o \
			reader_sch.o \
			reader_sm.o \
			reader_rcp.o

MAINCXXOBJ	= 	main.o \
			cpoptimizer.o

MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.c))
MAINSRC		+=	$(addprefix $(SRCDIR)/,$(MAINCXXOBJ:.o=.cpp))
MAINDEP		=	$(SRCDIR)/depend.main.$(OPT)

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))
MAINOBJFILES	+=	$(addprefix $(OBJDIR)/,$(MAINCXXOBJ))

ifeq ($(CPOPT),true)
FLAGS		+=	-DWITH_CPOPTIMIZER -DIL_STD -I$(LIBDIR)/cpinc -I$(LIBDIR)/cpxinc
FLAGS		+=	-I$(LIBDIR)/ilconcertinc
LDFLAGS         +=       $(LIBDIR)/libcp.a $(LIBDIR)/libconcert.a  $(LIBDIR)/libcplex.a -lpthread
endif

LASTSETTINGS	=	$(OBJDIR)/make.lastsettings

CPOPTDEP	:=	$(SRCDIR)/depend.cpoptimizer
CPOPTSRC	:=	$(shell cat $(CPOPTDEP))

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) -I$(SCIPDIR) lint/main-gcc.lnt +os\(lint.out\) -u -zero \
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

.PHONY: scip
scip:
		@$(MAKE) -C $(SCIPDIR) libs $^
		@$(MAKE) $^

.PHONY: doc
doc:
		@-(cd doc && $(LN_s) ../$(SCIPDIR)/doc/scip.css);
		@-(cd doc && $(LN_s) ../$(SCIPDIR)/doc/pictures/scippy.png);
		@-(cd doc && $(LN_s) ../$(SCIPDIR)/doc/pictures/miniscippy.png);
		@-(cd doc && $(LN_s) ../$(SCIPDIR)/doc/scipfooter.html footer.html);
		cd doc; $(DOXY) scheduler.dxy

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-rm -f $(LASTSETTINGS)
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

$(CHECKDIR):
		@-mkdir -p $(CHECKDIR)

.PHONY: test
test:		$(MAINFILE) $(CHECKDIR)
		@-(cd check && ln -fs ../$(SCIPDIR)/check/check.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck_cluster.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/check.awk);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/getlastprob.awk);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/configuration_set.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/configuration_logfiles.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/configuration_tmpfile_setup_scip.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/run.sh);
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) \
		$(CONTINUE) $(LOCK) "example" $(LPS) $(VALGRIND) $(CLIENTTMPDIR) $(REOPT) $(OPTCOMMAND) $(SETCUTOFF) $(MAXJOBS);

# include local targets
-include $(SCIPDIR)/make/local/make.targets

.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		@-rm -f $(LASTSETTINGS)
		@-(rm -f $(OBJDIR)/*.o && rmdir $(OBJDIR));
		@echo "-> remove main objective files"
endif
		@-rm -f $(MAINFILE) $(MAINLINK) $(MAINSHORTLINK)
		@echo "-> remove binary"

.PHONY: tags
tags:
		rm -f TAGS; ctags -e src/*.c src/*.h $(SCIPDIR)/src/scip/*.c $(SCIPDIR)/src/scip/*.h;

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).c|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(MAINDEP)'
		@echo `grep -l "WITH_CPOPTIMIZER" $(MAINSRC)` >$(CPOPTDEP)


-include $(LASTSETTINGS)

.PHONY: touchexternal
touchexternal:	$(CPOPTDEP)
ifneq ($(CPOPT),$(LAST_CPOPT))
		@-touch $(CPOPTSRC)
endif
		@-rm -f $(LASTSETTINGS)
		@echo "LAST_CPOPT=$(CPOPT)" >> $(LASTSETTINGS)


-include	$(MAINDEP)

$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) touchexternal $(MAINOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
