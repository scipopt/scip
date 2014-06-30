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
#@brief   Makefile for Java MIP solver example
#@author  Stefan Heinz


#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

JNISCIPDIR	=	../..
VERBOSE		=	false

#-----------------------------------------------------------------------------
# Java stettings
#-----------------------------------------------------------------------------

JAVAC		=	javac

CLASSDIR	=	classes
JAVASRCDIR	=	java


#-----------------------------------------------------------------------------
# SCIP JNI Library
#-----------------------------------------------------------------------------

MAINOBJ		=	JniKnapsack.java
JNICLASSES	=	$(addprefix $(CLASSDIR)/,$(MAINOBJ:.java=.class))

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(JNICLASSES)
endif

.PHONY: all
all:            $(JNICLASSES)

.PHONY: run
run:
		./run.sh

.PHONY: clean
clean:
		@-rm -f $(CLASSDIR)/*.class
		@echo "-> remove java class"

$(CLASSDIR):
		@-mkdir -p $(CLASSDIR)

$(CLASSDIR)/%.class:	$(JAVASRCDIR)/%.java $(CLASSDIR)
		@echo "-> compiling $@"
		$(JAVAC) -d $(CLASSDIR) -classpath $(JNISCIPDIR)/lib/scip.jar $<

#---- EOF --------------------------------------------------------------------
