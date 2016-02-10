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
#@brief   Makefile for SCIP JNI library
#@author  Stefan Heinz


VERSION 	= 	0.1

#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------

SCIPDIR         =       lib/scip
SOPLEXDIR	=       lib/soplex

#-----------------------------------------------------------------------------
# Force certain variable to specific values
#-----------------------------------------------------------------------------

override SHARED	=	true
override READLINE=	false
override ZIMPL	=	false
override ZLIB	=	false

#-----------------------------------------------------------------------------
# Java stettings
#-----------------------------------------------------------------------------

JAVAC		=	javac
JAVAH		=	javah
JAR		=	jar

CLASSDIR	=	classes
JAVASRCDIR	=	java

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project

FLAGS		+=	-I$(LIBDIR)/jniinc

#-----------------------------------------------------------------------------
# SCIP Java Library
#-----------------------------------------------------------------------------

JAVAPACKAGENAME	=	scip
JAVAPACKAGEJAR	=	$(LIBDIR)/$(JAVAPACKAGENAME).jar

#-----------------------------------------------------------------------------
# SCIP JNI Library
#-----------------------------------------------------------------------------

JNIPACKAGENAME	=	de.zib.jscip.nativ
JNILIBSHORTNAME	=	jscip
JNIPACKAGEJAR	=	$(LIBDIR)/$(JNILIBSHORTNAME).jar

JNILIBNAME	=	$(JNILIBSHORTNAME)-$(VERSION)

JNILIBOBJ	=	JniScipLibraryLoader.java \
			JniScipBranchAllfullstrong.java \
			JniScipBranchFullstrong.java \
			JniScipBranchInference.java \
			JniScipBranch.java \
			JniScipBranchLeastinf.java \
			JniScipBranchMostinf.java \
			JniScipBranchPscost.java \
			JniScipBranchRandom.java \
			JniScipBranchRelpscost.java \
			JniScipConflict.java \
			JniScipConsAbspower.java \
			JniScipConsAnd.java \
			JniScipConsBivariate.java \
			JniScipConsBounddisjunction.java \
			JniScipConsConjunction.java \
			JniScipConsCountsols.java \
			JniScipConsCumulative.java \
			JniScipConsDisjunction.java \
			JniScipConsIndicator.java \
			JniScipConsIntegral.java \
			JniScipCons.java \
			JniScipConsKnapsack.java \
			JniScipConsLinear.java \
			JniScipConsLinking.java \
			JniScipConsLogicor.java \
			JniScipConsNonlinear.java \
			JniScipConsOrbitope.java \
			JniScipConsOr.java \
			JniScipConsPseudoboolean.java \
			JniScipConsQuadratic.java \
			JniScipConsSetppc.java \
			JniScipConsSoc.java \
			JniScipConsSos1.java \
			JniScipConsSos2.java \
			JniScipConsSuperindicator.java \
			JniScipConsVarbound.java \
			JniScipConsXor.java \
			JniScipCutpool.java \
			JniScipDialogDefault.java \
			JniScipDialog.java \
			JniScipDispDefault.java \
			JniScipDisp.java \
			JniScipEvent.java \
			JniScipExpr.java \
			JniScipHeurActconsdiving.java \
			JniScipHeurClique.java \
			JniScipHeurCoefdiving.java \
			JniScipHeurCrossover.java \
			JniScipHeurDins.java \
			JniScipHeurFeaspump.java \
			JniScipHeurFixandinfer.java \
			JniScipHeurFracdiving.java \
			JniScipHeurGuideddiving.java \
			JniScipHeurIntdiving.java \
			JniScipHeurIntshifting.java \
			JniScipHeur.java \
			JniScipHeurLinesearchdiving.java \
			JniScipHeurLocalbranching.java \
			JniScipHeurMutation.java \
			JniScipHeurNlpdiving.java \
			JniScipHeurObjpscostdiving.java \
			JniScipHeurOctane.java \
			JniScipHeurOneopt.java \
			JniScipHeurPscostdiving.java \
			JniScipHeurRens.java \
			JniScipHeurRins.java \
			JniScipHeurRootsoldiving.java \
			JniScipHeurRounding.java \
			JniScipHeurShiftandpropagate.java \
			JniScipHeurShifting.java \
			JniScipHeurSimplerounding.java \
			JniScipHeurSubnlp.java \
			JniScipHeurTrivial.java \
			JniScipHeurTrysol.java \
			JniScipHeurTwoopt.java \
			JniScipHeurUndercover.java \
			JniScipHeurVbounds.java \
			JniScipHeurVeclendiving.java \
			JniScipHeurZeroobj.java \
			JniScipHeurZirounding.java \
			JniScipImplics.java \
			JniScip.java \
			JniScipLp.java \
			JniScipMessageDefault.java \
			JniScipMessage.java \
			JniScipMisc.java \
			JniScipNlpiIpopt.java \
			JniScipNlp.java \
			JniScipNodeselBfs.java \
			JniScipNodeselDfs.java \
			JniScipNodeselEstimate.java \
			JniScipNodeselHybridestim.java \
			JniScipNodesel.java \
			JniScipNodeselRestartdfs.java \
			JniScipParamset.java \
			JniScipPresolBoundshift.java \
			JniScipPresolComponents.java \
			JniScipPresolConvertinttobin.java \
			JniScipPresolDomcol.java \
			JniScipPresolGateextraction.java \
			JniScipPresolImplics.java \
			JniScipPresolInttobinary.java \
			JniScipPresol.java \
			JniScipPresolTrivial.java \
			JniScipPricer.java \
			JniScipPropDualfix.java \
			JniScipPropGenvbounds.java \
			JniScipProp.java \
			JniScipPropObbt.java \
			JniScipPropProbing.java \
			JniScipPropPseudoobj.java \
			JniScipPropRedcost.java \
			JniScipPropRootredcost.java \
			JniScipPropVbounds.java \
			JniScipReaderBnd.java \
			JniScipReaderCcg.java \
			JniScipReaderCip.java \
			JniScipReaderCnf.java \
			JniScipReaderFix.java \
			JniScipReaderFzn.java \
			JniScipReaderGms.java \
			JniScipReader.java \
			JniScipReaderLp.java \
			JniScipReaderMps.java \
			JniScipReaderOpb.java \
			JniScipReaderOsil.java \
			JniScipReaderPip.java \
			JniScipReaderPpm.java \
			JniScipReaderRlp.java \
			JniScipReaderSol.java \
			JniScipReaderWbo.java \
			JniScipReaderZpl.java \
			JniScipRelax.java \
			JniScipSepaCgmip.java \
			JniScipSepaClique.java \
			JniScipSepaClosecuts.java \
			JniScipSepaCmir.java \
			JniScipSepaFlowcover.java \
			JniScipSepaGomory.java \
			JniScipSepaImpliedbounds.java \
			JniScipSepaIntobj.java \
			JniScipSepa.java \
			JniScipSepaMcf.java \
			JniScipSepaOddcycle.java \
			JniScipSepaRapidlearning.java \
			JniScipSepaStrongcg.java \
			JniScipSepaZerohalf.java \
			JniScipSol.java \
			JniScipTree.java \
			JniScipVar.java

JNIENUMS	=	JniScipBoundchgtype.java \
			JniScipBoundtype.java \
			JniScipBranchdir.java \
			JniScipObjsense.java \
			JniScipParamemphasis.java \
			JniScipParamsetting.java \
			JniScipParamtype.java \
			JniScipResult.java \
			JniScipRetcode.java \
			JniScipSolorigin.java \
			JniScipStage.java \
			JniScipStatus.java \
			JniScipVarstatus.java \
			JniScipVartype.java


NATIVELIBOBJ	= 	NativeScipBranchAllfullstrong.java \
			NativeScipBranchFullstrong.java \
			NativeScipBranchInference.java \
			NativeScipBranch.java \
			NativeScipBranchLeastinf.java \
			NativeScipBranchMostinf.java \
			NativeScipBranchPscost.java \
			NativeScipBranchRandom.java \
			NativeScipBranchRelpscost.java \
			NativeScipConflict.java \
			NativeScipConsAbspower.java \
			NativeScipConsAnd.java \
			NativeScipConsBivariate.java \
			NativeScipConsBounddisjunction.java \
			NativeScipConsConjunction.java \
			NativeScipConsCountsols.java \
			NativeScipConsCumulative.java \
			NativeScipConsDisjunction.java \
			NativeScipConsIndicator.java \
			NativeScipConsIntegral.java \
			NativeScipCons.java \
			NativeScipConsKnapsack.java \
			NativeScipConsLinear.java \
			NativeScipConsLinking.java \
			NativeScipConsLogicor.java \
			NativeScipConsNonlinear.java \
			NativeScipConsOrbitope.java \
			NativeScipConsOr.java \
			NativeScipConsPseudoboolean.java \
			NativeScipConsQuadratic.java \
			NativeScipConsSetppc.java \
			NativeScipConsSoc.java \
			NativeScipConsSos1.java \
			NativeScipConsSos2.java \
			NativeScipConsSuperindicator.java \
			NativeScipConsVarbound.java \
			NativeScipConsXor.java \
			NativeScipCutpool.java \
			NativeScipDialogDefault.java \
			NativeScipDialog.java \
			NativeScipDispDefault.java \
			NativeScipDisp.java \
			NativeScipEvent.java \
			NativeScipException.java \
			NativeScipExpr.java \
			NativeScipHeurActconsdiving.java \
			NativeScipHeurClique.java \
			NativeScipHeurCoefdiving.java \
			NativeScipHeurCrossover.java \
			NativeScipHeurDins.java \
			NativeScipHeurFeaspump.java \
			NativeScipHeurFixandinfer.java \
			NativeScipHeurFracdiving.java \
			NativeScipHeurGuideddiving.java \
			NativeScipHeurIntdiving.java \
			NativeScipHeurIntshifting.java \
			NativeScipHeur.java \
			NativeScipHeurLinesearchdiving.java \
			NativeScipHeurLocalbranching.java \
			NativeScipHeurMutation.java \
			NativeScipHeurNlpdiving.java \
			NativeScipHeurObjpscostdiving.java \
			NativeScipHeurOctane.java \
			NativeScipHeurOneopt.java \
			NativeScipHeurPscostdiving.java \
			NativeScipHeurRens.java \
			NativeScipHeurRins.java \
			NativeScipHeurRootsoldiving.java \
			NativeScipHeurRounding.java \
			NativeScipHeurShiftandpropagate.java \
			NativeScipHeurShifting.java \
			NativeScipHeurSimplerounding.java \
			NativeScipHeurSubnlp.java \
			NativeScipHeurTrivial.java \
			NativeScipHeurTrysol.java \
			NativeScipHeurTwoopt.java \
			NativeScipHeurUndercover.java \
			NativeScipHeurVbounds.java \
			NativeScipHeurVeclendiving.java \
			NativeScipHeurZeroobj.java \
			NativeScipHeurZirounding.java \
			NativeScipImplics.java \
			NativeScip.java \
			NativeScipLp.java \
			NativeScipMessageDefault.java \
			NativeScipMessage.java \
			NativeScipMisc.java \
			NativeScipNlpiIpopt.java \
			NativeScipNlp.java \
			NativeScipNodeselBfs.java \
			NativeScipNodeselDfs.java \
			NativeScipNodeselEstimate.java \
			NativeScipNodeselHybridestim.java \
			NativeScipNodesel.java \
			NativeScipNodeselRestartdfs.java \
			NativeScipParamset.java \
			NativeScipPresolBoundshift.java \
			NativeScipPresolComponents.java \
			NativeScipPresolConvertinttobin.java \
			NativeScipPresolDomcol.java \
			NativeScipPresolGateextraction.java \
			NativeScipPresolImplics.java \
			NativeScipPresolInttobinary.java \
			NativeScipPresol.java \
			NativeScipPresolTrivial.java \
			NativeScipPricer.java \
			NativeScipPropDualfix.java \
			NativeScipPropGenvbounds.java \
			NativeScipProp.java \
			NativeScipPropObbt.java \
			NativeScipPropProbing.java \
			NativeScipPropPseudoobj.java \
			NativeScipPropRedcost.java \
			NativeScipPropRootredcost.java \
			NativeScipPropVbounds.java \
			NativeScipReaderBnd.java \
			NativeScipReaderCcg.java \
			NativeScipReaderCip.java \
			NativeScipReaderCnf.java \
			NativeScipReaderFix.java \
			NativeScipReaderFzn.java \
			NativeScipReaderGms.java \
			NativeScipReader.java \
			NativeScipReaderLp.java \
			NativeScipReaderMps.java \
			NativeScipReaderOpb.java \
			NativeScipReaderOsil.java \
			NativeScipReaderPip.java \
			NativeScipReaderPpm.java \
			NativeScipReaderRlp.java \
			NativeScipReaderSol.java \
			NativeScipReaderWbo.java \
			NativeScipReaderZpl.java \
			NativeScipRelax.java \
			NativeScipSepaCgmip.java \
			NativeScipSepaClique.java \
			NativeScipSepaClosecuts.java \
			NativeScipSepaCmir.java \
			NativeScipSepaFlowcover.java \
			NativeScipSepaGomory.java \
			NativeScipSepaImpliedbounds.java \
			NativeScipSepaIntobj.java \
			NativeScipSepa.java \
			NativeScipSepaMcf.java \
			NativeScipSepaOddcycle.java \
			NativeScipSepaRapidlearning.java \
			NativeScipSepaStrongcg.java \
			NativeScipSepaZerohalf.java \
			NativeScipSol.java \
			NativeScipTree.java \
			NativeScipVar.java

JNIPACKAGEDIR	=	$(subst .,/,$(JNIPACKAGENAME))
JNILIB		=	$(JNILIBNAME).$(BASE).$(LPS)
JNILIBFILE	=	$(LIBDIR)/lib$(JNILIB).$(LIBEXT)
JNILIBOBJFILES	=	$(addprefix $(OBJDIR)/,$(JNILIBOBJ:.java=.o))
JNILIBSRC	=	$(addprefix $(SRCDIR)/,$(JNILIBOBJ:.java=.c))
JNILIBHEADERS	=	$(addprefix $(SRCDIR)/,$(JNILIBOBJ:.java=.h))
JNILIBHEADERS	+=	$(addprefix $(SRCDIR)/,$(JNIENUMS:.java=.h))
JNILIBDEP	=	$(SRCDIR)/depend.lib$(JNILIBSHORTNAME).$(LPS).$(OPT)
JNILIBLINK	=	$(LIBDIR)/lib$(JNILIBSHORTNAME).$(BASE).$(LIBEXT)
JNILIBSHORTLINK =	$(LIBDIR)/lib$(JNILIBSHORTNAME).$(LIBEXT)

JNICLASSES	=	$(addprefix $(CLASSDIR)/$(JNIPACKAGEDIR)/jni/,$(JNILIBOBJ:.java=.class))
JNICLASSES	+=	$(addprefix $(CLASSDIR)/$(JNIPACKAGEDIR)/jni/,$(JNIENUMS:.java=.class))

JARCLASSES	=	$(addprefix -C $(CLASSDIR) $(JNIPACKAGEDIR)/jni/,$(JNILIBOBJ:.java=.class))
JARCLASSES	+=	$(addprefix -C $(CLASSDIR) $(JNIPACKAGEDIR)/jni/,$(JNIENUMS:.java=.class))
JARCLASSES	+=	$(addprefix -C $(CLASSDIR) $(JNIPACKAGEDIR)/,$(NATIVELIBOBJ:.java=.class))

JNILIBRARYLOADER =	$(JAVASRCDIR)/$(JNIPACKAGEDIR)/jni/JniScipLibraryLoader.java

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(JNILIBRARYLOADER) $(JNICLASSES) $(JAVAPACKAGEJAR) \
		$(JNILIBHEADERS) $(JNILIBFILE) $(JNILIBOBJFILES) $(JNILIBLINK) $(JNILIBSHORTLINK) depend
endif

.PHONY: all
all:            $(JNILIBRARYLOADER) $(JNICLASSES) $(JAVAPACKAGEJAR) $(SCIPDIR) $(JNILIBFILE) $(JNILIBLINK) $(JNILIBSHORTLINK)


.PHONY: scip
scip:
		@echo "Create SCIP shared libraries"
		@$(MAKE) -C $(SCIPDIR) libs SHARED=$(SHARED) READLINE=$(READLINE) ZLIB=$(ZLIB) ZIMPL=$(ZIMPL) $^
		@echo "Install SCIP libraries into $(LIBDIR) directory"
		@-install -m 755 -d $(LIBDIR)/
		@-install -m 644 -p -t $(LIBDIR)/ $(LPILIBFILE) $(NLPILIBFILE) $(SCIPLIBFILE)

.PHONY: soplex
soplex:
		@echo "Create SoPlex shared libraries"
		@$(MAKE) -C $(SOPLEXDIR)  SHARED=$(SHARED) ZLIB=$(ZLIB) $^
		@echo "Install SoPlex library into $(LIBDIR) directory"
		@-install -m 755 -d $(LIBDIR)/
		@-install -m 644 -p -t $(LIBDIR) $(SOPLEXDIR)/$(LIBDIR)/libsoplex.$(BASE).so

.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		@-(rm -f $(OBJDIR)/*.o && rmdir $(OBJDIR));
		@echo "-> remove objective files"
endif
		@-rm -f $(JNILIBSHORTLINK) $(JNILIBLINK) $(JNILIBFILE)
		@echo "-> remove $(JNILIBFILE)"
		@-rm -f $(JAVAPACKAGEJAR)
		@echo "-> remove $(JAVAPACKAGEJAR)"
		@-rm -f $(JNICLASSES)
		@echo "-> remove java class"
		@-rm -f $(JNILIBHEADERS)
		@echo "-> remove c headers"

.PHONY: depend
depend:
		@echo "-> create dependencies"
		$(SHELL) -ec '$(DCC) $(FLAGS) $(DFLAGS) $(JNILIBSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z_/]*\).c|$$\(LIBOBJDIR\)/\2.o: $(SRCDIR)/\2.c|g'\'' \
		>$(JNILIBDEP)'

$(JNILIBSHORTLINK):	$(JNILIBFILE)
		@rm -f $@
		@cd $(dir $@) && $(LN_s) $(notdir $(JNILIBFILE)) $(notdir $@)

$(JNILIBLINK):	$(JNILIBFILE)
		@rm -f $@
		cd $(dir $@) && $(LN_s) $(notdir $(JNILIBFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(CLASSDIR):
		@-mkdir -p $(CLASSDIR)

-include	$(JNILIBDEP)

$(JNILIBFILE):	$(OBJDIR) $(JNILIBHEADERS) $(JNILIBDEP) $(JNILIBOBJFILES)
		@echo "-> generating library $@"
		-rm -f $@
		$(LIBBUILD) $(LIBBUILDFLAGS) $(LIBBUILD_o)$@ \
		$(JNILIBOBJFILES) \
		$(LINKCC_L)$(LIBDIR)/ $(LINKCC_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
		$(LINKCC_l)$(LPILIB)$(LINKLIBSUFFIX) \
		$(LINKCC_l)$(NLPILIB)$(LINKLIBSUFFIX) \
		$(LPSLDFLAGS) $(LDFLAGS)
ifneq ($(RANLIB),)
		$(RANLIB) $@
endif

$(JNILIBRARYLOADER):
		@echo "-> generating library loader $@"
		@echo "package de.zib.jscip.nativ.jni;"            > $(JNILIBRARYLOADER)
		@echo ""                                           >> $(JNILIBRARYLOADER)
		@echo "/** AUTOMATED GENERATED FILE */"            >> $(JNILIBRARYLOADER)
		@echo "public class JniScipLibraryLoader {"        >> $(JNILIBRARYLOADER)
		@echo ""                                           >> $(JNILIBRARYLOADER)
		@echo "    public static void loadLibrary() {"     >> $(JNILIBRARYLOADER)
		@echo "        System.loadLibrary(\"$(JNILIB)\");" >> $(JNILIBRARYLOADER)
		@echo "    }"                                      >> $(JNILIBRARYLOADER)
		@echo ""                                           >> $(JNILIBRARYLOADER)
		@echo "    private JniScipLibraryLoader() {"       >> $(JNILIBRARYLOADER)
		@echo ""                                           >> $(JNILIBRARYLOADER)
		@echo "    }"                                      >> $(JNILIBRARYLOADER)
		@echo "}"                                          >> $(JNILIBRARYLOADER)


$(JAVAPACKAGEJAR): 	$(JNICLASSES)
		@echo "-> generating jar $@"
		$(JAR) -cvf $@ $(JARCLASSES)

$(OBJDIR)/%.o:	$(SRCDIR)/%.c $(SRCDIR)/%.h
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(SRCDIR)/%.h:	$(CLASSDIR)/$(JNIPACKAGEDIR)/jni/%.class
		@echo "-> create C header $@"
		@echo "-> compiler command: $(JAVAH) -force -jni -o $@ -classpath $(CLASSDIR) $(JNIPACKAGENAME).jni.$(notdir $(basename $@))"
		$(JAVAH) -force -jni -o $@ -classpath $(CLASSDIR) $(JNIPACKAGENAME).jni.$(notdir $(basename $@))

$(CLASSDIR)/%.class:	$(JAVASRCDIR)/%.java $(CLASSDIR)
		@echo "-> compiling $@"
		$(JAVAC) $< -d $(CLASSDIR) -sourcepath $(JAVASRCDIR)

#---- EOF --------------------------------------------------------------------
