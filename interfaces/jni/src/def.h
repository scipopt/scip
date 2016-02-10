/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   def.c
 * @brief  common defines and data types used in all packages of JNISCIP
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __JNISCIP_DEF_H__
#define __JNISCIP_DEF_H__

#include "JniScipRetcode.h"

#define JNIPACKAGENAME(x) de_zib_jscip_nativ_jni_ ## x

#define JNISCIP(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScip_ ## x

#define JNISCIPBOUNDTYPE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipBoundType_ ## x
#define JNISCIPBOUNDCHGTYPE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipBoundchgType_ ## x

#define JNISCIPBRANCH(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipBranch_ ## x
#define JNISCIPBRANCHALLFULLSTRONG(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipBranchAllfullstrong_ ## x
#define JNISCIPBRANCHDIR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipBranchDir_ ## x
#define JNISCIPBRANCHFULLSTRONG(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipBranchFullstrong_ ## x
#define JNISCIPBRANCHINFERENCE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipBranchInference_ ## x
#define JNISCIPBRANCHLEASTINF(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipBranchLeastinf_ ## x
#define JNISCIPBRANCHMOSTINF(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipBranchMostinf_ ## x
#define JNISCIPBRANCHPSCOST(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipBranchPscost_ ## x
#define JNISCIPBRANCHRANDOM(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipBranchRandom_ ## x
#define JNISCIPBRANCHRELPSCOST(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipBranchRelpscost_ ## x

#define JNISCIPCONFLICT(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConflict_ ## x

#define JNISCIPCONS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipCons_ ## x
#define JNISCIPCONSABSPOWER(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsAbspower_ ## x
#define JNISCIPCONSAND(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsAnd_ ## x
#define JNISCIPCONSBIVARIATE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsBivariate_ ## x
#define JNISCIPCONSBOUNDDISJUNCTION(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsBounddisjunction_ ## x
#define JNISCIPCONSCONJUNCTION(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsConjunction_ ## x
#define JNISCIPCONSCOUNTSOLS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsCountsols_ ## x
#define JNISCIPCONSCUMULATIVE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsCumulative_ ## x
#define JNISCIPCONSDISJUNCTION(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsDisjunction_ ## x
#define JNISCIPCONSINDICATOR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsIndicator_ ## x
#define JNISCIPCONSINTEGRAL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsIntegral_ ## x
#define JNISCIPCONSKNAPSACK(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsKnapsack_ ## x
#define JNISCIPCONSLINEAR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsLinear_ ## x
#define JNISCIPCONSLINKING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsLinking_ ## x
#define JNISCIPCONSLOGICOR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsLogicor_ ## x
#define JNISCIPCONSNONLINEAR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsNonlinear_ ## x
#define JNISCIPCONSOR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsOr_ ## x
#define JNISCIPCONSORBITOPE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsOrbitope_ ## x
#define JNISCIPCONSPSEUDOBOOLEAN(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsPseudoboolean_ ## x
#define JNISCIPCONSQUADRATIC(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsQuadratic_ ## x
#define JNISCIPCONSSETPPC(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsSetppc_ ## x
#define JNISCIPCONSSOC(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsSoc_ ## x
#define JNISCIPCONSSOS1(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsSos1_ ## x
#define JNISCIPCONSSOS2(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsSos2_ ## x
#define JNISCIPCONSSUPERINDICATOR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsSuperindicator_ ## x
#define JNISCIPCONSVARBOUND(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsVarbound_ ## x
#define JNISCIPCONSXOR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipConsXor_ ## x

#define JNISCIPCUTPOOL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipCutpool_ ## x

#define JNISCIPDIALOG(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipDialog_ ## x
#define JNISCIPDIALOGDEFAULT(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipDialogDefault_ ## x

#define JNISCIPDISP(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipDisp_ ## x
#define JNISCIPDISPDEFAULT(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipDispDefault_ ## x

#define JNISCIPERRORCODE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipErrorCode_ ## x

#define JNISCIPEVENT(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipEvent_ ## x
#define JNISCIPEVENTTYPES(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipEventTypes_ ## x
#define JNISCIPEVENTHDLR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipEventhdlr_ ## x

#define JNISCIPEXCEPTION(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipException_ ## x

#define JNISCIPEXPR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipExpr_ ## x

#define JNISCIPHEUR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeur_ ## x
#define JNISCIPHEURACTCONSDIVING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurActconsdiving_ ## x
#define JNISCIPHEURCLIQUE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurClique_ ## x
#define JNISCIPHEURCOEFDIVING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurCoefdiving_ ## x
#define JNISCIPHEURCROSSOVER(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurCrossover_ ## x
#define JNISCIPHEURDINS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurDins_ ## x
#define JNISCIPHEURFEASPUMP(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurFeaspump_ ## x
#define JNISCIPHEURFIXANDINFER(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurFixandinfer_ ## x
#define JNISCIPHEURFRACDIVING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurFracdiving_ ## x
#define JNISCIPHEURGUIDEDDIVING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurGuideddiving_ ## x
#define JNISCIPHEURINTDIVING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurIntdiving_ ## x
#define JNISCIPHEURINTSHIFTING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurIntshifting_ ## x
#define JNISCIPHEURLINESEARCHDIVING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurLinesearchdiving_ ## x
#define JNISCIPHEURLOCALBRANCHING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurLocalbranching_ ## x
#define JNISCIPHEURMUTATION(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurMutation_ ## x
#define JNISCIPHEURNLPDIVING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurNlpdiving_ ## x
#define JNISCIPHEUROBJPSCOSTDIVING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurObjpscostdiving_ ## x
#define JNISCIPHEUROCTANE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurOctane_ ## x
#define JNISCIPHEURONEOPT(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurOneopt_ ## x
#define JNISCIPHEURPSCOSTDIVING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurPscostdiving_ ## x
#define JNISCIPHEURRENS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurRens_ ## x
#define JNISCIPHEURRINS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurRins_ ## x
#define JNISCIPHEURROOTSOLDIVING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurRootsoldiving_ ## x
#define JNISCIPHEURROUNDING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurRounding_ ## x
#define JNISCIPHEURSHIFTANDPROPAGATE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurShiftandpropagate_ ## x
#define JNISCIPHEURSHIFTING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurShifting_ ## x
#define JNISCIPHEURSIMPLEROUNDING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurSimplerounding_ ## x
#define JNISCIPHEURSUBNLP(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurSubnlp_ ## x
#define JNISCIPHEURTRIVIAL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurTrivial_ ## x
#define JNISCIPHEURTRYSOL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurTrysol_ ## x
#define JNISCIPHEURTWOOPT(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurTwoopt_ ## x
#define JNISCIPHEURUNDERCOVER(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurUndercover_ ## x
#define JNISCIPHEURVBOUNDS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurVbounds_ ## x
#define JNISCIPHEURVECLENDIVING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurVeclendiving_ ## x
#define JNISCIPHEURZEROOBJ(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurZeroobj_ ## x
#define JNISCIPHEURZIROUNDING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipHeurZirounding_ ## x

#define JNISCIPIMPLICS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipImplics_ ## x

#define JNISCIPLP(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipLp_ ## x

#define JNISCIPMESSAGE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipMessage_ ## x
#define JNISCIPMESSAGEDEFAULT(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipMessageDefault_ ## x

#define JNISCIPMISC(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipMisc_ ## x

#define JNISCIPNLP(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipNlp_ ## x
#define JNISCIPNLPIIPOPT(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipNlpiIpopt_ ## x

#define JNISCIPNODESEL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipNodesel_ ## x
#define JNISCIPNODESELBFS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipNodeselBfs_ ## x
#define JNISCIPNODESELDFS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipNodeselDfs_ ## x
#define JNISCIPNODESELESTIMATE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipNodeselEstimate_ ## x
#define JNISCIPNODESELHYBRIDESTIM(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipNodeselHybridestim_ ## x
#define JNISCIPNODESELRESTARTDFS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipNodeselRestartdfs_ ## x

#define JNISCIPOBJSENSE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipObjSense_ ## x

#define JNISCIPPARAMCHAR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipParamChar_ ## x
#define JNISCIPPARAMEMPHASIS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipParamEmphasis_ ## x
#define JNISCIPPARAMINT(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipParamInt_ ## x
#define JNISCIPPARAMLONG(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipParamLong_ ## x
#define JNISCIPPARAMNAME(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipParamName_ ## x
#define JNISCIPPARAMREAL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipParamReal_ ## x
#define JNISCIPPARAMSETTING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipParamSetting_ ## x
#define JNISCIPPARAMSTRING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipParamString_ ## x
#define JNISCIPPARAMTYPE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipParamType_ ## x
#define JNISCIPPARAMSET(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipParamset_ ## x

#define JNISCIPPRESOL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPresol_ ## x
#define JNISCIPPRESOLBOUNDSHIFT(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPresolBoundshift_ ## x
#define JNISCIPPRESOLCOMPONENTS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPresolComponents_ ## x
#define JNISCIPPRESOLCONVERTINTTOBIN(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPresolConvertinttobin_ ## x
#define JNISCIPPRESOLDOMCOL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPresolDomcol_ ## x
#define JNISCIPPRESOLDUALFIX(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPresolDualfix_ ## x
#define JNISCIPPRESOLGATEEXTRACTION(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPresolGateextraction_ ## x
#define JNISCIPPRESOLIMPLICS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPresolImplics_ ## x
#define JNISCIPPRESOLINTTOBINARY(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPresolInttobinary_ ## x
#define JNISCIPPRESOLTRIVIAL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPresolTrivial_ ## x

#define JNISCIPPRICER(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPricer_ ## x

#define JNISCIPPROP(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipProp_ ## x
#define JNISCIPPROPDUALFIX(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPropDualfix_ ## x
#define JNISCIPPROPGENVBOUNDS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPropGenvbounds_ ## x
#define JNISCIPPROPOBBT(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPropObbt_ ## x
#define JNISCIPPROPPROBING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPropProbing_ ## x
#define JNISCIPPROPPSEUDOOBJ(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPropPseudoobj_ ## x
#define JNISCIPPROPREDCOST(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPropRedcost_ ## x
#define JNISCIPPROPROOTREDCOST(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPropRootredcost_ ## x
#define JNISCIPPROPVBOUNDS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipPropVbounds_ ## x

#define JNISCIPREADER(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReader_ ## x
#define JNISCIPREADERBND(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderBnd_ ## x
#define JNISCIPREADERCCG(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderCcg_ ## x
#define JNISCIPREADERCIP(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderCip_ ## x
#define JNISCIPREADERCNF(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderCnf_ ## x
#define JNISCIPREADERFIX(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderFix_ ## x
#define JNISCIPREADERFZN(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderFzn_ ## x
#define JNISCIPREADERGMS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderGms_ ## x
#define JNISCIPREADERLP(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderLp_ ## x
#define JNISCIPREADERMPS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderMps_ ## x
#define JNISCIPREADEROPB(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderOpb_ ## x
#define JNISCIPREADEROSIL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderOsil_ ## x
#define JNISCIPREADERPIP(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderPip_ ## x
#define JNISCIPREADERPPM(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderPpm_ ## x
#define JNISCIPREADERRLP(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderRlp_ ## x
#define JNISCIPREADERSOL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderSol_ ## x
#define JNISCIPREADERWBO(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderWbo_ ## x
#define JNISCIPREADERZPL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipReaderZpl_ ## x

#define JNISCIPRELAX(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipRelax_ ## x

#define JNISCIPRESULT(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipResult_ ## x

#define JNISCIPRETCODE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipRetcode_ ## x

#define JNISCIPSEPA(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepa_ ## x
#define JNISCIPSEPACGMIP(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaCgmip_ ## x
#define JNISCIPSEPACLIQUE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaClique_ ## x
#define JNISCIPSEPACLOSECUTS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaClosecuts_ ## x
#define JNISCIPSEPACMIR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaCmir_ ## x
#define JNISCIPSEPAFLOWCOVER(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaFlowcover_ ## x
#define JNISCIPSEPAGOMORY(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaGomory_ ## x
#define JNISCIPSEPAIMPLIEDBOUNDS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaImpliedbounds_ ## x
#define JNISCIPSEPAINTOBJ(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaIntobj_ ## x
#define JNISCIPSEPAMCF(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaMcf_ ## x
#define JNISCIPSEPAODDCYCLE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaOddcycle_ ## x
#define JNISCIPSEPARAPIDLEARNING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaRapidlearning_ ## x
#define JNISCIPSEPASTRONGCG(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaStrongcg_ ## x
#define JNISCIPSEPAZEROHALF(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSepaZerohalf_ ## x

#define JNISCIPSOL(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSol_ ## x
#define JNISCIPSOLORIGIN(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipSolOrigin_ ## x

#define JNISCIPSTAGE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipStage_ ## x

#define JNISCIPSTATUS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipStatus_ ## x

#define JNISCIPTCLIQUECOLORING(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipTcliqueColoring_ ## x
#define JNISCIPTCLIQUEDEF(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipTcliqueDef_ ## x

#define JNISCIPTREE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipTree_ ## x

#define JNISCIPVAR(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipVar_ ## x
#define JNISCIPVARSTATUS(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipVarStatus_ ## x
#define JNISCIPVARTYPE(x) JNICALL Java_de_zib_jscip_nativ_jni_JniScipVarType_ ## x

#define JNISCIP_ENUM jint

#define JNISCIP_CALL(x)                                                 \
   do                                                                   \
   {                                                                    \
      SCIP_RETCODE _restat_;                                            \
                                                                        \
      /* check if retcodes in SCIP and JNI SCIP match */                \
      assert(SCIP_OKAY == JNIPACKAGENAME(JniScipRetcode_SCIP_OKAY));         \
      assert(SCIP_ERROR == JNIPACKAGENAME(JniScipRetcode_SCIP_ERROR));       \
      assert(SCIP_NOMEMORY == JNIPACKAGENAME(JniScipRetcode_SCIP_NOMEMORY)); \
      assert(SCIP_READERROR == JNIPACKAGENAME(JniScipRetcode_SCIP_READERROR)); \
      assert(SCIP_WRITEERROR == JNIPACKAGENAME(JniScipRetcode_SCIP_WRITEERROR)); \
      assert(SCIP_NOFILE == JNIPACKAGENAME(JniScipRetcode_SCIP_NOFILE));     \
      assert(SCIP_FILECREATEERROR == JNIPACKAGENAME(JniScipRetcode_SCIP_FILECREATEERROR)); \
      assert(SCIP_LPERROR == JNIPACKAGENAME(JniScipRetcode_SCIP_LPERROR));   \
      assert(SCIP_NOPROBLEM == JNIPACKAGENAME(JniScipRetcode_SCIP_NOPROBLEM)); \
      assert(SCIP_INVALIDCALL == JNIPACKAGENAME(JniScipRetcode_SCIP_INVALIDCALL)); \
      assert(SCIP_INVALIDDATA == JNIPACKAGENAME(JniScipRetcode_SCIP_INVALIDDATA)); \
      assert(SCIP_INVALIDRESULT ==  JNIPACKAGENAME(JniScipRetcode_SCIP_INVALIDRESULT)); \
      assert(SCIP_PLUGINNOTFOUND == JNIPACKAGENAME(JniScipRetcode_SCIP_PLUGINNOTFOUND)); \
      assert(SCIP_PARAMETERUNKNOWN == JNIPACKAGENAME(JniScipRetcode_SCIP_PARAMETERUNKNOWN)); \
      assert(SCIP_PARAMETERWRONGTYPE == JNIPACKAGENAME(JniScipRetcode_SCIP_PARAMETERWRONGTYPE)); \
      assert(SCIP_PARAMETERWRONGVAL == JNIPACKAGENAME(JniScipRetcode_SCIP_PARAMETERWRONGVAL)); \
      assert(SCIP_KEYALREADYEXISTING == JNIPACKAGENAME(JniScipRetcode_SCIP_KEYALREADYEXISTING)); \
      assert(SCIP_MAXDEPTHLEVEL == JNIPACKAGENAME(JniScipRetcode_SCIP_MAXDEPTHLEVEL)); \
      assert(SCIP_BRANCHERROR == JNIPACKAGENAME(JniScipRetcode_SCIP_BRANCHERROR)); \
                                                                        \
      if( (_restat_ = (x)) != SCIP_OKAY )                               \
      {                                                                 \
         jclass errorclass;                                             \
         jobject error;                                                 \
         jmethodID methodID;                                            \
         char msg[124];                                                 \
                                                                        \
         (void)SCIPsnprintf(msg, 124, "Error <%d> in function call", _restat_); \
                                                                        \
         /* check if retcodes is on we know JNI SCIP match */           \
         assert(_restat_ == SCIP_ERROR                                  \
            || _restat_ == SCIP_NOMEMORY                                \
            || _restat_ == SCIP_READERROR                               \
            || _restat_ == SCIP_WRITEERROR                              \
            || _restat_ == SCIP_NOFILE                                  \
            || _restat_ == SCIP_FILECREATEERROR                         \
            || _restat_ == SCIP_LPERROR                                 \
            || _restat_ == SCIP_NOPROBLEM                               \
            || _restat_ == SCIP_INVALIDCALL                             \
            || _restat_ == SCIP_INVALIDDATA                             \
            || _restat_ == SCIP_INVALIDRESULT                           \
            || _restat_ == SCIP_PLUGINNOTFOUND                          \
            || _restat_ == SCIP_PARAMETERUNKNOWN                        \
            || _restat_ == SCIP_PARAMETERWRONGTYPE                      \
            || _restat_ == SCIP_PARAMETERWRONGVAL                       \
            || _restat_ == SCIP_KEYALREADYEXISTING                      \
            || _restat_ == SCIP_MAXDEPTHLEVEL);                         \
                                                                        \
         SCIPerrorMessage("%s\n", msg);                                 \
                                                                        \
         errorclass = (*env)->FindClass(env, "de/zib/jscip/nativ/NativeScipException"); \
         methodID = (*env)->GetMethodID(env, errorclass, "<init>", "(I)V"); \
         error = (*env)->NewObject(env, errorclass, methodID, _restat_); \
         (*env)->Throw(env, error);                                     \
      }                                                                 \
   }                                                                    \
   while( FALSE )

#endif
