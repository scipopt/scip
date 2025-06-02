/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   certificate.cpp
 * @brief  methods for certificate output
 * @author Ambros Gleixner
 * @author Daniel Steffy
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <map>

#include "lpiexact/lpiexact.h"
#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/lp.h"
#include "scip/lpexact.h"
#include "scip/misc.h"
#include "scip/pub_cons.h"
#include "scip/pub_lpexact.h"
#include "scip/pub_misc.h"
#include "scip/pub_message.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/prob.h"
#include "scip/cuts.h"
#include "scip/cons_exactlinear.h"
#include "scip/scip_certificate.h"
#include "scip/scip_exact.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/struct_certificate.h"
#include "scip/struct_lpexact.h"
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/struct_var.h"
#include "scip/var.h"
#include "scip/certificate.h"

#define SCIP_HASHSIZE_CERTIFICATE    500 /**< size of hash map for certificate -> nodesdata mapping used for certificate output */
#define SCIP_MB_TO_CHAR_RATE   1048576.0 /**< conversion rate from MB to characters */

/** updates file size and returns whether maximum file size has been reached */
static
SCIP_Bool checkAndUpdateFilesize(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Real             nchars              /**< number of characters printed */
   )
{
   if( certificate->filesize < certificate->maxfilesize )
      certificate->filesize += nchars/(SCIP_MB_TO_CHAR_RATE);
   if( certificate->filesize < certificate->maxfilesize )
      return TRUE;
   return FALSE;
}

/** checks whether node is a left node or not */
static
SCIP_Bool certificateIsLeftNode(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_NODE*            node                /**< node from branch and bound tree */
   )
{
   SCIP_CERTNODEDATA* nodedata;
   SCIP_CERTNODEDATA* nodedataparent;

   assert(certificate != NULL);
   assert(SCIPcertificateIsEnabled(certificate));

   assert(node != NULL);
   assert(SCIPnodeGetType(node) != SCIP_NODETYPE_PROBINGNODE);

   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);

   if( SCIPnodeGetParent(node) == NULL )
      return FALSE;

   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedataparent = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, SCIPnodeGetParent(node));

   assert(nodedata->assumptionindex_self != -1);
   if( nodedataparent->assumptionindex_left == nodedata->assumptionindex_self )
      return TRUE;
   else
      return FALSE;
}

/** prints variable bound assumption into certificate
 *
 *  @return index of this bound in the certificate file
 */
static
SCIP_Longint printBoundAssumption(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_VAR*             var,                /**< variable to print assumption for */
   SCIP_RATIONAL*        boundval,           /**< value of the bound */
   SCIP_BOUNDTYPE        boundtype           /**< is it the upper bound? */
   )
{
   assert(certificate != NULL);
   assert(SCIPcertificateIsEnabled(certificate));

#ifndef NDEBUG
   certificate->lastinfo->isbound = TRUE;
   certificate->lastinfo->boundtype = boundtype;
   certificate->lastinfo->varindex = SCIPvarGetCertificateIndex(var);
   certificate->lastinfo->isglobal = FALSE;
   certificate->lastinfo->certificateindex = certificate->indexcounter;
   SCIPrationalSetRational(certificate->lastinfo->boundval, boundval);
#endif

   /** @todo it could be better to separate the printing from insertion of variable bound */
   SCIPcertificatePrintProofMessage(certificate, "A%lld %c ", certificate->indexcounter, (boundtype == SCIP_BOUNDTYPE_LOWER) ? 'G' : 'L');

   SCIP_CALL( SCIPcertificatePrintProofRational(certificate, boundval) );
   SCIPcertificatePrintProofMessage(certificate, " 1 %d 1 { asm } -1\n", SCIPvarGetCertificateIndex(var));
   certificate->indexcounter++;

   return certificate->indexcounter - 1;
}

/** free nodedata of corresponding node */
static
SCIP_RETCODE certificateFreeNodeData(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_NODE*            node                /**< focus node */
   )
{
   SCIP_CERTNODEDATA* nodedata;

   assert(node != NULL);
   assert(certificate != NULL);
   assert(SCIPcertificateIsEnabled(certificate));

   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);
   SCIPrationalFreeBlock(certificate->blkmem, &nodedata->derbound_left);
   SCIPrationalFreeBlock(certificate->blkmem, &nodedata->derbound_right);
   SCIPrationalFreeBlock(certificate->blkmem, &nodedata->derbound_self);
   BMSfreeBlockMemory(certificate->blkmem, &nodedata);
   SCIP_CALL( SCIPhashmapRemove(certificate->nodedatahash, node) );

   return SCIP_OKAY;
}

/** prints dual bound to proof section and increments indexcounter */
static
SCIP_RETCODE certificatePrintDualbound(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   const char*           linename,           /**< name of the unsplitting line */
   SCIP_RATIONAL*        lowerbound,         /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   int                   len,                /**< number of dual multipiers */
   SCIP_Longint*         ind,                /**< index array */
   SCIP_RATIONAL**       val                 /**< array of dual multipliers */
   )
{
   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   certificate->indexcounter++;
   certificate->lastinfo->isbound = FALSE;

   if( linename == NULL )
   {
      SCIPcertificatePrintProofMessage(certificate, "DualBound_%d ", certificate->indexcounter - 1);
   }
   else
   {
      SCIPcertificatePrintProofMessage(certificate, "%s ", linename);
   }

   if( SCIPrationalIsInfinity(lowerbound) )
   {
      SCIPcertificatePrintProofMessage(certificate, "G 1 0");
   }
   else
   {
      SCIPcertificatePrintProofMessage(certificate, "G ");
      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, lowerbound) );
      SCIPcertificatePrintProofMessage(certificate, " ");
      SCIPcertificatePrintProofMessage(certificate, "OBJ");
   }

   if( val == NULL )
   {
      SCIPcertificatePrintProofMessage(certificate, " { lin ... } -1\n");
   }
   else
   {
      int i;

      SCIPcertificatePrintProofMessage(certificate, " { lin %d", len);

      for( i = 0; i < len; i++ )
      {
         /** @todo perform line breaking before exceeding maximum line length */
         assert(!SCIPrationalIsAbsInfinity(val[i]));
         SCIPcertificatePrintProofMessage(certificate, " %d ", ind[i]);
         SCIP_CALL( SCIPcertificatePrintProofRational(certificate, val[i]) );
      }
      SCIPcertificatePrintProofMessage(certificate, " } -1\n");
   }

   /* print rounding derivation */
   if( !SCIPrationalIsNegInfinity(lowerbound) && certificate->objintegral && !SCIPrationalIsIntegral(lowerbound) )
   {
      certificate->indexcounter++;
      certificate->lastinfo->isbound = FALSE;

      SCIPcertificatePrintProofMessage(certificate, "R%d G ", certificate->indexcounter - 1);
      SCIPrationalRoundInteger(lowerbound, lowerbound, SCIP_R_ROUND_UPWARDS);

      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, lowerbound) );

      SCIPcertificatePrintProofMessage(certificate, " ");
      SCIPcertificatePrintProofMessage(certificate, "OBJ");
      SCIPcertificatePrintProofMessage(certificate, " { rnd 1 %d 1 } -1\n", certificate->indexcounter - 2);
   }

   return SCIP_OKAY;
}

/** prints the best solution found */
static
SCIP_RETCODE certificatePrintSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_SOL*             sol                 /**< solution to be printed */
   )
{
   SCIP_VAR** vars;
   SCIP_RATIONAL** vals;
   int nvars;
   int nnonz;
   int i;

   /* check if certificate output should be created */
   if( certificate->origfile == NULL )
      return SCIP_OKAY;

   assert(scip != NULL);

   if( sol == NULL )
   {
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "SOL 0\n");
      return SCIP_OKAY;
   }
   else if( !SCIPsolIsExact(sol) )
   {
      SCIP_CALL( SCIPmakeSolExact(scip, sol) );
   }

   /* get variables and number of the transformed problem */
   if( isorigfile )
   {
      SCIP_CALL( SCIPgetOrigVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   }
   else
   {
      SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   }

   SCIP_CALL( SCIPrationalCreateBufferArray(SCIPbuffer(scip), &vals, nvars) );

   /* get number of non-zero coefficient in the solution */
   nnonz = 0;
   for( i = 0; i < nvars; i++)
   {
      SCIPsolGetValExact(vals[i], sol, scip->set, scip->stat, vars[i]);
      if( !SCIPrationalIsZero(vals[i]) )
         nnonz++;
   }

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "SOL 1\nbest %d", nnonz);

   for( i = 0; i < nvars; i++ )
   {
      if( !SCIPrationalIsZero(vals[i]) )
      {
         /* print the solution into certificate */
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, " %d ", SCIPvarGetCertificateIndex(vars[i]));
         SCIP_CALL( SCIPcertificatePrintProblemRational(certificate, isorigfile, vals[i]) );
      }
   }
   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "\n");

   SCIPrationalFreeBufferArray(SCIPbuffer(scip), &vals, nvars);

   return SCIP_OKAY;
}

/** updates the current derived bound of the node with newbound, if newbound is better */
SCIP_RETCODE SCIPcertificateUpdateBoundData(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_NODE*            node,               /**< node data structure */
   SCIP_Longint          fileindex,          /**< index of new bound's proof */
   SCIP_RATIONAL*        newbound            /**< value of new bound */
   )
{
   SCIP_CERTNODEDATA* nodedata;

   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);

   /* do nothing if newbound is not better than the current bound */
   if( SCIPrationalIsLT(newbound, nodedata->derbound_self) )
      return SCIP_OKAY;

   nodedata->inheritedbound = FALSE;
   nodedata->derindex_self = fileindex;
   SCIPrationalSetRational(nodedata->derbound_self, newbound);

   return SCIP_OKAY;
}

/** creates certificate data structure */
SCIP_RETCODE SCIPcertificateCreate(
   SCIP_CERTIFICATE**    certificate,        /**< pointer to store the certificate information */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   SCIP_ALLOC( BMSallocMemory(certificate) );

   (*certificate)->messagehdlr = messagehdlr;
   (*certificate)->lastinfo = NULL;
   (*certificate)->blkmem = NULL;
   (*certificate)->indexcounter = 0;
   (*certificate)->indexcounter_ori = 0;
   (*certificate)->conscounter = 0;
   (*certificate)->origfile = NULL;
   (*certificate)->transfile = NULL;
   (*certificate)->origfilename = NULL;
   (*certificate)->derivationfile = NULL;
   (*certificate)->derivationfilename = NULL;
   (*certificate)->filesize = 0.0;
   (*certificate)->maxfilesize = SCIP_REAL_MAX;
   (*certificate)->rowdatahash = NULL;
   (*certificate)->naggrinfos = 0;
   (*certificate)->nmirinfos = 0;
   (*certificate)->aggrinfosize = 0;
   (*certificate)->mirinfosize = 0;
   (*certificate)->nodedatahash = NULL;
   (*certificate)->rootbound = NULL;
   (*certificate)->finalbound = NULL;
   (*certificate)->derindex_root = -1;
   (*certificate)->rootinfeas = FALSE;
   (*certificate)->objintegral = FALSE;
   (*certificate)->workingmirinfo = FALSE;
   (*certificate)->workingaggrinfo = FALSE;
   (*certificate)->vals = NULL;
   (*certificate)->valssize = 0;
   (*certificate)->aggrinfo = NULL;
   (*certificate)->mirinfo = NULL;
   (*certificate)->transfile_initialized = FALSE;

   return SCIP_OKAY;
}

/** frees certificate data structure */
void SCIPcertificateFree(
   SCIP_CERTIFICATE**    certificate         /**< pointer to store the certificate information */
   )
{
   assert(certificate != NULL);
   assert(*certificate != NULL);
   assert((*certificate)->origfile == NULL);
   assert((*certificate)->transfile == NULL);
   assert((*certificate)->derivationfile == NULL);

   BMSfreeMemory(certificate);
}

/* @todo replace scip pointer by set->scip */
/** initializes certificate information and creates files for certificate output */
SCIP_RETCODE SCIPcertificateInit(
   SCIP*                 scip,               /**< scip data structure */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   int filenamelen;
   int bufferlen;
   int nvars;
   int nintvars;
   int nbinvars;
   int nboundconss;
   int ncertcons;
   int ntransvars;
   int j;
   char* name = NULL;
   SCIP_VAR** vars;
   SCIP_VAR** transvars;
   SCIP_CONS** conss;
   SCIP_RATIONAL* lb;
   SCIP_RATIONAL* ub;

   assert(certificate != NULL);
   assert(set != NULL);
   assert(set->certificate_filename != NULL);
   assert(certificate->derivationfile == NULL);
   assert(certificate->nodedatahash == NULL);
   assert(certificate->rowdatahash == NULL);

   if( !(set->exact_enable) || (set->certificate_filename[0] == '-' && set->certificate_filename[1] == '\0') )
      return SCIP_OKAY;

   filenamelen = (int) strlen(set->certificate_filename);
   SCIP_CALL( SCIPsetAllocBufferArray(set, &name, filenamelen + 1) );
   BMScopyMemoryArray(name, set->certificate_filename, filenamelen);
   name[filenamelen] = '\0';

   SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      "storing certificate information in file <%s>\n", set->certificate_filename);

   certificate->transfile = SCIPfopen(set->certificate_filename, "wT");
   certificate->maxfilesize = set->certificate_maxfilesize;

   bufferlen = (int) strlen(name);
   SCIP_ALLOC( BMSallocMemoryArray(&certificate->derivationfilename, filenamelen+5) );
   SCIP_ALLOC( BMSallocMemoryArray(&certificate->origfilename, filenamelen+5) );
   BMScopyMemoryArray(certificate->derivationfilename, name, bufferlen);
   BMScopyMemoryArray(certificate->origfilename, name, bufferlen);
   certificate->derivationfilename[bufferlen] = '_';
   certificate->derivationfilename[bufferlen+1] = 'd';
   certificate->derivationfilename[bufferlen+2] = 'e';
   certificate->derivationfilename[bufferlen+3] = 'r';
   certificate->origfilename[bufferlen] = '_';
   certificate->origfilename[bufferlen+1] = 'o';
   certificate->origfilename[bufferlen+2] = 'r';
   certificate->origfilename[bufferlen+3] = 'i';
   certificate->derivationfilename[bufferlen+4] = '\0';
   certificate->origfilename[bufferlen+4] = '\0';
   certificate->derivationfile = SCIPfopen(certificate->derivationfilename, "wT");
   certificate->origfile = SCIPfopen(certificate->origfilename, "wT");

   if( certificate->transfile == NULL || certificate->origfile == NULL || certificate->derivationfile == NULL )
   {
      SCIPerrorMessage("error creating file <%s> and auxiliary certificate files\n", set->certificate_filename);
      SCIPprintSysError(set->certificate_filename);
      return SCIP_FILECREATEERROR;
   }

   /* initialisation of hashmaps and hashtables */
   SCIP_CALL( SCIPhashmapCreate(&certificate->nodedatahash, blkmem, SCIP_HASHSIZE_CERTIFICATE) );
   SCIP_CALL( SCIPhashmapCreate(&certificate->rowdatahash, blkmem, SCIP_HASHSIZE_CERTIFICATE) );
   SCIP_CALL( SCIPhashmapCreate(&certificate->aggrinfohash, blkmem, SCIP_HASHSIZE_CERTIFICATE) );
   SCIP_CALL( SCIPhashmapCreate(&certificate->mirinfohash, blkmem, SCIP_HASHSIZE_CERTIFICATE) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(certificate->aggrinfo), SCIP_HASHSIZE_CERTIFICATE) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(certificate->mirinfo), SCIP_HASHSIZE_CERTIFICATE) );
   certificate->aggrinfosize = SCIP_HASHSIZE_CERTIFICATE;
   certificate->mirinfosize = SCIP_HASHSIZE_CERTIFICATE;

   certificate->blkmem = blkmem;
   SCIPsetFreeBufferArray(set, &name);

   SCIP_CALL( SCIPgetOrigVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(scip, &transvars, &ntransvars, NULL, NULL, NULL, NULL) );
   nboundconss = 0;
   for( j = 0 ; j < nvars ; j++ )
   {
      lb = SCIPvarGetLbGlobalExact(vars[j]);
      ub = SCIPvarGetUbGlobalExact(vars[j]);
      if( !SCIPrationalIsAbsInfinity(lb) )
         nboundconss++;
      if( !SCIPrationalIsAbsInfinity(ub) )
         nboundconss++;
   }

   /* print the Version Header into certificate */
   SCIPcertificatePrintVersionHeader(certificate, TRUE);

   /* print the Variable Header into certificate */
   SCIPcertificatePrintVarHeader(certificate, TRUE, nvars);
   for( j = 0; j < nvars; j++ )
   {
      const char* varname;

      varname = SCIPvarGetName(vars[j]);
      SCIPvarSetCertificateIndex(vars[j], j);
      if( strstr(varname, " ") != NULL || strstr(varname, "\t") != NULL || strstr(varname, "\n") != NULL
         || strstr(varname, "\v") != NULL || strstr(varname, "\f") != NULL || strstr(varname, "\r") != NULL )
      {
         SCIPerrorMessage("Variable name <%s> cannot be printed to certificate file because it contains whitespace.\n",
            varname);
         return SCIP_ERROR;
      }

      SCIPcertificatePrintProblemMessage(certificate, TRUE, "%s\n", varname);
   }

   /* print the Integer Variable Header into certificate */
   SCIPcertificatePrintIntHeader(certificate, TRUE, nintvars + nbinvars);
   for( j = 0; j < nvars; j++ )
   {
      if( SCIPvarGetType(vars[j]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(vars[j]) == SCIP_VARTYPE_INTEGER )
      {
         SCIPcertificatePrintProblemMessage(certificate, TRUE, "%d \n", SCIPvarGetCertificateIndex(vars[j]));
      }
   }

   {
      SCIP_RATIONAL** objcoefs;
      SCIP_CALL( SCIPrationalCreateBufferArray(SCIPbuffer(scip), &objcoefs, nvars) );

      for( j = 0; j < nvars; j++)
         SCIPrationalSetRational(objcoefs[j], SCIPvarGetObjExact(vars[j]));

      /* print the objective function into certificate header */
      SCIP_CALL( SCIPcertificateSetAndPrintObjective(certificate, TRUE, blkmem, objcoefs, nvars) );

      SCIPrationalFreeBufferArray(SCIPbuffer(scip), &objcoefs, nvars);
   }

   conss = SCIPgetOrigConss(scip);
   ncertcons = 0;
   for( j = 0; j < SCIPgetNOrigConss(scip); j++ )
   {
      SCIP_CONS* cons;
      SCIP_CONSHDLR* conshdlr;

      cons = conss[j];
      conshdlr = SCIPconsGetHdlr(cons);

      if( strcmp(SCIPconshdlrGetName(conshdlr), "exactlinear") == 0 )
      {
         lb = SCIPgetLhsExactLinear(scip, cons);
         ub = SCIPgetRhsExactLinear(scip, cons);

         if( !SCIPrationalIsEQ(lb, ub) && !SCIPrationalIsAbsInfinity(lb) && !SCIPrationalIsAbsInfinity(ub) )
            ncertcons += 2;
         else
            ncertcons += 1;
      }
      else
      {
         SCIPerrorMessage("Cannot print certificate for non-exact constraints \n");
         SCIPABORT();
         return SCIP_ERROR;
      }
   }

   SCIPcertificatePrintConsHeader(certificate, TRUE, ncertcons, nboundconss);

   for( j = 0; j < SCIPgetNOrigConss(scip); j++ )
   {
      SCIP_CONS* cons;
      cons = conss[j];
      SCIP_CALL( SCIPcertifyConsOrigExactLinear(scip, SCIPconsGetHdlr(cons), cons) );
   }

   for( j = 0; j < nvars; j++ )
   {
      if( !SCIPrationalIsAbsInfinity(SCIPvarGetLbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, TRUE, NULL, vars[j], SCIPvarGetLbGlobalExact(vars[j]), FALSE) );
      }
      if( !SCIPrationalIsAbsInfinity(SCIPvarGetUbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, TRUE, NULL, vars[j], SCIPvarGetUbGlobalExact(vars[j]), TRUE) );
      }
   }

   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &certificate->rootbound) );
   SCIP_CALL( SCIPrationalCreateBlock(blkmem, &certificate->finalbound) );
   certificate->valssize = SCIPgetNVars(scip) + SCIPgetNConss(scip);
   SCIP_CALL( SCIPrationalCreateBlockArray(SCIPblkmem(scip), &(certificate->vals), certificate->valssize) );

   return SCIP_OKAY;
}

/** initializes certificate information and creates files for certificate output */
SCIP_RETCODE SCIPcertificateInitTransFile(
   SCIP*                 scip                /**< scip data structure */
   )
{
   int nvars;
   int nintvars;
   int nbinvars;
   int nboundconss;
   int ncertcons;
   int j;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_RATIONAL* lb;
   SCIP_RATIONAL* ub;
   SCIP_CERTIFICATE* certificate;
   BMS_BLKMEM* blkmem;
   SCIP_Bool cutoff;

   assert(scip != NULL);
   assert(scip->set->certificate_filename != NULL);

   certificate = SCIPgetCertificate(scip);
   blkmem = SCIPblkmem(scip);
   cutoff = FALSE;

   assert(certificate != NULL);

   if( certificate->transfile_initialized )
      return SCIP_OKAY;

   /* the transfile is constructed using the (exact) LP, so make sure this is constructed here */
   if( !SCIPisLPConstructed(scip) )
   {
      SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   }

   if( certificate->transfile == NULL || certificate->origfile == NULL || certificate->derivationfile == NULL )
   {
      SCIPerrorMessage("error creating file <%s> and auxiliary certificate files\n", scip->set->certificate_filename);
      SCIPprintSysError(scip->set->certificate_filename);
      return SCIP_FILECREATEERROR;
   }

   certificate->transfile_initialized = TRUE;
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   nboundconss = 0;
   for( j = 0 ; j < nvars ; j++ )
   {
      lb = SCIPvarGetLbGlobalExact(vars[j]);
      ub = SCIPvarGetUbGlobalExact(vars[j]);
      if( !SCIPrationalIsAbsInfinity(lb) )
         nboundconss++;
      if( !SCIPrationalIsAbsInfinity(ub) )
         nboundconss++;
   }

   /* print the Version Header into certificate */
   SCIPcertificatePrintVersionHeader(certificate, FALSE);

   /* print the Variable Header into certificate */
   SCIPcertificatePrintVarHeader(certificate, FALSE, nvars);
   for( j = 0; j < nvars; j++ )
   {
      const char* varname;

      varname = SCIPvarGetName(vars[j]);
      SCIPvarSetCertificateIndex(vars[j], j);
      if( strstr(varname, " ") != NULL || strstr(varname, "\t") != NULL || strstr(varname, "\n") != NULL
         || strstr(varname, "\v") != NULL || strstr(varname, "\f") != NULL || strstr(varname, "\r") != NULL )
      {
         SCIPerrorMessage("Variable name <%s> cannot be printed to certificate file because it contains whitespace.\n",
            varname);
         return SCIP_ERROR;
      }

      SCIPcertificatePrintProblemMessage(certificate, FALSE, "%s\n", varname);
   }

   /* print the Integer Variable Header into certificate */
   SCIPcertificatePrintIntHeader(certificate, FALSE, nintvars + nbinvars);
   for( j = 0; j < nvars; j++ )
   {
      if( SCIPvarGetType(vars[j]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(vars[j]) == SCIP_VARTYPE_INTEGER )
      {
         SCIPcertificatePrintProblemMessage(certificate, FALSE, "%d \n", SCIPvarGetCertificateIndex(vars[j]));
      }
   }

   {
      SCIP_RATIONAL** objcoefs;
      SCIP_CALL( SCIPrationalCreateBufferArray(SCIPbuffer(scip), &objcoefs, nvars) );

      for( j = 0; j < nvars; j++)
         SCIPrationalSetRational(objcoefs[j], SCIPvarGetObjExact(vars[j]));

      /* print the objective function into certificate header */
      SCIP_CALL( SCIPcertificateSetAndPrintObjective(certificate, FALSE, blkmem, objcoefs, nvars) );

      SCIPrationalFreeBufferArray(SCIPbuffer(scip), &objcoefs, nvars);
   }

   conss = SCIPgetConss(scip);
   ncertcons = 0;
   for( j = 0; j < SCIPgetNConss(scip); j++ )
   {
      SCIP_CONS* cons;
      SCIP_CONSHDLR* conshdlr;

      cons = conss[j];
      conshdlr = SCIPconsGetHdlr(cons);

      SCIPdebug(SCIPprintCons(scip, conss[j], NULL));

      if( strcmp(SCIPconshdlrGetName(conshdlr), "exactlinear") == 0 )
      {
         lb = SCIPgetLhsExactLinear(scip, cons);
         ub = SCIPgetRhsExactLinear(scip, cons);

         if( !SCIPrationalIsEQ(lb, ub) && !SCIPrationalIsAbsInfinity(lb) && !SCIPrationalIsAbsInfinity(ub) )
         {
            SCIPdebugMessage("constraint is a ranged constraint \n");
            ncertcons += 2;
         }
         else
         {
            SCIPdebugMessage("constraint only has one side \n");
            ncertcons += 1;
         }
      }
      else
      {
         SCIPerrorMessage("Cannot print certificate for non-exact constraints \n");
         SCIPABORT();
         return SCIP_ERROR;
      }
   }

   SCIPcertificatePrintConsHeader(certificate, FALSE, ncertcons, nboundconss);

   for( j = 0; j < SCIPgetNConss(scip); j++ )
   {
      SCIP_CONS* cons;
      cons = conss[j];
      SCIP_CALL( SCIPconsPrintCertificateExactLinear(scip, cons) );
   }

   for( j = 0; j < nvars; j++ )
   {
      if( !SCIPrationalIsAbsInfinity(SCIPvarGetLbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, FALSE, NULL, vars[j], SCIPvarGetLbGlobalExact(vars[j]), FALSE) );
         SCIPvarSetLbCertificateIndexGlobal(vars[j], certificate->indexcounter - 1);
         SCIPvarSetLbCertificateIndexLocal(vars[j], certificate->indexcounter - 1);
      }
      if( !SCIPrationalIsAbsInfinity(SCIPvarGetUbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, FALSE, NULL, vars[j], SCIPvarGetUbGlobalExact(vars[j]), TRUE) );
         SCIPvarSetUbCertificateIndexGlobal(vars[j], certificate->indexcounter - 1);
         SCIPvarSetUbCertificateIndexLocal(vars[j], certificate->indexcounter - 1);
      }
   }

   return SCIP_OKAY;
}

/** concatenates the certificate and the _der file and deletes the _der file */
static
void concatenateCertificate(
   SCIP_CERTIFICATE*     certificate         /**< the certificate pointer */
   )
{
   SCIP_FILE* derivationfile;
   char buffer[SCIP_MAXSTRLEN];
   size_t size;

   derivationfile = SCIPfopen(certificate->derivationfilename, "r");

   /* append the derivation file to the problem file */
   while( 0 != (size = SCIPfread(buffer, sizeof(char), SCIP_MAXSTRLEN, derivationfile)) )
      (void) SCIPfwrite(buffer, sizeof(char), size, certificate->transfile);

   SCIPfclose(derivationfile);

   /* delete the derivation file */
   (void) remove(certificate->derivationfilename);
}

/** closes the certificate output files */
SCIP_RETCODE SCIPcertificateExit(
   SCIP*                 scip                /**< scip data structure */
   )
{
   assert(scip != NULL);

   SCIP_CERTIFICATE* certificate = SCIPgetCertificate(scip);
   SCIP_MESSAGEHDLR* messagehdlr = SCIPgetMessagehdlr(scip);
   SCIP_SET* set = scip->set;

   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   assert(certificate != NULL);

   if( certificate->origfile != NULL )
   {
      SCIP_Bool printingaborted = checkAndUpdateFilesize(certificate, 0) ? FALSE : TRUE;

      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
         "closing certificate file (wrote approx. %.1f MB%s)\n", certificate->filesize,
         printingaborted ? ", aborted printing after reaching max. file size" : "");

      if( printingaborted )
      {
         (void) SCIPfprintf(certificate->origfile, "\n# ... aborted printing: max. file size reached.\n");
         (void) SCIPfprintf(certificate->transfile, "\n# ... aborted printing: max. file size reached.\n");
         (void) SCIPfprintf(certificate->derivationfile, "\n# ... aborted printing: max. file size reached.\n");
      }

      if( certificate->derivationfile != NULL )
      {
         SCIPfclose(certificate->derivationfile);
         certificate->derivationfile = NULL;
         concatenateCertificate(certificate);
         /* if the file is empty (e.g. because we detected infeasibility in presolving) we delete it */
         if( certificate->indexcounter == 0 )
         {
            SCIPdebugMessage("derivation file is empty; deleting it");
            (void) remove(set->certificate_filename);
         }
      }
      SCIPfclose(certificate->origfile);
      SCIPfclose(certificate->transfile);
      certificate->origfile = NULL;
      certificate->transfile = NULL;

      BMSfreeMemoryArray(&certificate->derivationfilename);
      BMSfreeMemoryArray(&certificate->origfilename);

      if( certificate->lastinfo != NULL )
      {
         SCIPrationalFreeBlock(certificate->blkmem, &certificate->lastinfo->boundval);
         BMSfreeBlockMemory(certificate->blkmem, &certificate->lastinfo);
      }

      if( certificate->rowdatahash)
         SCIPhashmapFree(&certificate->rowdatahash);

      if( certificate->nodedatahash )
      {
         assert(SCIPhashmapIsEmpty(certificate->nodedatahash));
         SCIPhashmapFree(&certificate->nodedatahash);
      }
      if( certificate->aggrinfohash )
      {
         SCIP_CALL( SCIPcertificateClearAggrinfo(scip) );
      }
      if( certificate->mirinfohash )
      {
         SCIP_CALL( SCIPcertificateClearMirinfo(scip) );
      }

      SCIPrationalFreeBlock(certificate->blkmem, &certificate->rootbound);
      SCIPrationalFreeBlock(certificate->blkmem, &certificate->finalbound);
      SCIPrationalFreeBlockArray(certificate->blkmem, &certificate->vals, certificate->valssize);
   }

   return SCIP_OKAY;
}

/** returns certificate data structure */
SCIP_CERTIFICATE* SCIPgetCertificate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);

   return scip->stat->certificate;
}

/** returns whether the certificate output is activated */
SCIP_Bool SCIPcertificateIsEnabled(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   if( certificate != NULL && certificate->transfile != NULL && certificate->origfile != NULL && certificate->derivationfile != NULL )
      return TRUE;
   return FALSE;
}

/** returns current certificate file size in MB */
SCIP_Real SCIPcertificateGetFilesize(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   if( !SCIPcertificateIsEnabled(certificate) )
      return 0.0;
   else
      return certificate->filesize;
}

/** returns current certificate index (return -1 if certificate not active) */
SCIP_Longint SCIPcertificateGetCurrentIndex(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   if( !SCIPcertificateIsEnabled(certificate) )
      return -1;
   else
      return certificate->indexcounter;
}

#ifndef NDEBUG
/** checks if information is consistent with printed certificate line */
SCIP_Bool SCIPcertificateEnsureLastBoundInfoConsistent(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_VAR*             var,                /**< variable that gets changed */
   SCIP_BOUNDTYPE        boundtype,          /**< lb or ub changed? */
   SCIP_Real             newbound,           /**< new bound */
   SCIP_Bool             needsglobal         /**< if the bound needs to be global */
   )
{
   SCIP_Bool consistent;

   if( !SCIPcertificateIsEnabled(certificate) )
      return TRUE;

   assert(certificate != NULL);

   consistent = certificate->lastinfo->isbound;
   consistent = consistent && certificate->lastinfo->varindex == SCIPvarGetCertificateIndex(var); /*lint !e1785*/
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      consistent = consistent && certificate->lastinfo->boundtype == SCIP_BOUNDTYPE_LOWER;  /*lint !e1785*/
      consistent = consistent && SCIPrationalRoundReal(certificate->lastinfo->boundval, SCIP_R_ROUND_DOWNWARDS) >= newbound; /*lint !e1785*/
   }
   else
   {
      consistent = consistent && certificate->lastinfo->boundtype == SCIP_BOUNDTYPE_UPPER; /*lint !e1785*/
      consistent = consistent && SCIPrationalRoundReal(certificate->lastinfo->boundval, SCIP_R_ROUND_UPWARDS) <= newbound; /*lint !e1785*/
   }
   consistent = consistent && (!needsglobal || certificate->lastinfo->isglobal); /*lint !e1785*/
   consistent = consistent && certificate->lastinfo->certificateindex == certificate->indexcounter - 1; /*lint !e1785*/

   return consistent;
}
#endif

/** sets the objective function used when printing dual bounds */
SCIP_RETCODE SCIPcertificateSetAndPrintObjective(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the line be printed to the origfile or the transfile */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_RATIONAL**       coefs,              /**< objective function coefficients */
   int                   nvars               /**< number of variables */
   )
{
   char obj[SCIP_MAXSTRLEN - 2];
   int nnonz;
   int i;

   assert(coefs != NULL);

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   /* create a hash table for the variable bounds (we work on the tranformed problem) */
   if( isorigfile )
   {
      /* create working memory for bound struct */
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &certificate->lastinfo) );
      SCIP_CALL( SCIPrationalCreateBlock(blkmem, &certificate->lastinfo->boundval) );
   }

   nnonz = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( !SCIPrationalIsZero(coefs[i]) )
         nnonz++;
   }

   (void) SCIPsnprintf(obj, SCIP_MAXSTRLEN - 2, "OBJ min\n %d ", nnonz);

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s ", obj);

   for( i = 0; i < nvars; i++ )
   {
      if( !SCIPrationalIsZero(coefs[i]) )
      {
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s%d ", (i > 0 ? " " : ""), i );
         SCIP_CALL( SCIPcertificatePrintProblemRational(certificate, isorigfile, coefs[i]) );
      }
   }

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "\n");

   return SCIP_OKAY;
}

/** prints the last part of the certificate header (RTP range/sol, ...) */
SCIP_RETCODE SCIPcertificatePrintResult(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             isorigfile,         /**< should the line be printed to the origfile or the transfile */
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   SCIP_RATIONAL* primalbound;
   SCIP_RATIONAL* dualbound;
   SCIP_SOL* bestsol;

   assert(scip != NULL);

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &primalbound) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &dualbound) );

   /* check status first: OPTIMAL / INFEAS / NOT OPTIMAL */
   if( SCIPisInRestart(scip) )
      return SCIP_ERROR;

   if( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL )
   {
      bestsol = SCIPgetBestSol(scip);

      if( !SCIPsolIsExact(bestsol) )
      {
         SCIP_CALL( SCIPmakeSolExact(scip, bestsol) );
      }

      if( isorigfile )
      {
         SCIP_CALL( SCIPretransformSolExact(scip, bestsol) );
      }

      if( isorigfile )
         SCIPgetPrimalboundExact(scip, primalbound);
      else
      {
         SCIPgetUpperboundExact(scip, primalbound);
      }

      assert(!SCIPrationalIsAbsInfinity(primalbound));

      /* for the orig file we only print the primal bound, since the derivation happens in the transformed problem */
      if( isorigfile )
      {
         SCIPrationalSetNegInfinity(dualbound);
         /* print RTP range (same when optimal solution found) */
         SCIP_CALL( SCIPcertificatePrintRtpRange(certificate, isorigfile, dualbound, primalbound) );
      }
      else
         SCIP_CALL( SCIPcertificatePrintRtpRange(certificate, isorigfile, primalbound, primalbound) );

      /* print optimal solution into certificate */
      SCIP_CALL( certificatePrintSol(scip, isorigfile, certificate, bestsol) );
   }
   else if( SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE )
   {
      SCIPcertificatePrintRtpInfeas(certificate, isorigfile);
      SCIP_CALL(certificatePrintSol(scip, isorigfile, certificate, NULL) );
   }
   else
   {
      /* two cases to distinguish: a primal bound has been found or not */
      if( SCIPisPrimalboundSol(scip) )
      {
         bestsol = SCIPgetBestSol(scip);
         if( isorigfile )
         {
            SCIP_CALL( SCIPretransformSolExact(scip, bestsol) );
         }

         if( isorigfile )
            SCIPgetPrimalboundExact(scip, primalbound);
         else
            SCIPgetUpperboundExact(scip, primalbound);
      }
      else
      {
         bestsol = NULL;
         SCIPrationalSetInfinity(primalbound);
      }

      if( isorigfile )
         SCIPrationalSetNegInfinity(dualbound);
      else
         SCIPrationalSetRational(dualbound, certificate->finalbound);

      SCIP_CALL( SCIPcertificatePrintRtpRange(certificate, isorigfile, dualbound, primalbound) );
      SCIP_CALL( certificatePrintSol(scip, isorigfile, certificate, bestsol) );
   }

   SCIPcertificatePrintDerHeader(certificate, isorigfile);

   SCIPrationalFreeBuffer(set->buffer, &dualbound);
   SCIPrationalFreeBuffer(set->buffer, &primalbound);

   return SCIP_OKAY;
}

/** prints the last part of the certificate header (RTP range/sol, ...) */
SCIP_RETCODE SCIPcertificateSaveFinalbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   assert(scip != NULL);

   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   SCIPgetLowerboundExact(scip, certificate->finalbound);

   return SCIP_OKAY;
}

/** prints a string to the problem section of the certificate file */
void SCIPcertificatePrintProblemMessage(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the line be printed to the origfile or the transfile */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[3 * SCIP_MAXSTRLEN];

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return;

   va_start(ap, formatstr);
   (void) vsnprintf(buffer, 3 * SCIP_MAXSTRLEN, formatstr, ap);

   if( checkAndUpdateFilesize(certificate, strlen(buffer)) )
   {
      if( isorigfile )
         (void) SCIPfprintf(certificate->origfile, "%s", buffer);
      else
         (void) SCIPfprintf(certificate->transfile, "%s", buffer);
   }
   va_end(ap);
}

/** prints a string to the proof section of the certificate file */
void SCIPcertificatePrintProofMessage(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[3 * SCIP_MAXSTRLEN];

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return;

   va_start(ap, formatstr);
   (void) vsnprintf(buffer, 3 * SCIP_MAXSTRLEN, formatstr, ap);

   if( checkAndUpdateFilesize(certificate, strlen(buffer)) )
      (void) SCIPfprintf(certificate->derivationfile, "%s", buffer);

   va_end(ap);
}


/** prints a rational number to the problem section of the certificate file */
SCIP_RETCODE SCIPcertificatePrintProblemRational(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the line be printed to the origfile or the transfile */
   SCIP_RATIONAL*        val                 /**< rational number to print */
   )
{
   int len = SCIPrationalStrLen(val) + 1;
   char* buffer = NULL;

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
     return SCIP_OKAY;

   SCIP_ALLOC( BMSallocMemoryArray(&buffer, len) );
   (void) SCIPrationalToString(val, buffer, len);
   if( checkAndUpdateFilesize(certificate, strlen(buffer)) )
   {
      if( isorigfile )
         (void) SCIPfputs(buffer, certificate->origfile);
      else
         (void) SCIPfputs(buffer, certificate->transfile);
   }
   BMSfreeMemoryArray(&buffer);

   return SCIP_OKAY;
}


/** prints a rational number to the proof section of the certificate file */
SCIP_RETCODE SCIPcertificatePrintProofRational(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_RATIONAL*        val                 /**< rational number to print */
   )
{
   int len = SCIPrationalStrLen(val) + 1;
   char* buffer = NULL;

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
     return SCIP_OKAY;;

   SCIP_ALLOC( BMSallocMemoryArray(&buffer, len) );
   (void) SCIPrationalToString(val, buffer, len);

   if( checkAndUpdateFilesize(certificate, strlen(buffer)) )
      (void) SCIPfputs(buffer, certificate->derivationfile);

   BMSfreeMemoryArray(&buffer);

   return SCIP_OKAY;
}

/** prints a comment to the problem section of the certificate file */
void SCIPcertificatePrintProblemComment(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the line be printed to the origfile or the transfile */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[3 * SCIP_MAXSTRLEN];

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return;

   (void) SCIPfprintf(certificate->origfile, "# ");

   va_start(ap, formatstr);
   (void) vsnprintf(buffer, 3 * SCIP_MAXSTRLEN, formatstr, ap);

   if( checkAndUpdateFilesize(certificate, 2 + strlen(buffer)) )
   {
      if( isorigfile )
         (void) SCIPfprintf(certificate->origfile, "%s", buffer);
      else
         (void) SCIPfprintf(certificate->transfile, "%s", buffer);
   }
   va_end(ap);
}

/** prints a comment to the proof section of the certificate file */
void SCIPcertificatePrintProofComment(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[3 * SCIP_MAXSTRLEN];

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return;

   (void) SCIPfprintf(certificate->derivationfile, "# ");

   va_start(ap, formatstr);
   (void) vsnprintf(buffer, 3 * SCIP_MAXSTRLEN, formatstr, ap);

   if( checkAndUpdateFilesize(certificate, 2 + strlen(buffer)) )
      (void) SCIPfprintf(certificate->derivationfile, "%s", buffer);

   va_end(ap);
}

/** prints version header */
void SCIPcertificatePrintVersionHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile          /**< should the line be printed to the origfile or the transfile */
   )
{
   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return;

   assert(certificate != NULL);

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "VER 1.0 \n");
}

/** prints variable section header */
void SCIPcertificatePrintVarHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the line be printed to the origfile or the transfile */
   int                   nvars               /**< number of variables */
   )
{
   assert(nvars >= 0);

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return;

   assert(certificate != NULL);

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "VAR %d \n", nvars);
}

/** prints integer section header */
void SCIPcertificatePrintIntHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the line be printed to the origfile or the transfile */
   int                   nints               /**< number of integer variables */
   )
{
   assert(nints >= 0);

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return;

   assert(certificate != NULL);

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "INT %d\n", nints);
}

/** prints constraint section header */
void SCIPcertificatePrintConsHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the line be printed to the origfile or the transfile */
   int                   nconss,             /**< number of all constraints */
   int                   nboundconss         /**< number of bound constraints */
   )
{
   assert(nconss >= 0);

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return;

   assert(certificate != NULL);

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "CON %d %d\n", nconss + nboundconss, nboundconss);
}

/** prints derivation section header */
void SCIPcertificatePrintDerHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile          /**< should the line be printed to the origfile or the transfile */
   )
{
   SCIP_Longint nders;

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return;

   assert(certificate != NULL);

   if( !isorigfile )
      nders = certificate->indexcounter - certificate->conscounter;
   else
      nders = 0;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "DER %d\n", nders);
}

/** prints constraint */
SCIP_RETCODE SCIPcertificatePrintCons(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the line be printed to the origfile or the transfile */
   const char*           consname,           /**< name of the constraint */
   const char            sense,              /**< sense of the constraint, i.e., G, L, or E */
   SCIP_RATIONAL*        side,               /**< left/right-hand side */
   int                   len,                /**< number of nonzeros */
   int*                  ind,                /**< index array */
   SCIP_RATIONAL**       val                 /**< coefficient array */
   )
{
   SCIP_Longint index;
   int i;

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   index = isorigfile ? certificate->indexcounter_ori : certificate->indexcounter;

   if( consname == NULL )
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "C%d %c ", index, sense);
   else
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s %c ", consname, sense);

   SCIP_CALL( SCIPcertificatePrintProblemRational(certificate, isorigfile, side) );

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, " %d", len);

   for( i = 0; i < len; i++ )
   {
      /** @todo perform line breaking before exceeding maximum line length */
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, " %d ", ind[i]);
      SCIP_CALL( SCIPcertificatePrintProblemRational(certificate, isorigfile, val[i]) );
   }
   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "\n");

   if( isorigfile )
      certificate->indexcounter_ori++;
   else
   {
      certificate->indexcounter++;
      certificate->lastinfo->isbound = FALSE;
   }

   if( !isorigfile )
   {
      certificate->conscounter++;
   }

   return SCIP_OKAY;
}

/** prints a line for an exact row to the certificate (without derivation)
 *
 *  @param alternativerhs is used instead of the real rhs of the row (infinity if real rhs should be used).
 *  This is necessary for integer cut where the rhs was rounded down from the original rhs
 */
static
SCIP_RETCODE certificatePrintRow(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate structure */
   SCIP_ROWEXACT*        rowexact,           /**< exact SCIP row */
   SCIP_Real             alternativerhs      /**< rhs to be used instead or rowexact->rhs (infinity to disable this) */
   )
{
   SCIP_ROW* row;
   SCIP_RATIONAL* rhs;
   int i;

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   assert(rowexact != NULL);

   row = SCIProwExactGetRow(rowexact);
   assert(SCIProwGetNNonz(row) == SCIProwExactGetNNonz(rowexact));

   SCIPcertificatePrintProofMessage(certificate, "L%d_%s %c ", certificate->indexcounter, row->name, 'L');

   /* if we rounded the rhs in cut tightening we need to first verify the cut without it */
   if( SCIPsetIsInfinity(set, alternativerhs) )
   {
      rhs = SCIProwExactGetRhs(rowexact);
      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, rhs) );
   }
   else
   {
      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &rhs) );
      SCIPrationalSetReal(rhs, alternativerhs);
      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, rhs) );
      SCIPrationalFreeBuffer(set->buffer, &rhs);
   }

   SCIPcertificatePrintProofMessage(certificate, " %d", SCIProwGetNNonz(row));

   for( i = 0; i < SCIProwGetNNonz(row); i++ )
   {
      SCIP_RATIONAL* val;
      int varindex;
      /** @todo perform line breaking before exceeding maximum line length */
      assert(rowexact->cols[i]->fpcol->var->index == rowexact->fprow->cols[i]->var->index);

      varindex = SCIPvarGetCertificateIndex(SCIPcolExactGetVar(SCIProwExactGetCols(rowexact)[i]));
      val = SCIProwExactGetVals(rowexact)[i];

      SCIPcertificatePrintProofMessage(certificate, " %d ", varindex);
      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, val) );
   }
   certificate->indexcounter++;
   certificate->lastinfo->isbound = FALSE;

   return SCIP_OKAY;
}

/** prints mir split for the specified aggrrow */
static
SCIP_RETCODE certificatePrintMirSplit(
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_ROW*             row                 /**< row that split should be printed for */
   )
{
   SCIP_MIRINFO* mirinfo;
   SCIP_RATIONAL* val;
   SCIP_VAR** vars;
   int i, j;
   SCIP_Real slackrhs;
   std::map<int, SCIP_Real> coefs;

   assert(SCIPcertificateIsEnabled(certificate));
   assert(SCIPhashmapExists(certificate->mirinfohash, (void*) row));

   mirinfo = (SCIP_MIRINFO*) SCIPhashmapGetImage(certificate->mirinfohash, (void*) row);

   vars = SCIPprobGetVars(prob);
   slackrhs = 0;

   SCIPrationalDebugMessage("printing split disjunction <= %q \\/ >= %q+1  \n", mirinfo->rhs, mirinfo->rhs);

   for( i = 0; i < mirinfo->nsplitvars; ++i )
   {
      if( mirinfo->splitcoefficients[i] == 0 )
         continue;
      coefs[SCIPvarGetCertificateIndex(vars[mirinfo->varinds[i]])] += mirinfo->splitcoefficients[i];
      SCIPdebugMessage("+%g%s", mirinfo->splitcoefficients[i], SCIPvarGetName(vars[mirinfo->varinds[i]]));
   }
   for( i = 0; i < mirinfo->nslacks; ++i )
   {
      SCIP_ROW* slackrow = mirinfo->slackrows[i];
      SCIP_Real slackval = mirinfo->slackcoefficients[i];

      if( slackval == 0 )
         continue;

      for( j = 0; j < SCIProwGetNNonz(slackrow); ++j)
      {
         SCIP_VAR* var = SCIPcolGetVar(SCIProwGetCols(slackrow)[j]);
         int varidx = SCIPvarGetCertificateIndex(var);
         SCIP_Real rowcoef = SCIProwGetVals(slackrow)[j];

         assert(SCIPrealIsExactlyIntegral(rowcoef * slackval));
         assert(SCIPvarIsBinary(var) || SCIPvarIsIntegral(var));
         coefs[varidx] += rowcoef  * slackval;
         assert(SCIPrealIsExactlyIntegral(coefs[varidx]));
      }
      if( mirinfo->slacksign[i] == 1 )
      {
         assert(SCIPrealIsExactlyIntegral(SCIProwGetRhs(slackrow) - SCIProwGetConstant(slackrow)));
         slackrhs += (SCIProwGetRhs(slackrow) - SCIProwGetConstant(slackrow)) * slackval;

         assert(SCIPrealIsExactlyIntegral(slackrhs));
         SCIPdebugMessage("+%gsrhs_%s", mirinfo->slackcoefficients[i], SCIProwGetName(mirinfo->slackrows[i]));
      }
      else
      {
         assert(SCIPrealIsExactlyIntegral(SCIProwGetLhs(slackrow) - SCIProwGetConstant(slackrow)));
         slackrhs += (SCIProwGetLhs(slackrow) - SCIProwGetConstant(slackrow)) * slackval;

         assert(SCIPrealIsExactlyIntegral(slackrhs));
         SCIPdebugMessage("+%gslhs_%s", mirinfo->slackcoefficients[i], SCIProwGetName(mirinfo->slackrows[i]));
      }
   }

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &val) );

   /* transform the split back into original variable space -> undo the bound transformations */
   SCIPcertificatePrintProofMessage(certificate, "A%d_split %c ", certificate->indexcounter, 'L');

   /* add rhs change from integer slacks and print rhs */
   SCIPrationalAddReal(val, mirinfo->rhs, slackrhs);
   SCIP_CALL( SCIPcertificatePrintProofRational(certificate, val) );

   SCIPcertificatePrintProofMessage(certificate, " %d", coefs.size());

   for( auto & coef : coefs )
   {
      /** @todo perform line breaking before exceeding maximum line length */
      int varindex = coef.first;
      SCIPrationalSetReal(val, coef.second);

      assert(SCIPrationalIsIntegral(val));

      SCIPcertificatePrintProofMessage(certificate, " %d ", varindex);
      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, val) );
   }

   SCIPcertificatePrintProofMessage(certificate, " { asm } -1 \n");

   certificate->indexcounter++;
   certificate->lastinfo->isbound = FALSE;

   SCIPcertificatePrintProofMessage(certificate, "A%d_split %c ", certificate->indexcounter, 'G');

   SCIPrationalAddReal(val, mirinfo->rhs, slackrhs + 1.0);
   SCIP_CALL( SCIPcertificatePrintProofRational(certificate, val) );

   SCIPcertificatePrintProofMessage(certificate, " %d", coefs.size());

   for( auto & coef : coefs )
   {
      /** @todo perform line breaking before exceeding maximum line length */
      int varindex = coef.first;
      SCIPrationalSetReal(val, coef.second);

      assert(SCIPrationalIsIntegral(val));

      SCIPcertificatePrintProofMessage(certificate, " %d ", varindex);
      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, val) );
   }

   SCIPcertificatePrintProofMessage(certificate, " { asm } -1 \n");

   certificate->indexcounter++;
   certificate->lastinfo->isbound = FALSE;

   SCIPrationalFreeBuffer(set->buffer, &val);

   return SCIP_OKAY;
}

#ifdef SCIP_DISABLED_CODE
/** prints all local bounds that differ from their global bounds as the bounds to take into account */
static
SCIP_RETCODE certificatePrintIncompleteDerStart(
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_ROW**            rows,               /**< rows to be considered */
   int                   nrows,              /**< number of rows to be considered */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_Bool             local               /**< TRUE if the cut is only valid locally */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(!SCIPcertificateIsEnabled(certificate));

   vars = SCIPprobGetVars(prob);
   nvars = SCIPprobGetNVars(prob);

   SCIPcertificatePrintProofMessage(certificate, " { lin incomplete ");

   /* add non-global bounds */
   for( i = 0; i < nvars; i++ )
   {
      SCIP_VAR* var = vars[i];
      SCIP_Longint index;

      index = local ? SCIPvarGetLbCertificateIndexLocal(var) : SCIPvarGetLbCertificateIndexGlobal(var);
      if( index != -1 )
      {
	      assert(index <= certificate->indexcounter);
	      SCIPcertificatePrintProofMessage(certificate, " %d ", index);
      }
      index = local ? SCIPvarGetUbCertificateIndexLocal(var) : SCIPvarGetUbCertificateIndexGlobal(var);
      if( index != -1 )
      {
	      assert(index <= certificate->indexcounter);
	      SCIPcertificatePrintProofMessage(certificate, " %d ", index);
      }
   }

   /* add non-default rows */
   for( i = 0; i < nrows; ++i )
   {
      SCIP_ROWEXACT* rowexact = SCIProwGetRowExact(rows[i]);
      SCIP_Longint key;

      if( !local && SCIProwIsLocal(rows[i]) )
         continue;

      if( !SCIPhashmapExists(certificate->rowdatahash, (void*) rowexact) )
         continue;

      assert(rowexact != NULL);
      assert(SCIPhashmapExists(certificate->rowdatahash, (void*) rowexact));
      key = SCIPhashmapGetImageLong(certificate->rowdatahash, (void*) rowexact);
      SCIPcertificatePrintProofMessage(certificate, " %d ", key);
      /* for a ranged row, we need both sides to be safe */
      if( !SCIPrationalIsAbsInfinity(rowexact->lhs) && !SCIPrationalIsAbsInfinity(rowexact->rhs) && !SCIPrationalIsEQ(rowexact->lhs, rowexact->rhs) )
      {
         SCIPcertificatePrintProofMessage(certificate, " %d ", key + 1);
      }
   }

   return SCIP_OKAY;
}
#endif

/** prints all local bounds that differ from their global bounds as the bounds to take into account */
static
SCIP_RETCODE certificatePrintWeakDerStart(
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_Bool             local               /**< TRUE if the cut is only valid locally */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;
   int nboundentries;

   assert(SCIPcertificateIsEnabled(certificate));

   nboundentries = 0;

   vars = SCIPprobGetVars(prob);
   nvars = SCIPprobGetNVars(prob);

   /* count the number of needed entries */
   for( i = 0; i < nvars && local; i++ )
   {
      SCIP_VAR* var = vars[i];
      if( !SCIPrationalIsEQ(var->exactdata->glbdom.lb, var->exactdata->locdom.lb) )
         nboundentries++;
      if( !SCIPrationalIsEQ(var->exactdata->glbdom.ub, var->exactdata->locdom.ub) )
         nboundentries++;

      assert(!SCIPrationalIsEQ(var->exactdata->glbdom.lb, var->exactdata->locdom.lb) == (var->glbdom.lb != var->locdom.lb));
      assert(!SCIPrationalIsEQ(var->exactdata->glbdom.ub, var->exactdata->locdom.ub) == (var->glbdom.ub != var->locdom.ub));
   }

   SCIPcertificatePrintProofMessage(certificate, " { lin weak { %d", nboundentries);

   for( i = 0; i < nvars && nboundentries > 0; i++ )
   {
      SCIP_VAR* var = vars[i];

      if( !SCIPrationalIsEQ(var->exactdata->glbdom.lb, var->exactdata->locdom.lb) )
      {
         SCIP_Longint index;
         SCIP_RATIONAL* boundval;

         index = SCIPvarGetLbCertificateIndexLocal(var);
         boundval = SCIPvarGetLbLocalExact(var);
         SCIPcertificatePrintProofMessage(certificate, " L %d %d ", SCIPvarGetCertificateIndex(var), index);
         SCIP_CALL( SCIPcertificatePrintProofRational(certificate, boundval) );
      }
      if( !SCIPrationalIsEQ(var->exactdata->glbdom.ub, var->exactdata->locdom.ub) )
      {
         SCIP_Longint index;
         SCIP_RATIONAL* boundval;

         index = SCIPvarGetUbCertificateIndexLocal(var);
         boundval = SCIPvarGetUbLocalExact(var);
         SCIPcertificatePrintProofMessage(certificate, " U %d %d ", SCIPvarGetCertificateIndex(var), index);
         SCIP_CALL( SCIPcertificatePrintProofRational(certificate, boundval) );
      }
   }

   SCIPcertificatePrintProofMessage(certificate, " } ");

   return SCIP_OKAY;
}

/** prints verification of row as a MIR cut (viewed as a split cut) */
SCIP_RETCODE SCIPcertificatePrintMirCut(
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_LP*              lp,                 /**< SCIP lp data structure */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_ROW*             row,                /**< the row to be printed */
   const char            sense               /**< sense of the constraint, i.e., G, L, or E */
   )
{
   SCIP_ROWEXACT* rowexact;
   SCIP_RATIONAL* tmpval;
   SCIP_RATIONAL* value;
   SCIP_RATIONAL* oneminusf0;
   SCIP_AGGREGATIONINFO* aggrinfo;
   SCIP_Longint leftdisjunctionindex;
   SCIP_Longint rightdisjunctionindex;
   SCIP_MIRINFO* mirinfo;
   int i;

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   assert(row != NULL);
   assert(sense == 'L'); /* for now only this case is needed */
   assert(SCIPhashmapExists(certificate->aggrinfohash, (void*) row));

   rowexact = SCIProwGetRowExact(row);

   /* only do something if row does not already exist*/
   if( SCIPhashmapExists(certificate->rowdatahash, (void*) rowexact) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmpval) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &oneminusf0) );
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &value) );

   SCIPdebugMessage("Printing MIR-certifiaction for row ");
   SCIPdebug(SCIProwExactPrint(rowexact,  set->scip->messagehdlr, NULL));

   /* get aggregation info and print aggregation row to certificate */
   aggrinfo = (SCIP_AGGREGATIONINFO*) SCIPhashmapGetImage(certificate->aggrinfohash, (void*) row);
   mirinfo = (SCIP_MIRINFO*) SCIPhashmapGetImage(certificate->mirinfohash, (void*) row);

   /* compute the correct split from the aggregation row, and print the two assumptions (\xi \le \lfloor \beta \rfloor), and  (\xi \ge \lfloor \beta + 1 \rfloor) */
   SCIP_CALL( certificatePrintMirSplit(set, prob, certificate, row) );

   leftdisjunctionindex = certificate->indexcounter - 2;
   rightdisjunctionindex = certificate->indexcounter - 1;

   /* if this aggregation depends on another no yet certified MIR cut, we need to print that first */
   for( i = 0; i < aggrinfo->naggrrows; i++ )
   {
      SCIP_ROWEXACT* aggrrow;
      aggrrow = SCIProwGetRowExact(aggrinfo->aggrrows[i]);
        assert(aggrrow != NULL);
      if( !SCIPhashmapExists(certificate->rowdatahash, (void*) aggrrow) )
      {
         SCIP_CALL( SCIPcertificatePrintMirCut(set, lp, certificate, prob, aggrinfo->aggrrows[i], 'L') );
      }
   }

   /* print possibly missing derivations to the certifiacte */
   for( i = 0; i < aggrinfo->nnegslackrows; i++ )
   {
      SCIP_ROWEXACT* slackrow;
      slackrow = SCIProwGetRowExact(aggrinfo->negslackrows[i]);

      if( !SCIPhashmapExists(certificate->rowdatahash, (void*) slackrow) )
      {
         SCIP_CALL( SCIPcertificatePrintMirCut(set, lp, certificate, prob, aggrinfo->negslackrows[i], 'L') );
      }
   }

   /* print possibly missing derivations to the certifiacte */
   for( i = 0; i < mirinfo->nslacks; i++ )
   {
      SCIP_ROWEXACT* slackrow;
      slackrow = SCIProwGetRowExact(mirinfo->slackrows[i]);

      if( !SCIPhashmapExists(certificate->rowdatahash, (void*) slackrow) )
      {
         SCIP_CALL( SCIPcertificatePrintMirCut(set, lp, certificate, prob, mirinfo->slackrows[i], 'L') );
      }
   }

   /* print the mir cut with proof 1 * (\xi \le \lfloor \beta \rfloor) - (1/1-f)(\nu \ge 0) */
   /* we dont need the \nu \ge 0 part since it will be taken care of by the vipr completion part */
   assert(rowexact != NULL);
   {
      SCIP_CALL( certificatePrintRow(set, certificate, rowexact, mirinfo->unroundedrhs) );

      /* calculate 1-f0 */
      SCIPrationalSetReal(oneminusf0, 1.0);
      SCIPrationalDiff(oneminusf0, oneminusf0, mirinfo->frac);

      SCIP_CALL( certificatePrintWeakDerStart(certificate, prob, SCIProwIsLocal(row)) );

      /* 1 * (\xi \le \lfloor \beta \rfloor) we also have to add the correct multipliers for the negative slacks that were used here */
      SCIPcertificatePrintProofMessage(certificate, "%d %d ", 1 + aggrinfo->nnegslackrows + mirinfo->nslacks, leftdisjunctionindex);
      /* multiply with scaling parameter that was used during cut computation */
      SCIPrationalSetReal(tmpval, mirinfo->scale);
      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, tmpval) );

      SCIPdebugMessage("Verifying left part of split disjunction, multipliers 1 and 1/%g \n", SCIPrationalGetReal(oneminusf0));

      SCIPdebugMessage("Correcting for negative continous slacks ( needed for v >= 0 part ) \n");
      for( i = 0; i < aggrinfo->nnegslackrows; i++ )
      {
         SCIP_Longint key;
         SCIP_ROWEXACT* slackrow;
         slackrow = SCIProwGetRowExact(aggrinfo->negslackrows[i]);

         SCIPdebugMessage("adding (weight/(1-f0)) %g times row: ", aggrinfo->substfactor[i]);
         SCIPdebug(SCIProwExactPrint(slackrow, set->scip->messagehdlr, NULL));
         assert(slackrow != NULL);
         assert(SCIPhashmapExists(certificate->rowdatahash, (void*) slackrow));
         SCIPrationalSetReal(tmpval, aggrinfo->substfactor[i]);

         key = SCIPcertificateGetRowIndex(certificate, slackrow, SCIPrationalIsPositive(tmpval));

         SCIPcertificatePrintProofMessage(certificate, " %d ", key);
         SCIPrationalMultReal(tmpval, tmpval, mirinfo->scale);
         SCIP_CALL( SCIPcertificatePrintProofRational(certificate, tmpval) );
      }

      SCIPdebugMessage("Correcting for integer slacks  ( needed for v >= 0 part ) \n");
      for( i = 0; i < mirinfo->nslacks; i++ )
      {
         SCIP_Longint key;
         SCIP_ROWEXACT* slackrow;
         slackrow = SCIProwGetRowExact(mirinfo->slackrows[i]);
         SCIP_Longint upar;

         assert(slackrow != NULL);
         assert(SCIPhashmapExists(certificate->rowdatahash, (void*) slackrow));
         if( mirinfo->slackroundeddown[i] )
            SCIPrationalSetReal(tmpval, 0);
         else
         {
            /* value = weight * scale * sign -> compute (1 - fr)/(1-f0) where fr is the fractionality of value */
            SCIPrationalSetReal(tmpval, mirinfo->slackweight[i]);
            SCIPrationalMultReal(tmpval, tmpval, mirinfo->slackscale[i]);
            SCIPrationalMultReal(tmpval, tmpval, mirinfo->slacksign[i]);
            (void) SCIPrationalRoundLong(&upar, tmpval, SCIP_R_ROUND_UPWARDS);
            SCIPrationalDiffReal(tmpval, tmpval, upar);
            SCIPrationalNegate(tmpval, tmpval);
            SCIPrationalMultReal(tmpval, tmpval, mirinfo->slacksign[i]);
            SCIPrationalDiv(tmpval, tmpval, oneminusf0);
            SCIPrationalDebugMessage("tmpval 1 %g \n", SCIPrationalGetReal(tmpval));
            SCIPrationalSetReal(tmpval, mirinfo->slackusedcoef[i]);;
            SCIPrationalDiffReal(tmpval, tmpval, mirinfo->slackcoefficients[i]);
            SCIPrationalDebugMessage("tmpval 2 %g (slacksign %d, splitcoef %g, cutval %g) \n", SCIPrationalGetReal(tmpval), mirinfo->slacksign[i], mirinfo->slackcoefficients[i], mirinfo->slackusedcoef[i]);
         }

         key = SCIPcertificateGetRowIndex(certificate, slackrow, SCIPrationalIsPositive(tmpval));

         SCIPcertificatePrintProofMessage(certificate, " %d ", key);
         SCIPrationalMultReal(tmpval, tmpval, mirinfo->scale);
         SCIP_CALL( SCIPcertificatePrintProofRational(certificate, tmpval) );
      }

      SCIPcertificatePrintProofMessage(certificate, " } -1 \n");

      SCIPdebugMessage("Verifying right part of split disjunction, multipliers -f/(1-f) and 1/1-f \n");

      /* print the mir cut with proof (-f/1-f) * (\xi \ge \lfloor \beta + 1 \rfloor) + (1/1-f)(\xi - \nu \le \beta) */
      SCIP_CALL( certificatePrintRow(set, certificate, rowexact, mirinfo->unroundedrhs) );

      SCIP_CALL( certificatePrintWeakDerStart(certificate, prob, SCIProwIsLocal(row)) );
      SCIPcertificatePrintProofMessage(certificate, " %d ", 1 + aggrinfo->naggrrows + aggrinfo->nnegslackrows + mirinfo->nslacks);

      /* (-f/1-f) * (\xi \ge \lfloor \beta + 1 \rfloor) */
      SCIPcertificatePrintProofMessage(certificate, "%d ", rightdisjunctionindex);
      SCIPrationalSetRational(tmpval, mirinfo->frac);
      SCIPrationalNegate(tmpval, tmpval); /* -f */
      SCIPrationalDiv(tmpval, tmpval, oneminusf0); /* -f/(1-f) */
      /* multiply with scaling factor that was used in cut derivation */
      SCIPrationalMultReal(tmpval, tmpval, mirinfo->scale);
      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, tmpval) );

      SCIPrationalDivReal(tmpval, tmpval, mirinfo->scale);

      SCIPdebugMessage("Correcting for negative continous slacks \n");
      /* we also have to add the correct multipliers for the negative continuous slacks that were used here */
      for( i = 0; i < aggrinfo->nnegslackrows; i++ )
      {
         SCIP_Longint key;
         SCIP_ROWEXACT* slackrow;
         slackrow = SCIProwGetRowExact(aggrinfo->negslackrows[i]);

         SCIPrationalSetReal(value, -aggrinfo->negslackweights[i]);
         SCIPrationalDiv(value, value, oneminusf0);
         SCIPrationalAddReal(value, value, aggrinfo->substfactor[i]);
         SCIPrationalDebugMessage("adding %q times row (negative continous slacks) (%g aggweight %g substfactor): ", value, -aggrinfo->negslackweights[i] / SCIPrationalGetReal(oneminusf0), aggrinfo->substfactor[i]);
         SCIPdebug(SCIProwExactPrint(slackrow, set->scip->messagehdlr, NULL));

         assert(slackrow != NULL);
         assert(SCIPhashmapExists(certificate->rowdatahash, (void*) slackrow));

         key = SCIPcertificateGetRowIndex(certificate, slackrow, SCIPrationalIsPositive(value));

         SCIPcertificatePrintProofMessage(certificate, " %d ", key);
         SCIPrationalMultReal(value, value, mirinfo->scale);
         SCIP_CALL( SCIPcertificatePrintProofRational(certificate, value) );
      }

      SCIPdebugMessage("Correcting for %d integer slacks \n", mirinfo->nrounddownslacks);
      /* we also have to add the correct multipliers for the integer slacks that were used here */
      for( i = 0; i < mirinfo->nslacks; i++ )
      {
         SCIP_Longint key;
         SCIP_ROWEXACT* slackrow;
         slackrow = SCIProwGetRowExact(mirinfo->slackrows[i]);

         if( mirinfo->slackroundeddown[i] )
         {
            SCIPdebugMessage("Rounded down intger slack on row %s \n", mirinfo->slackrows[i]->name);
            SCIPrationalSetReal(value, mirinfo->slackweight[i]);
            SCIPrationalMultReal(value, value, mirinfo->slackscale[i]);
            SCIPrationalMultReal(value, value, mirinfo->slacksign[i]);
            SCIPrationalDiffReal(value, value, mirinfo->slackcoefficients[i] * (-mirinfo->slacksign[i]));
            SCIPrationalMultReal(value, value, mirinfo->slacksign[i]);
            SCIPrationalDiv(value, value, oneminusf0);
         }
         else
         {
            SCIPdebugMessage("Rounded up intger slack on row %s \n", mirinfo->slackrows[i]->name);
            SCIPrationalSetReal(value, mirinfo->slackweight[i]);
            SCIPrationalMultReal(value, value, mirinfo->slackscale[i]);
            SCIPrationalMultReal(value, value, mirinfo->slacksign[i]);
            SCIPrationalDiffReal(value, value, (mirinfo->slackcoefficients[i] * -mirinfo->slacksign[i]) - 1); /* fr exactly */
            SCIPrationalDiff(value, value, mirinfo->frac); /* fr - f0 */
            SCIPrationalDiv(value, value, oneminusf0); /* (fr-f0) / (1-f0) */
            SCIPrationalAddReal(value, value, (mirinfo->slackcoefficients[i] * -mirinfo->slacksign[i]) - 1); /* (down(ar) + (fr-f0) / (1-f0) */
            SCIPrationalMultReal(value, value, mirinfo->slacksign[i]);
            SCIPrationalDebugMessage("Exact coefficient %q(%g), used coefficient %g\n", value, SCIPrationalGetReal(value), mirinfo->slackusedcoef[i]);

            SCIPrationalAddReal(value, value, mirinfo->slackusedcoef[i]);
         }

         SCIPrationalDebugMessage("adding %q(%g) times row: ", value, SCIPrationalGetReal(value));
         SCIPdebug(SCIProwExactPrint(slackrow, set->scip->messagehdlr, NULL));

         assert(slackrow != NULL);
         assert(SCIPhashmapExists(certificate->rowdatahash, (void*) slackrow));

         key = SCIPcertificateGetRowIndex(certificate, slackrow, SCIPrationalIsPositive(value));

         SCIPcertificatePrintProofMessage(certificate, " %d ", key);
         SCIPrationalMultReal(value, value, mirinfo->scale);
         SCIP_CALL( SCIPcertificatePrintProofRational(certificate, value) );
      }

      SCIPdebugMessage("Adding %d aggregation rows \n", aggrinfo->naggrrows);
      /* we also have to add the correct multipliers for the aggregation rows that were used here */
      for( i = 0; i < aggrinfo->naggrrows; i++ )
      {
         SCIP_Longint key;
         SCIP_ROWEXACT* aggrrow;
         aggrrow = SCIProwGetRowExact(aggrinfo->aggrrows[i]);

         SCIPrationalSetReal(value, aggrinfo->weights[i]);
         SCIPrationalDiv(value, value, oneminusf0);

         SCIPdebugMessage("adding (%g/%g) = %g times row: ", aggrinfo->weights[i], SCIPrationalGetReal(oneminusf0), SCIPrationalGetReal(value));
         SCIPdebug(SCIProwExactPrint(aggrrow, set->scip->messagehdlr, NULL));

         assert(aggrrow != NULL);
         assert(SCIPhashmapExists(certificate->rowdatahash, (void*) aggrrow));

         key = SCIPcertificateGetRowIndex(certificate, aggrrow, SCIPrationalIsPositive(value));

         SCIPcertificatePrintProofMessage(certificate, " %d ", key);
         SCIPrationalMultReal(value, value, mirinfo->scale);
         SCIP_CALL( SCIPcertificatePrintProofRational(certificate, value) );
      }

      SCIPcertificatePrintProofMessage(certificate, " } -1\n");
   }
#ifdef SCIP_DISABLED_CODE
   /* this was an alternative piece of code to the block right above */
   {
      SCIP_CALL( certificatePrintRow(set, certificate, rowexact, mirinfo->unroundedrhs) );
      SCIP_CALL( certificatePrintIncompleteDerStart(certificate, SCIPlpGetRows(lp), SCIPlpGetNRows(lp), prob, SCIProwIsLocal(row)) );
      SCIPcertificatePrintProofMessage(certificate, " %d } -1\n", leftdisjunctionindex);

      SCIP_CALL( certificatePrintRow(set, certificate, rowexact, mirinfo->unroundedrhs) );
      SCIP_CALL( certificatePrintIncompleteDerStart(certificate, SCIPlpGetRows(lp), SCIPlpGetNRows(lp), prob, SCIProwIsLocal(row)) );
      SCIPcertificatePrintProofMessage(certificate, " %d } -1\n", rightdisjunctionindex);
   }
#endif

   /* print the unsplitting of the split disjunction */
   SCIP_CALL( certificatePrintRow(set, certificate, rowexact, mirinfo->unroundedrhs) );
   SCIPcertificatePrintProofMessage(certificate, " { uns %d %d  %d %d  } -1\n", certificate->indexcounter - 3, leftdisjunctionindex,
         certificate->indexcounter - 2, rightdisjunctionindex);

   /* if we rounded the rhs we need to still certify that part */
   if( !SCIPsetIsInfinity(set, mirinfo->unroundedrhs) )
   {
      /* print the row with the rounded rhs */
      SCIP_CALL( certificatePrintRow(set, certificate, rowexact, SCIPsetInfinity(set)) );
      SCIPcertificatePrintProofMessage(certificate, " { rnd 1 %d 1 } -1\n", certificate->indexcounter - 2);
   }
   SCIP_CALL( SCIPhashmapInsertLong(certificate->rowdatahash, SCIProwGetRowExact(row), certificate->indexcounter - 1) );


   SCIP_CALL( SCIPcertificateFreeAggrInfo(set, certificate, lp, aggrinfo, row) );
   SCIP_CALL( SCIPcertificateFreeMirInfo(set, certificate, lp, mirinfo, row) );

   SCIPrationalFreeBuffer(set->buffer, &value);
   SCIPrationalFreeBuffer(set->buffer, &oneminusf0);
   SCIPrationalFreeBuffer(set->buffer, &tmpval);

   return SCIP_OKAY;
}

/** prints a variable bound to the problem section of the certificate file and returns line index */
SCIP_RETCODE SCIPcertificatePrintBoundCons(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the line be printed to the origfile or the transfile */
   const char*           boundname,          /**< name of the bound constraint */
   SCIP_VAR*             var,                /**< variable to print the bound cons for */
   SCIP_RATIONAL*        boundval,           /**< value of the bound */
   SCIP_Bool             isupper             /**< is it the upper bound? */
   )
{
   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   if( !isorigfile )
   {
      certificate->indexcounter++;
      certificate->conscounter++;

#ifndef NDEBUG
      certificate->lastinfo->isbound = TRUE;
      certificate->lastinfo->boundtype = isupper ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER;
      certificate->lastinfo->varindex = SCIPvarGetCertificateIndex(var);
      certificate->lastinfo->isglobal = TRUE;
      certificate->lastinfo->certificateindex = certificate->indexcounter - 1;
      SCIPrationalSetRational(certificate->lastinfo->boundval, boundval);
#endif

      if( isupper )
         SCIPvarSetUbCertificateIndexGlobal(var, certificate->indexcounter - 1);
      else
         SCIPvarSetLbCertificateIndexGlobal(var, certificate->indexcounter - 1);

      if( boundname == NULL )
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, "B%d %c ", certificate->indexcounter - 1, (isupper ? 'L' : 'G'));
      else
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s %c ", boundname, (isupper ? 'L' : 'G'));
   }
   else
   {
      certificate->indexcounter_ori++;
      if( boundname == NULL )
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, "B%d %c ", certificate->indexcounter_ori - 1, (isupper ? 'L' : 'G'));
      else
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s %c ", boundname, (isupper ? 'L' : 'G'));
   }

   SCIP_CALL( SCIPcertificatePrintProblemRational(certificate, isorigfile, boundval) );
   SCIPcertificatePrintProblemMessage(certificate, isorigfile, " 1 %d 1\n", SCIPvarGetCertificateIndex(var));

   return SCIP_OKAY;
}

/** installs updated node data in parent node */
SCIP_RETCODE SCIPcertificateUpdateParentData(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_NODE*            node,               /**< node data structure */
   SCIP_Longint          fileindex,          /**< index of new bound */
   SCIP_RATIONAL*        newbound            /**< pointer to value of new bound, NULL if infeasible */
   )
{
   SCIP_CERTNODEDATA* nodedataparent;
   SCIP_CERTNODEDATA* nodedata;

   assert(node != NULL);
   assert(fileindex >= 0);

   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   /* retrieve node data */
   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);

   if( newbound != NULL && SCIPrationalIsLT(newbound, nodedata->derbound_self) )
      return SCIP_OKAY;

   /* if the node is the root node, then only update the index and bound */
   if( SCIPnodeGetParent(node) == NULL )
   {
      certificate->derindex_root = fileindex;
      if( newbound == NULL )
         certificate->rootinfeas = TRUE;
      else
         SCIPrationalSetRational(certificate->rootbound, newbound);
      return SCIP_OKAY;
   }

   /* retrieve parent node data */
   assert(SCIPhashmapExists(certificate->nodedatahash, SCIPnodeGetParent(node)));
   nodedataparent = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, SCIPnodeGetParent(node));

   /* First we check whether the node is left/right child of its parent node; left/rightfilled tells us if a bound has
    * already been derived for this node. We only install the new bound if the node is not already marked as infeasible
    * and the new bound is better than the current bound if it is filled.
    */
   if( certificateIsLeftNode(certificate, node) )
   {
      if( newbound != NULL && !nodedataparent->leftinfeas && (!nodedataparent->leftfilled || SCIPrationalIsGT(newbound, nodedataparent->derbound_left)) )
      {
         nodedataparent->derindex_left = fileindex;
         SCIPrationalSetRational(nodedataparent->derbound_left, newbound);
      }
      if( newbound == NULL || SCIPrationalIsInfinity(newbound) )
      {
         nodedataparent->derindex_left = fileindex;
         nodedataparent->leftinfeas = TRUE;
      }
      nodedataparent->leftfilled = TRUE;
   }
   else
   {
      if( newbound != NULL && !nodedataparent->rightinfeas && (!nodedataparent->rightfilled || SCIPrationalIsGT(newbound, nodedataparent->derbound_right)) )
      {
         nodedataparent->derindex_right = fileindex;
         SCIPrationalSetRational(nodedataparent->derbound_right, newbound);
      }
      if( newbound == NULL || SCIPrationalIsInfinity(newbound) )
      {
         nodedataparent->rightinfeas = TRUE;
         nodedataparent->derindex_right = fileindex;
      }
      nodedataparent->rightfilled = TRUE;
   }

   return SCIP_OKAY;
}

/** prints a dual bound from an exact lp solution */
SCIP_RETCODE SCIPcertificatePrintDualboundExactLP(
   SCIP_CERTIFICATE*     certificate,        /**< scip certificate struct */
   SCIP_LPEXACT*         lpexact,            /**< the exact lp */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_NODE*            node,               /**< the current node */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Bool             usefarkas           /**< should an infeasibility proof be printed? */
   )
{
   SCIP_RATIONAL** vals;
   SCIP_RATIONAL* val;
   SCIP_RATIONAL* lowerbound;
   SCIP_RATIONAL* farkasrhs;
   SCIP_RATIONAL* tmp;
   SCIP_Longint* ind = NULL;
   int len;
   int i;
   SCIP_Longint key;

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   SCIPdebugMessage("Printing dual bound from exact LP. Certificate index %lld \n", certificate->indexcounter);

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmp) );
   SCIPlpExactGetObjval(lpexact, set, tmp);

   /* only print line if bound improved */
   if( !usefarkas && SCIPrationalIsLT(tmp, SCIPnodeGetLowerboundExact(node)) )
   {
      SCIPrationalFreeBuffer(set->buffer, &tmp);
      return SCIP_OKAY;
   }

   /* init the flag objintegral if at root node */
   if( SCIPnodeGetParent(node) == NULL )
      certificate->objintegral = SCIPprobIsObjIntegral(prob);

   assert(lpexact != NULL);
   assert(certificate->transfile != NULL);

   /* if needed extend vals array */
   if( lpexact->ncols + lpexact->nrows > certificate->valssize )
   {
      SCIP_CALL( SCIPrationalReallocBlockArray(certificate->blkmem, &(certificate->vals), certificate->valssize, lpexact->nrows + lpexact->ncols + 50) );
      certificate->valssize =  lpexact->ncols + lpexact->nrows + 50;
   }

   vals = certificate->vals;

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &farkasrhs) );

   SCIP_CALL( SCIPsetAllocBufferArray(set, &ind, (lpexact->nrows + lpexact->ncols)) );

   len = 0;
   for( i = 0; i < lpexact->ncols; ++i )
   {
      SCIP_COLEXACT* col = lpexact->cols[i];
      SCIP_VAR* var = SCIPcolExactGetVar(col);

      if( usefarkas )
         val = col->farkascoef;
      else
         val = col->redcost;

      assert(!SCIPrationalIsAbsInfinity(val));

      if( !SCIPrationalIsZero(val) )
      {
         if( usefarkas )
            SCIPrationalNegate(vals[len], val);
         else
            SCIPrationalSetRational(vals[len], val);

         ind[len] = SCIPrationalIsNegative(vals[len]) ? SCIPvarGetUbCertificateIndexLocal(var) : SCIPvarGetLbCertificateIndexLocal(var);

         SCIPrationalDebugMessage("Column %d for var %s has index %l and farkas coef %q \n", col->index, col->var->name, ind[len], val);

         /* update farkasrhs */
         if( usefarkas )
         {
            val = SCIPrationalIsNegative(vals[len]) ? SCIPvarGetUbLocalExact(var) : SCIPvarGetLbLocalExact(var);
            SCIPrationalAddProd(farkasrhs, vals[len], val);
         }
         len++;
      }
   }

   for( i = 0; i < lpexact->nrows; ++i )
   {
      SCIP_ROWEXACT* row;
      row = lpexact->rows[i];
      val = usefarkas ? row->dualfarkas : row->dualsol;
      assert(!SCIPrationalIsAbsInfinity(val));

      if( !SCIPrationalIsZero(val) )
      {
         SCIPrationalSetRational(vals[len], val);
         key = SCIPhashmapGetImageLong(certificate->rowdatahash, (void*) row);

         if( key == SCIP_LONGINT_MAX && SCIProwGetOrigintype(SCIProwExactGetRow(row)) == SCIP_ROWORIGINTYPE_SEPA )
         {
            SCIP_CALL( SCIPcertificatePrintMirCut(set, lpexact->fplp, certificate, prob, SCIProwExactGetRow(row), 'L') );
            key = SCIPhashmapGetImageLong(certificate->rowdatahash, (void*) row);
         }
         else if( key == SCIP_LONGINT_MAX && SCIProwExactGetNNonz(row) == 1 )
         {
            SCIP_VAR* var = SCIPcolExactGetVar(SCIProwExactGetCols(row)[0]);
            key = SCIPrationalIsNegative(val) ? SCIPvarGetUbCertificateIndexLocal(var) : SCIPvarGetLbCertificateIndexLocal(var);
         }

         assert(key != SCIP_LONGINT_MAX);

         ind[len] = key;
         /* if we have a ranged row, and the dual corresponds to the upper bound,
          * the index for the rhs-constraint is one larger in the certificate */
         if( !SCIPrationalIsEQ(row->lhs, row->rhs) && !SCIPrationalIsAbsInfinity(row->lhs) && SCIPrationalIsNegative(val) )
             ind[len] += 1;

         SCIPrationalDebugMessage("Row (index %d, %s has index %l and farkas coef %q ", row->index, row->fprow->name, ind[len], val);
         SCIPdebug(SCIProwExactPrint(row, set->scip->messagehdlr, NULL) );

         /* update farkasrhs */
         if( usefarkas )
         {
            val = SCIPrationalIsPositive(vals[len]) ? row->lhs : row->rhs;
            SCIPrationalAddProd(farkasrhs, vals[len], val);
            SCIPrationalDiffProd(farkasrhs, vals[len], row->constant);
         }

         len++;
      }
   }

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &lowerbound) );
   if( usefarkas )
   {
      SCIPrationalSetInfinity(lowerbound);
   }
   else
   {
      /* vipr does not accept infinity, so in case of objlimit, get the objval from the lpi */
      if( !SCIPrationalIsInfinity(lpexact->lpobjval) )
         SCIPrationalSetRational(lowerbound, lpexact->lpobjval);
      else
      {
         SCIP_CALL( SCIPlpiExactGetObjval(lpexact->lpiexact, lowerbound) );
      }
   }

   SCIP_CALL( certificatePrintDualbound(certificate, NULL, lowerbound, len, ind, vals) );
   SCIP_CALL( SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1, lowerbound) );
   SCIP_CALL( SCIPcertificateUpdateBoundData(certificate, node, certificate->indexcounter - 1, lowerbound) );

   SCIPrationalFreeBuffer(set->buffer, &lowerbound);
   SCIPsetFreeBufferArray(set, &ind);
   SCIPrationalFreeBuffer(set->buffer, &farkasrhs);
   SCIPrationalFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** prints a dual bound from the pseudo solution
 *
 *  in case of a bound change (branching), this happens before the bound change is processed;
 *  therefore we add the option to give on varindex, boundchgindex pair to pass directly to the method
 */
SCIP_RETCODE SCIPcertificatePrintDualboundPseudo(
   SCIP_CERTIFICATE*     certificate,        /**< scip certificate struct */
   SCIP_LPEXACT*         lpexact,            /**< the exact lp */
   SCIP_NODE*            node,               /**< current node */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Bool             lowerchanged,       /**< do the modified indices address a change in lb or ub? */
   int                   modifiedvarindex,   /**< index of modified variable, or -1 */
   SCIP_Longint          boundchangeindex,   /**< index of unprocessed bound change in the certificate, or -1 */
   SCIP_Real             psval               /**< the pseudo obj value (or inf to use exact lp value) */
   )
{
   SCIP_VAR** vars;
   SCIP_RATIONAL* pseudoobjval;
   SCIP_RATIONAL** bounds;
   SCIP_Longint* dualind = NULL;
   int nvars;
   int duallen;
   int i;
   int nnonzeros;

   assert(certificate != NULL);
   assert(prob != NULL);
   assert(set != NULL);
   assert((modifiedvarindex >= 0 && boundchangeindex >= 0) || (modifiedvarindex == -1 && boundchangeindex == -1) );

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   /* only print if bound is finite and improving */
   if( SCIPsetIsInfinity(set, -psval) || psval < SCIPnodeGetLowerbound(node) )
      return SCIP_OKAY;

   /* if at root node set, objintegral flag */
   if( SCIPnodeGetParent(node) == NULL )
      certificate->objintegral = SCIPprobIsObjIntegral(prob);

   vars = SCIPprobGetVars(prob);
   nvars = SCIPprobGetNVars(prob);
   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &pseudoobjval) );

   /* infinity means we use the exact lp value */
   if( SCIPsetIsInfinity(set, psval) )
      SCIPlpExactGetPseudoObjval(lpexact, set, pseudoobjval);
   else
      SCIPrationalSetReal(pseudoobjval, psval);
   duallen = SCIPprobGetNObjVars(prob, set);
   SCIP_CALL( SCIPrationalCreateBufferArray(set->buffer, &bounds, duallen) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualind, duallen) );
   /* computes pseudo objective value (with bound change if necessary) */
   /*lint --e{838}*/

   nnonzeros = 0;
   for( i = 0; i < nvars; i++ )
   {
      SCIP_RATIONAL* obj = SCIPvarGetObjExact(vars[i]);
      if( !SCIPrationalIsZero(obj) )
      {
         SCIPrationalSetRational(bounds[nnonzeros], obj);

         assert(!SCIPrationalIsAbsInfinity(bounds[nnonzeros]));

         /* retrieve the line in the certificate of the bound */
         if( SCIPrationalIsPositive(obj) )
         {
            dualind[nnonzeros] = SCIPvarGetLbCertificateIndexLocal(vars[i]);
            if( lowerchanged && modifiedvarindex == SCIPvarGetCertificateIndex(vars[i]) )
               dualind[nnonzeros] = boundchangeindex;
         }
         else
         {
            dualind[nnonzeros] = SCIPvarGetUbCertificateIndexLocal(vars[i]);
            if( !lowerchanged && modifiedvarindex == SCIPvarGetCertificateIndex(vars[i]) )
               dualind[nnonzeros] = boundchangeindex;
         }

         nnonzeros++;
      }
   }
   assert(nnonzeros == duallen);

   /* print pseudo solution into certificate file */
   SCIP_CALL( certificatePrintDualbound(certificate, NULL, pseudoobjval, duallen, dualind, bounds) );

   SCIP_CALL( SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1,
      pseudoobjval) );

   SCIP_CALL( SCIPcertificateUpdateBoundData(certificate, node, certificate->indexcounter - 1, pseudoobjval) );

   SCIPsetFreeBufferArray(set, &dualind);
   SCIPrationalFreeBufferArray(set->buffer, &bounds, nnonzeros);
   SCIPrationalFreeBuffer(set->buffer, &pseudoobjval);

   return SCIP_OKAY;
}

/** prints the bound that a node inherits from its parent to the certificate */
SCIP_RETCODE SCIPcertificatePrintInheritedBound(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_NODE*            node                /**< node data */
   )
{
   SCIP_CERTNODEDATA* nodedata;
   SCIP_RATIONAL* lowerbound;

   assert(node != NULL);

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   /* certificate is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* get the current node data */
   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*) SCIPhashmapGetImage(certificate->nodedatahash, node);

   if( nodedata->inheritedbound && nodedata->assumptionindex_self != - 1 )
   {
      SCIP_Longint ind[1];
      SCIP_RATIONAL* val;

      ind[0] = nodedata->derindex_self;

      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &lowerbound) );
      SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &val) );

      SCIPrationalSetRational(lowerbound, nodedata->derbound_self);
      SCIPrationalSetFraction(val, 1LL, 1LL);

      SCIP_CALL( certificatePrintDualbound(certificate, NULL, lowerbound, 1, ind, &val) );
      SCIP_CALL( SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1, lowerbound) );

      SCIPrationalFreeBuffer(set->buffer, &lowerbound);
      SCIPrationalFreeBuffer(set->buffer, &val);
   }

   return SCIP_OKAY;
}

/** returns the index for a row in the certificate
 *
 *  @todo let this method return LONG_MAX if row is not in the hashmap; add method to check existence, and to insert an
 *        element, and use these throughout the SCIP core
 */
SCIP_Longint SCIPcertificateGetRowIndex(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_ROWEXACT*        row,                /**< row to consider */
   SCIP_Bool             rhs                 /**< whether we want the index for the rhs or the lhs */
   )
{
   SCIP_Longint ret = SCIPhashmapGetImageLong(certificate->rowdatahash, row);
   assert( ret != SCIP_LONGINT_MAX );
   /* for ranged rows, the key always corresponds to the >= part of the row;
         therefore we need to increase it by one to get the correct key */
   if( !SCIPrationalIsAbsInfinity(row->rhs) && !SCIPrationalIsAbsInfinity(row->lhs) && !SCIPrationalIsEQ(row->lhs, row->rhs) && rhs)
      ret += 1;
   return ret;
}

/** updates the parent certificate node data when branching */
SCIP_RETCODE SCIPcertificateUpdateBranchingData(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_LP*              lp,                 /**< LP informations */
   SCIP_NODE*            node,               /**< node data */
   SCIP_VAR*             branchvar,          /**< the variable that gets branched on */
   SCIP_BOUNDTYPE        boundtype,          /**< the bounding type */
   SCIP_Real             newbound            /**< the new bound */
   )
{
   SCIP_CERTNODEDATA* nodedataparent;
   SCIP_CERTNODEDATA* nodedata;
   SCIP_RATIONAL* branchbound;

   assert(node != NULL);
   assert(stat != NULL);

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   /* certificate is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE || lp->diving )
      return SCIP_OKAY;

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &branchbound) );
   SCIPrationalSetReal(branchbound, newbound);

   assert(SCIPrationalIsIntegral(branchbound));

   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);

   if( branchvar != NULL )
   {
      nodedata->assumptionindex_self = printBoundAssumption(certificate, branchvar,
         branchbound, boundtype);
   }

   if( SCIPnodeGetParent(node) != NULL && certificate->transfile != NULL )
   {
      nodedataparent = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, SCIPnodeGetParent(node));

      if( branchvar != NULL )
      {
         if( boundtype == SCIP_BOUNDTYPE_UPPER )
            nodedataparent->assumptionindex_right = nodedata->assumptionindex_self;
         if( boundtype == SCIP_BOUNDTYPE_LOWER )
            nodedataparent->assumptionindex_left = nodedata->assumptionindex_self;
      }
   }

   SCIPrationalFreeBuffer(set->buffer, &branchbound);

   return SCIP_OKAY;
}

/** create a new node data structure for the current node */
SCIP_RETCODE SCIPcertificateNewNodeData(
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   )
{
   SCIP_CERTNODEDATA* nodedata = NULL;

   assert(stat != NULL );
   assert(node != NULL );

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   /* certificate is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocBlockMemory(certificate->blkmem, &nodedata) );

   nodedata->derindex_left = -1;
   nodedata->derindex_right = -1;
   nodedata->assumptionindex_left = -1;
   nodedata->assumptionindex_right = -1;
   SCIP_CALL( SCIPrationalCreateString(certificate->blkmem, &nodedata->derbound_left, "-inf") );
   SCIP_CALL( SCIPrationalCreateString(certificate->blkmem, &nodedata->derbound_right, "-inf") );
   SCIP_CALL( SCIPrationalCreateString(certificate->blkmem, &nodedata->derbound_self, "-inf") );
   nodedata->assumptionindex_self = -1;
   nodedata->leftinfeas = FALSE;
   nodedata->leftfilled = FALSE;
   nodedata->rightinfeas = FALSE;
   nodedata->rightfilled = FALSE;
   nodedata->inheritedbound = TRUE;
   nodedata->derindex_self = -1;
   if( SCIPnodeGetParent(node) != NULL )
   {
      SCIP_NODE* parent = SCIPnodeGetParent(node);
      SCIP_CERTNODEDATA* parentdata;
      assert(SCIPhashmapExists(certificate->nodedatahash, parent));
      parentdata = (SCIP_CERTNODEDATA*) SCIPhashmapGetImage(certificate->nodedatahash, (void*) parent);
      assert(parentdata != NULL);

      nodedata->derindex_self = parentdata->derindex_self;
      SCIPrationalSetRational(nodedata->derbound_self, parentdata->derbound_self);
   }

   /* link the node to its nodedata in the corresponding hashmap */
   SCIP_CALL( SCIPhashmapSetImage(certificate->nodedatahash, node, (void*)nodedata) );

   return SCIP_OKAY;
}

/** prints cutoff bound for objective value **/
SCIP_RETCODE SCIPcertificatePrintCutoffBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_RATIONAL*        bound,              /**< the bound */
   SCIP_Longint*         certificateline     /**< save the line index */
   )
{
   SCIP_RATIONAL* newbound;

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &newbound) );
   if( SCIPisObjIntegral(scip) && ! SCIPrationalIsIntegral(bound) )
      SCIPrationalRoundInteger(newbound, bound, SCIP_R_ROUND_DOWNWARDS);
   else
      SCIPrationalSetRational(newbound, bound);

   SCIPcertificatePrintProofMessage(certificate, "O%d L ", certificate->indexcounter);
   SCIP_CALL( SCIPcertificatePrintProofRational(certificate, newbound) );
   SCIPcertificatePrintProofMessage(certificate, " OBJ { sol } -1\n");

   certificate->indexcounter++;
   *certificateline = certificate->indexcounter - 1;
   SCIPrationalFreeBuffer(SCIPbuffer(scip), &newbound);

   return SCIP_OKAY;
}

/** create a new node data structure for the current node */
SCIP_RETCODE SCIPcertificatePrintAggrrow(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_PROB*            prob,               /**< SCIP problem data */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_AGGRROW*         aggrrow,            /**< agrrrow that results from the aggregation */
   SCIP_ROW**            aggrrows,           /**< array of rows used fo the aggregation */
   SCIP_Real*            weights,            /**< array of weights */
   int                   naggrrows,          /**< length of the arrays */
   SCIP_Bool             local,              /**< true if local bound information can be used */
   SCIP_Longint*         certificateline     /**< pointer to store the certificate line index or NULL */
   )
{
   int i;
   SCIP_RATIONAL* tmpval;
   SCIP_ROWEXACT* rowexact;
   SCIP_VAR** vars;

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPrationalCreateBuffer(set->buffer, &tmpval) );
   vars = SCIPprobGetVars(prob);

   SCIPdebugMessage("printing certificate for aggrrow: ");
   SCIPdebug(SCIPaggrRowPrint(set->scip, aggrrow, NULL));

   SCIPcertificatePrintProofMessage(certificate, "AggrRow_%d %c ", certificate->indexcounter, 'L');

   SCIPrationalSetReal(tmpval, SCIPaggrRowGetRhs(aggrrow));

   SCIP_CALL( SCIPcertificatePrintProofRational(certificate, tmpval) );

   SCIPcertificatePrintProofMessage(certificate, " %d", SCIPaggrRowGetNNz(aggrrow));

   for( i = 0; i < SCIPaggrRowGetNNz(aggrrow); i++ )
   {
      int varindex;

      /** @todo perform line breaking before exceeding maximum line length */
      varindex = SCIPvarGetCertificateIndex(vars[SCIPaggrRowGetInds(aggrrow)[i]]);
      SCIPrationalSetReal(tmpval, SCIPaggrRowGetValueSafely(aggrrow, i));

      SCIPcertificatePrintProofMessage(certificate, " %d ", varindex);
      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, tmpval) );
   }

   SCIP_CALL( certificatePrintWeakDerStart(certificate, prob, local) );
   SCIPcertificatePrintProofMessage(certificate, " %d", naggrrows);
   for( i = 0; i < naggrrows; i++ )
   {
      SCIP_Longint key;
      rowexact = SCIProwGetRowExact(aggrrows[i]);

      SCIPdebugMessage("adding %g times row: ", weights[i]);
      SCIPdebug(SCIProwExactPrint(rowexact, set->scip->messagehdlr, NULL));
      SCIPrationalSetReal(tmpval, weights[i]);

      assert(rowexact != NULL);
      assert(SCIPhashmapExists(certificate->rowdatahash, (void*) rowexact));

      key = SCIPcertificateGetRowIndex(certificate, rowexact, SCIPrationalIsPositive(tmpval));

      SCIPcertificatePrintProofMessage(certificate, " %d ", key);
      SCIP_CALL( SCIPcertificatePrintProofRational(certificate, tmpval) );
   }

   SCIPcertificatePrintProofMessage(certificate, " } -1\n");

   if( certificateline != NULL )
      *certificateline = certificate->indexcounter;
   certificate->indexcounter++;
   certificate->lastinfo->isbound = FALSE;

   SCIPrationalFreeBuffer(set->buffer, &tmpval);

   return SCIP_OKAY;
}

/** free all aggregation information */
SCIP_RETCODE SCIPcertificateClearAggrinfo(
   SCIP*                 scip                /**< global SCIP data structure */
   )
{
   SCIP_Longint i;
   SCIP_CERTIFICATE* certificate;

   certificate = SCIPgetCertificate(scip);

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   if( certificate == NULL || certificate->aggrinfo == NULL )
      return SCIP_OKAY;

   assert(certificate != NULL);

   for( i = certificate->naggrinfos - 1; i >= 0; i-- )
   {
      SCIP_CALL( SCIPcertificateFreeAggrInfo(scip->set, certificate, scip->lp, certificate->aggrinfo[i], NULL) );
   }

   BMSfreeBlockMemoryArray(certificate->blkmem, &certificate->aggrinfo, certificate->aggrinfosize);
   SCIP_CALL( SCIPhashmapRemoveAll(certificate->aggrinfohash) );
   SCIPhashmapFree(&certificate->aggrinfohash);
   certificate->naggrinfos = 0;
   certificate->aggrinfohash = NULL;

   return SCIP_OKAY;
}

/** free all mir information */
SCIP_RETCODE SCIPcertificateClearMirinfo(
   SCIP*                 scip                /**< global SCIP data structure */
   )
{
   int i, j;
   SCIP_CERTIFICATE* certificate;

   certificate = SCIPgetCertificate(scip);

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   if( certificate == NULL || certificate->mirinfo == NULL )
      return SCIP_OKAY;

   assert(certificate != NULL);

   for( i = 0; i < certificate->nmirinfos; i++ )
   {
      for( j = 0; j < certificate->mirinfo[i]->nslacks; ++j )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &(certificate->mirinfo[i]->slackrows[j])) );
      }
      SCIPrationalFreeBlock(certificate->blkmem, &certificate->mirinfo[i]->frac);
      SCIPrationalFreeBlock(certificate->blkmem, &certificate->mirinfo[i]->rhs);
      BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->varinds), certificate->mirinfo[i]->nsplitvars);
      BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->splitcoefficients), certificate->mirinfo[i]->nsplitvars);
      BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->upperused), certificate->mirinfo[i]->nsplitvars);
      BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->localbdused), certificate->mirinfo[i]->nsplitvars);
      BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slackrows), certificate->mirinfo[i]->nslacks);
      BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slacksign), certificate->mirinfo[i]->nslacks);
      BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slackcoefficients), certificate->mirinfo[i]->nslacks);
      BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slackweight), certificate->mirinfo[i]->nslacks);
      BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slackroundeddown), certificate->mirinfo[i]->nslacks);
      BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slackscale), certificate->mirinfo[i]->nslacks);
      BMSfreeBlockMemoryArray(certificate->blkmem, &(certificate->mirinfo[i]->slackusedcoef), certificate->mirinfo[i]->nslacks);
      BMSfreeBlockMemory(certificate->blkmem, &certificate->mirinfo[i]);
   }

   certificate->nmirinfos = 0;
   BMSfreeBlockMemoryArray(certificate->blkmem, &certificate->mirinfo, certificate->mirinfosize);
   SCIP_CALL( SCIPhashmapRemoveAll(certificate->mirinfohash) );
   SCIPhashmapFree(&certificate->mirinfohash);

   return SCIP_OKAY;
}

/** free aggregation information for row */
SCIP_RETCODE SCIPcertificateFreeAggrInfo(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate structure */
   SCIP_LP*              lp,                 /**< SCIP lp data structure */
   SCIP_AGGREGATIONINFO* aggrinfo,           /**< SCIP aggregation info */
   SCIP_ROW*             row                 /**< row that should be freed, or NULL if not needed */
   )
{
   SCIP_Longint arraypos;
   int i;

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   /* remove the (no longer needed aggrinfo), move last element to now freed spot */
   arraypos = aggrinfo->arpos;

   if( row != NULL )
   {
      SCIP_CALL( SCIPhashmapRemove(certificate->aggrinfohash, (void*) row) );
   }

   for( i = 0; i < aggrinfo->naggrrows; i++ )
   {
      SCIP_CALL( SCIProwRelease(&(aggrinfo->aggrrows[i]), certificate->blkmem, set, lp) );
   }
   for( i = 0; i < aggrinfo->nnegslackrows; i++ )
   {
      SCIP_CALL( SCIProwRelease(&(aggrinfo->negslackrows[i]), certificate->blkmem, set, lp) );
   }

   BMSfreeBlockMemoryArray(certificate->blkmem, &(aggrinfo->substfactor), aggrinfo->nnegslackrows);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(aggrinfo->negslackweights), aggrinfo->nnegslackrows);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(aggrinfo->negslackrows), aggrinfo->nnegslackrows);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(aggrinfo->weights), aggrinfo->naggrrows);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(aggrinfo->aggrrows), aggrinfo->naggrrows);
   SCIPaggrRowFree(set->scip, &(aggrinfo->aggrrow));
   BMSfreeBlockMemory(certificate->blkmem, &aggrinfo);
   if( arraypos != certificate->naggrinfos - 1 )
   {
      certificate->aggrinfo[arraypos] = certificate->aggrinfo[certificate->naggrinfos - 1];
      certificate->aggrinfo[arraypos]->arpos = arraypos;
   }
   certificate->naggrinfos--;

   return SCIP_OKAY;
}

/** free mir information for row */
SCIP_RETCODE SCIPcertificateFreeMirInfo(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate structure */
   SCIP_LP*              lp,                 /**< SCIP lp data structure */
   SCIP_MIRINFO*         mirinfo,            /**< SCIP mir info */
   SCIP_ROW*             row                 /**< row that should be freed, or NULL if not needed */
   )
{
   SCIP_Longint arraypos;
   int i;

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   /* remove the mirinfo, move last element to the now freed up one */
   arraypos = mirinfo->arpos;
   SCIP_CALL( SCIPhashmapRemove(certificate->mirinfohash, (void*) row) );
   for( i = 0; i < mirinfo->nslacks; ++i)
   {
      SCIP_CALL( SCIProwRelease(&(mirinfo->slackrows[i]), certificate->blkmem, set, lp) );
   }

   SCIPrationalFreeBlock(certificate->blkmem, &(mirinfo->rhs));
   SCIPrationalFreeBlock(certificate->blkmem, &(mirinfo->frac));
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->splitcoefficients), mirinfo->nsplitvars);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->varinds), mirinfo->nsplitvars);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->upperused), mirinfo->nsplitvars);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->localbdused), mirinfo->nsplitvars);

   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slackrows), mirinfo->nslacks);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slacksign), mirinfo->nslacks);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slackcoefficients), mirinfo->nslacks);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slackweight), mirinfo->nslacks);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slackscale), mirinfo->nslacks);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slackusedcoef), mirinfo->nslacks);
   BMSfreeBlockMemoryArray(certificate->blkmem, &(mirinfo->slackroundeddown), mirinfo->nslacks);

   BMSfreeBlockMemory(certificate->blkmem, &mirinfo);
   if( arraypos != certificate->nmirinfos - 1 )
   {
      certificate->mirinfo[arraypos] = certificate->mirinfo[certificate->nmirinfos - 1];
      certificate->mirinfo[arraypos]->arpos = arraypos;
   }
   certificate->nmirinfos--;

   return SCIP_OKAY;
}

/** free information that is possibly still stored about this row in the certificate structure */
SCIP_RETCODE SCIPcertificateFreeRowInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< a SCIP row */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_AGGREGATIONINFO* aggrinfo;
   SCIP_MIRINFO* mirinfo;

   if( !SCIPisExact(scip) || !SCIPisCertified(scip) )
      return SCIP_OKAY;

   certificate = SCIPgetCertificate(scip);

   /* only do something if row does not already exist*/
   if( SCIPhashmapExists(certificate->rowdatahash, (void*) SCIProwGetRowExact(row)) )
   {
      SCIP_CALL( SCIPhashmapRemove(certificate->rowdatahash, (void*) SCIProwGetRowExact(row)) );
      return SCIP_OKAY;
   }

   if( certificate->workingaggrinfo || certificate->workingmirinfo )
      return SCIP_OKAY;

   SCIPdebugMessage("Removing information stored in certificate for row \n");

   if( (certificate->aggrinfohash != NULL) && SCIPhashmapExists(certificate->aggrinfohash, (void*) row) )
   {
      aggrinfo = (SCIP_AGGREGATIONINFO*) SCIPhashmapGetImage(certificate->aggrinfohash, (void*) row);
      SCIP_CALL( SCIPcertificateFreeAggrInfo(scip->set, certificate, scip->lp, aggrinfo, row) );
   }

   if( (certificate->mirinfohash != NULL) && SCIPhashmapExists(certificate->mirinfohash, (void*) row) )
   {
      mirinfo = (SCIP_MIRINFO*) SCIPhashmapGetImage(certificate->mirinfohash, (void*) row);
      SCIP_CALL( SCIPcertificateFreeMirInfo(scip->set, certificate, scip->lp, mirinfo, row) );
   }

   return SCIP_OKAY;
}

/** create a new node data structure for the current node */
SCIP_RETCODE SCIPcertificateNewAggrInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< agrrrow that results from the aggregation */
   SCIP_ROW**            aggrrows,           /**< array of rows used fo the aggregation */
   SCIP_Real*            weights,            /**< array of weights */
   int                   naggrrows,          /**< length of the arrays */
   SCIP_ROW**            negslackrows,       /**< array of rows that are added implicitly with negative slack */
   SCIP_Real*            negslackweights,    /**< array of negative slack weights */
   int                   nnegslackrows       /**< length of the negative slack array */
   )
{
   int i;
   SCIP_AGGREGATIONINFO* info = NULL;
   SCIP_CERTIFICATE* certificate;

   assert(scip != NULL );

   certificate = SCIPgetCertificate(scip);

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   assert(certificate != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &info) );

   SCIPdebugMessage("adding aggrinfo, with %d rows to certficate \n", naggrrows);

   SCIP_CALL( SCIPaggrRowCopy(scip, &(info->aggrrow), aggrrow) );

   info->naggrrows = naggrrows;
   info->nnegslackrows = nnegslackrows;
   info->fileindex = certificate->indexcounter - 1;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(info->aggrrows), naggrrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(info->weights), naggrrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(info->negslackrows), nnegslackrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(info->negslackweights), nnegslackrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(info->substfactor), nnegslackrows) );
   for( i = 0; i < naggrrows; i++ )
   {
      SCIPdebugMessage("adding row %s with weight %g to aggrinfo \n", SCIProwGetName(aggrrows[i]), weights[i]);
      info->aggrrows[i] = aggrrows[i];
      info->weights[i] = weights[i];
      SCIProwCapture(aggrrows[i]);
   }
   for( i = 0; i < nnegslackrows; i++ )
   {
      SCIPdebugMessage("adding (implicitly) row %s with negative slack multiplier %g to aggrinfo \n", SCIProwGetName(negslackrows[i]), negslackweights[i]);
      info->negslackrows[i] = negslackrows[i];
      info->negslackweights[i] = negslackweights[i];
      SCIProwCapture(negslackrows[i]);
   }

   /* link the node to its nodedata in the corresponding hashmap */
   if( certificate->aggrinfosize == certificate->naggrinfos )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &certificate->aggrinfo, certificate->aggrinfosize, certificate->aggrinfosize + 100) );
      certificate->aggrinfosize += 100;
   }
   certificate->aggrinfo[certificate->naggrinfos] = info;
   info->arpos = certificate->naggrinfos;
   certificate->naggrinfos++;
   certificate->workingaggrinfo = TRUE;

   return SCIP_OKAY;
}

/** create a new split info structure for the current cut */
SCIP_RETCODE SCIPcertificateNewMirInfo(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_MIRINFO* mirinfo;
   SCIP_CERTIFICATE* certificate;

   assert(scip != NULL );

   certificate = SCIPgetCertificate(scip);

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   if( certificate->mirinfosize == certificate->nmirinfos )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &certificate->mirinfo, certificate->mirinfosize, certificate->mirinfosize + 100) );
      certificate->mirinfosize += 100;
   }

   SCIP_CALL( SCIPallocBlockMemory(scip, &(certificate->mirinfo[certificate->nmirinfos])) );

   mirinfo =  certificate->mirinfo[certificate->nmirinfos];

   SCIPdebugMessage("adding mirinfo, with to certficate \n");

   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip), &(mirinfo->rhs)) );
   SCIP_CALL( SCIPrationalCreateBlock(SCIPblkmem(scip), &(mirinfo->frac)) );
   mirinfo->nlocalvars = 0;
   mirinfo->nslacks = 0;
   mirinfo->nrounddownslacks = 0;
   mirinfo->nsplitvars = SCIPgetNVars(scip);
   mirinfo->arpos = certificate->nmirinfos;
   mirinfo->scale = 1.0;
   mirinfo->unroundedrhs = SCIPinfinity(scip);

   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(mirinfo->splitcoefficients), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(mirinfo->upperused), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(mirinfo->localbdused), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slackrows), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slackcoefficients), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slacksign), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slackweight), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slackscale), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slackusedcoef), SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mirinfo->slackroundeddown), SCIPgetNVars(scip)) );
   mirinfo->varinds = NULL;

   certificate->nmirinfos++;
   certificate->workingmirinfo = TRUE;

   return SCIP_OKAY;
}

/** prints unsplitting information to proof section */
SCIP_RETCODE SCIPcertificatePrintUnsplitting(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_NODE*            node                /**< node data */
   )
{
   SCIP_CERTNODEDATA* nodedata;
   SCIP_RATIONAL* lowerbound;
   SCIP_Bool infeas;

   assert(node != NULL);

   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   /* certificate is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* get the current node data */
   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*) SCIPhashmapGetImage(certificate->nodedatahash, node);
   (void) SCIPrationalCreateBuffer(set->buffer, &lowerbound);
   infeas = FALSE;

   assert(nodedata != NULL);

   if( nodedata->leftfilled && nodedata->rightfilled )
   {
      if( nodedata->leftinfeas && nodedata->rightinfeas )
      {
         infeas = TRUE;
         SCIPrationalSetInfinity(lowerbound);
      }
      else if( nodedata->leftinfeas )
         SCIPrationalSetRational(lowerbound, nodedata->derbound_right);
      else if( nodedata->rightinfeas )
         SCIPrationalSetRational(lowerbound, nodedata->derbound_left);
      else
         SCIPrationalMin(lowerbound, nodedata->derbound_left, nodedata->derbound_right);

      if( SCIPrationalIsInfinity(nodedata->derbound_left) && SCIPrationalIsInfinity(nodedata->derbound_right) )
         infeas = TRUE;

      certificate->indexcounter++;
      certificate->lastinfo->isbound = FALSE;

      SCIPcertificatePrintProofMessage(certificate, "UnsplitNode%d_%d ", SCIPnodeGetNumber(node), certificate->indexcounter - 1);

      if( infeas )
      {
         SCIPcertificatePrintProofMessage(certificate, "G 1 0");
      }
      else
      {
         SCIPcertificatePrintProofMessage(certificate, "G ");
         SCIP_CALL( SCIPcertificatePrintProofRational(certificate, lowerbound) );
         SCIPcertificatePrintProofMessage(certificate, " ");
         SCIPcertificatePrintProofMessage(certificate, "OBJ");
      }

      assert(nodedata->derindex_right < certificate->indexcounter - 1);
      assert(nodedata->derindex_left < certificate->indexcounter - 1);
      assert(SCIPrationalIsGE(lowerbound, SCIPnodeGetLowerboundExact(node)));
      assert(SCIPrationalIsGEReal(lowerbound, SCIPnodeGetLowerbound(node)));

      SCIPcertificatePrintProofMessage(certificate, " { uns %d %d  %d %d  } -1\n", nodedata->derindex_left, nodedata->assumptionindex_left,
         nodedata->derindex_right, nodedata->assumptionindex_right);

      SCIP_CALL( SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1, lowerbound) );
   }
   else
   {
      SCIPdebugMessage("Node %lld is a leaf! \n", SCIPnodeGetNumber(node));
      /* if a leaf has an inherited bound, we need to print a bound for it and update the parent data
         don't do it if we interrupted the solve, e.g. due to timeout */
      if( nodedata->inheritedbound && nodedata->assumptionindex_self != - 1 )
      {
         SCIP_Longint ind[1];
         SCIP_RATIONAL* val;

         ind[0] = nodedata->derindex_self;

         (void) SCIPrationalCreateBuffer(set->buffer, &val);

         SCIPrationalSetRational(lowerbound, nodedata->derbound_self);
         SCIPrationalSetFraction(val, 1LL, 1LL);

         SCIP_CALL( certificatePrintDualbound(certificate, NULL, lowerbound, 1, ind, &val) );
         SCIP_CALL( SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1, lowerbound) );

         SCIPrationalFreeBuffer(set->buffer, &val);
      }
   }

   SCIP_CALL( certificateFreeNodeData(certificate, node) );

   SCIPrationalFreeBuffer(set->buffer, &lowerbound);

   return SCIP_OKAY;
}

/** prints RTP section with lowerbound and upperbound range */
SCIP_RETCODE SCIPcertificatePrintRtpRange(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   SCIP_RATIONAL*        lowerbound,         /**< pointer to lower bound on the objective */
   SCIP_RATIONAL*        upperbound          /**< pointer to upper bound on the objective */
   )
{
   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return SCIP_OKAY;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "RTP range ");
   if( SCIPrationalIsNegInfinity(lowerbound) )
   {
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "-inf");
   }
   else
   {
      SCIP_CALL( SCIPcertificatePrintProblemRational(certificate, isorigfile, lowerbound) );
   }

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, " ");
   if( SCIPrationalIsInfinity(upperbound) )
   {
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "inf");
   }
   else
   {
      SCIP_CALL( SCIPcertificatePrintProblemRational(certificate, isorigfile, upperbound) );
   }
   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "\n");

   return SCIP_OKAY;
 }

/** prints RTP section for infeasibility */
void SCIPcertificatePrintRtpInfeas(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_Bool             isorigfile          /**< should the original solution be printed or in transformed space */
   )
{
   /* check whether certificate output should be created */
   if( !SCIPcertificateIsEnabled(certificate) )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "RTP infeas\n");
 }

/** sets the last bound index for the certificate */
SCIP_RETCODE SCIPcertificateSetLastBoundIndex(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_Longint          index               /**< index of new bound */
   )
{
   assert(index >= 0);

   if( SCIPcertificateIsEnabled(certificate) )
      certificate->lastboundindex = index;

   return SCIP_OKAY;
}

/** returns the last bound index for the certificate */
SCIP_Longint SCIPcertificateGetLastBoundIndex(
   SCIP_CERTIFICATE*     certificate         /**< certificate data structure */
   )
{
   return SCIPcertificateIsEnabled(certificate) ? certificate->lastboundindex : -1;
}

/** prints a proof that boundchange is leads to infeasibility */
SCIP_RETCODE SCIPcertificatePrintCutoffConflictingBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_VAR*             var,                /**< variable */
   SCIP_RATIONAL*        lb,                 /**< lower bound */
   SCIP_RATIONAL*        ub,                 /**< upper bound */
   SCIP_Longint          lbindex,            /**< index of the lower bound */
   SCIP_Longint          ubindex             /**< index of the upper bound */
   )
{
   SCIP_RATIONAL* lowerbound;

   if( !SCIPisCertified(scip) )
      return SCIP_OKAY;

   assert(certificate != NULL);

   switch( var->varstatus )
   {
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_ORIGINAL:
         SCIPABORT();
         return SCIP_ERROR;
      case SCIP_VARSTATUS_NEGATED:
         if( lb != NULL )
         {
            SCIPrationalMultReal(lb, lb, -1);
            SCIPrationalAddReal(lb, lb, 1.0);
         }
         if( ub != NULL )
         {
            SCIPrationalMultReal(ub, ub, -1);
            SCIPrationalAddReal(ub, ub, 1.0);
         }
         assert( SCIPvarGetNegationConstant(var) == 1 );
         SCIP_CALL( SCIPcertificatePrintCutoffConflictingBounds(scip, certificate, var->negatedvar, ub, lb, ubindex, lbindex) );
         if( lb != NULL )
         {
            SCIPrationalAddReal(lb, lb, -1.0);
            SCIPrationalMultReal(lb, lb, -1);
         }
         if( ub != NULL )
         {
            SCIPrationalAddReal(ub, ub, -1.0);
            SCIPrationalMultReal(ub, ub, -1);
         }
         return SCIP_OKAY;
         break;
      case SCIP_VARSTATUS_AGGREGATED:
         {
            if( lb != NULL )
               SCIPrationalDiv(lb, lb, var->exactdata->aggregate.scalar);
            if( ub != NULL )
               SCIPrationalDiv(ub, ub, var->exactdata->aggregate.scalar);

            assert(SCIPrationalIsZero(var->exactdata->aggregate.constant));
            SCIP_Bool swapBounds = SCIPrationalIsPositive(var->exactdata->aggregate.scalar) ? FALSE: TRUE;
            SCIP_CALL( SCIPcertificatePrintCutoffConflictingBounds(scip, certificate, var->data.aggregate.var, swapBounds ? ub : lb, swapBounds ? lb : ub,  swapBounds ? ubindex : lbindex, swapBounds ? lbindex : ubindex) );
            if( lb != NULL )
               SCIPrationalMult(lb, lb, var->exactdata->aggregate.scalar);
            if( ub != NULL )
               SCIPrationalMult(ub, ub, var->exactdata->aggregate.scalar);
            return SCIP_OKAY;
         }
         break;
      case SCIP_VARSTATUS_COLUMN:
         break;
      default:
         SCIPABORT();
         return SCIP_ERROR;
   }

   if( lb == NULL )
   {
      lb = SCIPvarGetLbLocalExact(var);
      lbindex = SCIPvarGetLbCertificateIndexLocal(var);
   }
   if( ub == NULL )
   {
      ub = SCIPvarGetUbLocalExact(var);
      ubindex = SCIPvarGetUbCertificateIndexLocal(var);
   }
   assert( SCIPrationalIsGT(lb, ub) );
   SCIP_CALL(SCIPrationalCreateBuffer(SCIPbuffer(scip), &lowerbound));

   SCIPcertificatePrintProofMessage(certificate, "BoundConflict%d ", certificate->indexcounter);
   SCIPcertificatePrintProofMessage(certificate, "G ");
   SCIPrationalDiff(lowerbound, lb, ub);
   SCIP_CALL( SCIPcertificatePrintProofRational(certificate, lowerbound) );
   SCIPcertificatePrintProofMessage(certificate, " 0 { lin 2 %d 1 %d -1 } -1\n", lbindex, ubindex);

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &lowerbound);

   SCIP_CALL( SCIPcertificateUpdateParentData(certificate, SCIPgetCurrentNode(scip), certificate->indexcounter, NULL) );
   certificate->indexcounter++;

   return SCIP_OKAY;
}

/** prints a proof for a new global bound */
SCIP_RETCODE SCIPcertificatePrintGlobalBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_VAR*             var,                /**< variable */
   SCIP_BOUNDTYPE        boundtype,          /**< Whether we have an upper bound or a lower bound */
   SCIP_RATIONAL*        value,              /**< value of the bound */
   SCIP_Longint          certificateindex    /**< index in the certificate */
   )
{
   SCIP_RETCODE res;

   if( !SCIPisCertified(scip) )
      return SCIP_OKAY;

   assert(certificate != NULL);

   switch( var->varstatus )
   {
      case SCIP_VARSTATUS_FIXED:
      case SCIP_VARSTATUS_LOOSE:
      case SCIP_VARSTATUS_ORIGINAL:
         SCIPABORT();
         return SCIP_ERROR;
      case SCIP_VARSTATUS_NEGATED:
         SCIPrationalMultReal(value, value, -1);
         assert( SCIPvarGetNegationConstant(var) == 1 );
         SCIPrationalAddReal(value, value, 1.0);
         res = SCIPcertificatePrintGlobalBound(scip, certificate, var->negatedvar, boundtype == SCIP_BOUNDTYPE_UPPER ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER, value, certificateindex);
         SCIPrationalAddReal(value, value, -1.0);
         SCIPrationalMultReal(value, value, -1);
         return res;
         break;
      case SCIP_VARSTATUS_AGGREGATED:
         SCIPrationalDiv(value, value, var->exactdata->aggregate.scalar);
         assert(SCIPrationalIsZero(var->exactdata->aggregate.constant));
         res = SCIPcertificatePrintGlobalBound(scip, certificate, var->data.aggregate.var, (boundtype == SCIP_BOUNDTYPE_UPPER) == SCIPrationalIsPositive(var->exactdata->aggregate.scalar) ? SCIP_BOUNDTYPE_UPPER: SCIP_BOUNDTYPE_LOWER, value, certificateindex);
         SCIPrationalMult(value, value, var->exactdata->aggregate.scalar);
         return res;
         break;
      case SCIP_VARSTATUS_COLUMN:
         break;
      default:
         SCIPABORT();
         return SCIP_ERROR;
   }

#ifndef NDEBUG
   certificate->lastinfo->isbound = TRUE;
   certificate->lastinfo->boundtype = boundtype;
   certificate->lastinfo->varindex = SCIPvarGetCertificateIndex(var);
   certificate->lastinfo->isglobal = TRUE;
   certificate->lastinfo->certificateindex = certificate->indexcounter;
   SCIPrationalSetRational(certificate->lastinfo->boundval, value);
#endif

   SCIPcertificatePrintProofMessage(certificate, "GlobalBound_%d %c ", certificate->indexcounter,
      boundtype == SCIP_BOUNDTYPE_UPPER ? 'L' : 'G');
   SCIP_CALL( SCIPcertificatePrintProofRational(certificate, value) );
   SCIPcertificatePrintProofMessage(certificate, " 1 %d 1 ", SCIPvarGetCertificateIndex(var));
   SCIPcertificatePrintProofMessage(certificate, "{ lin 1 %d 1 } -1 global\n", certificateindex);

   certificate->indexcounter++;
   return SCIP_OKAY;
}

/* prints information for constraint to certificate file */
SCIP_RETCODE SCIPconsPrintCertificateExactLinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CERTIFICATE* certificate;
   SCIP_ROWEXACT* row;
   SCIP_RATIONAL* correctedside;
   int* varsindex = NULL;
   int i;
   SCIP_Longint image;
   SCIP_RATIONAL* lhs;
   SCIP_RATIONAL* rhs;

   /*lint --e{715}*/
   assert(scip != NULL);
   assert(cons != NULL);

   /* print constraint into certificate output */
   if( !SCIPisCertified(scip) )
      return SCIP_OKAY;
   certificate = SCIPgetCertificate(scip);
   row = SCIPgetRowExactExactLinear(scip, cons);

   lhs = SCIProwExactGetLhs(row);
   rhs = SCIProwExactGetRhs(row);

   assert(row != NULL);

   image = SCIPhashmapGetImageLong(certificate->rowdatahash, row);
   /* add row to hashmap */
   if( image != SCIP_LONGINT_MAX )
   {
      SCIPmessageFPrintWarning(scip->messagehdlr, "%lu \n", (size_t) SCIPhashmapGetImage(certificate->rowdatahash, row));
      SCIPerrorMessage("Duplicate row in certificate row hashmap\n");
      SCIPABORT();
      return SCIP_ERROR;
   }
   else
   {
      SCIP_CALL( SCIPhashmapInsertLong(certificate->rowdatahash, row, certificate->indexcounter) );
      SCIP_CALL( SCIPhashmapInsertLong(certificate->rowdatahash, cons, certificate->indexcounter) );
      assert(SCIPhashmapExists(certificate->rowdatahash, row));
      assert(SCIPhashmapExists(certificate->rowdatahash, cons));
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varsindex, SCIProwExactGetNNonz(row)) );
   for( i = 0; i < SCIProwExactGetNNonz(row); ++i )
      varsindex[i] = SCIPvarGetCertificateIndex(SCIPcolExactGetVar(SCIProwExactGetCols(row)[i]));

   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &correctedside) );

   /* print constraint */
   if( SCIPrationalIsEQ(lhs, rhs) )
   {
      assert(!SCIPrationalIsAbsInfinity(lhs));
      SCIPrationalDiff(correctedside, SCIProwExactGetLhs(row), SCIProwExactGetConstant(row));
      SCIP_CALL( SCIPcertificatePrintCons(certificate, FALSE, NULL, 'E', correctedside, SCIProwExactGetNNonz(row), varsindex, SCIProwExactGetVals(row)) );
   }
   else
   {
      if( !SCIPrationalIsNegInfinity(lhs) )
      {
         SCIPrationalDiff(correctedside, SCIProwExactGetLhs(row), SCIProwExactGetConstant(row));
         SCIP_CALL( SCIPcertificatePrintCons(certificate, FALSE, NULL, 'G', correctedside, SCIProwExactGetNNonz(row), varsindex, SCIProwExactGetVals(row)) );
      }
      if( !SCIPrationalIsInfinity(rhs) )
      {
         SCIPrationalDiff(correctedside, SCIProwExactGetRhs(row), SCIProwExactGetConstant(row));
         SCIP_CALL( SCIPcertificatePrintCons(certificate, FALSE, NULL, 'L', correctedside, SCIProwExactGetNNonz(row), varsindex, SCIProwExactGetVals(row)) );
      }
   }

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &correctedside);
   SCIPfreeBufferArray(scip, &varsindex);

   return SCIP_OKAY;
}

/** returns the index of the given constraint in the certificate */
SCIP_Longint SCIPcertificateGetConsIndex(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_RATIONAL*        lhs,                /**< lhs of the constraint */
   SCIP_RATIONAL*        rhs,                /**< rhs of the constraint */
   SCIP_Bool             useRhs              /**< whether to return the index of the rhs or lhs */
   )
{
   SCIP_Longint ret;

   ret = SCIPhashmapGetImageLong(certificate->rowdatahash, cons);

   assert( ret != SCIP_LONGINT_MAX );
   if( !SCIPrationalIsAbsInfinity(rhs) && !SCIPrationalIsAbsInfinity(lhs) && !SCIPrationalIsEQ(lhs, rhs) && useRhs)
      ret += 1;

   assert(ret >= 0);

   return ret;
}
