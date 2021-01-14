/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   certificate.c
 * @brief  methods for certificate output
 * @author Ambros Gleixner
 * @author Daniel Steffy
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "blockmemshell/memory.h"
#include "scip/cons_exactlp.h"
#include "scip/scip.h"
#include "scip/set.h"
#include "scip/lp.h"
#include "scip/lpexact.h"
#include "scip/pub_lpexact.h"
#include "lpiexact/lpiexact.h"
#include "scip/pub_misc.h"
#include "scip/prob.h"
#include "scip/certificate.h"
#include "scip/struct_certificate.h"
#include "scip/sol.h"
#include "scip/struct_scip.h"
#include "scip/var.h"

#define SCIP_HASHSIZE_CERTIFICATE    500 /**< size of hash map for certificate -> nodesdata mapping used for certificate output */

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyVarbound)
{  /*lint --e{715}*/
   /* the key is the element itself */
   return elem;
}

/** returns TRUE iff both keys are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqVarbound)
{
   assert(key1 != NULL);
   assert(key2 != NULL);

   if( ((SCIP_CERTIFICATEBOUND*)key1)->isupper != ((SCIP_CERTIFICATEBOUND*)key2)->isupper )
      return FALSE;

   if( ((SCIP_CERTIFICATEBOUND*)key1)->varindex != ((SCIP_CERTIFICATEBOUND*)key2)->varindex )
      return FALSE;

   if( !RatIsEqual(((SCIP_CERTIFICATEBOUND*)key1)->boundval, ((SCIP_CERTIFICATEBOUND*)key2)->boundval) )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValVarbound)
{
   SCIP_CERTIFICATEBOUND* bound;

   bound = (SCIP_CERTIFICATEBOUND*)key;

   return ((((unsigned int)(10*RatApproxReal(bound->boundval))) << 22) + (bound->varindex << 2) + (unsigned int)bound->isupper);
}

/** updates file size */
static
void updateFilesize(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Real             nchars              /**< number of characters printed */
   )
{
   certificate->filesize += nchars/1048576.0;
}

static
int getVarBoundFileIndex(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_COLEXACT*        col,                /**< exact column */
   SCIP_Rational*        multiplier          /**< dual multiplier (farkascoef or redcost) */
   )
{
   void* image;

   if( RatIsNegative(multiplier) )
   {
      certificate->workbound->isupper = TRUE;
      RatSet(certificate->workbound->boundval, SCIPcolExactGetUb(col));
   }
   else
   {
      certificate->workbound->isupper = FALSE;
      RatSet(certificate->workbound->boundval, SCIPcolExactGetLb(col));
   }
   certificate->workbound->varindex = SCIPvarGetCertificateIndex(SCIPcolGetVar(col->fpcol));

   image = SCIPhashtableRetrieve(certificate->varboundtable, certificate->workbound);

   if( image == NULL )
   {
      SCIPerrorMessage("col not in varbound-table \n");
      SCIPABORT();
      return -1;
   }

   return ((SCIP_CERTIFICATEBOUND*)image)->fileindex;
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

/** print variable bound assumption into certificate if not already present and add it to varboundtable,
 *  return index of this bound in the certificate file
 */
static
SCIP_Longint printBoundAssumption(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   int                   varindex,           /**< index of the variable */
   SCIP_Rational*        boundval,           /**< value of the bound */
   SCIP_BOUNDTYPE        boundtype           /**< is it the upper bound? */
   )
{
   SCIP_CERTIFICATEBOUND* image;
   SCIP_CERTIFICATEBOUND* foundbound;

   /* check whether output should be created */
   if ( certificate->origfile == NULL )
      return SCIP_OKAY;

   /* install bound information in working structure */
   certificate->workbound->fileindex = certificate->indexcounter;
   certificate->workbound->varindex = varindex;
   RatSet(certificate->workbound->boundval, boundval);
   certificate->workbound->isupper = boundtype == SCIP_BOUNDTYPE_UPPER ? TRUE : FALSE;

   /* check if an assumption with this current working structure exists, NULL if not.
    * If the assumption already exists, then we only check that it has the good values. If the assumption does not exists,
    * then we insert the new assumption into the hashtable storing all variable bound assumptions and finally print it
    */
   image = SCIPhashtableRetrieve(certificate->varboundtable, certificate->workbound);

   if( image != NULL )
   {
      foundbound = (SCIP_CERTIFICATEBOUND*)image;

      SCIPdebugMessage("Found bound assumption at line %" SCIP_LONGINT_FORMAT ": <variable %d> %s approx. %g\n",
         foundbound->fileindex, varindex, (certificate->workbound->isupper ? "<=" : ">="), RatApproxReal(boundval));

      assert(foundbound->fileindex >= 0);
      assert(foundbound->varindex == varindex);
      assert(RatIsEqual(foundbound->boundval, boundval));
      assert(foundbound->isupper == certificate->workbound->isupper);

      return foundbound->fileindex;
   }
   else
   {
      SCIP_CERTIFICATEBOUND* insertbound;

      SCIPdebugMessage("Print bound assumption at line %" SCIP_LONGINT_FORMAT ": <variable %d> %s approx. %g\n",
         certificate->workbound->fileindex, varindex, (certificate->workbound->isupper  ? "<=" : ">="), RatApproxReal(boundval));

      /* insertion */
      SCIP_ALLOC( BMSduplicateBlockMemory(certificate->blkmem, &insertbound, certificate->workbound) );
      SCIP_CALL( RatCreateBlock(certificate->blkmem, &insertbound->boundval) );
      RatSet(insertbound->boundval, boundval);

      /* ensure size and insert boundval in array to be able to free it at the end */
      if( certificate->nboundvals >= certificate->boundvalsize )
      {
         BMSreallocBlockMemoryArray(certificate->blkmem, &certificate->boundvals, certificate->boundvalsize, certificate->boundvalsize + 100);
         certificate->boundvalsize += 100;
      }
      certificate->boundvals[certificate->nboundvals] = insertbound;
      certificate->nboundvals++;

      SCIP_CALL( SCIPhashtableInsert(certificate->varboundtable, insertbound) );

      /** @todo: it could be better to seperate the printing from insertion of variable bound */
      SCIPcertificatePrintProofMessage(certificate, "A%lld %c ", insertbound->fileindex, (certificate->workbound->isupper ? 'L' : 'G'));

      SCIPcertificatePrintProofRational(certificate, boundval, 10);
      SCIPcertificatePrintProofMessage(certificate, " 1 %d 1 { asm } -1\n", varindex);
      certificate->indexcounter++;

      return insertbound->fileindex;
   }
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

   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);
   RatFreeBlock(certificate->blkmem, &nodedata->derbound_left);
   RatFreeBlock(certificate->blkmem, &nodedata->derbound_right);
   RatFreeBlock(certificate->blkmem, &nodedata->derbound_inherit);
   BMSfreeBlockMemory(certificate->blkmem, &nodedata);
   SCIP_CALL( SCIPhashmapRemove(certificate->nodedatahash, node) );

   return SCIP_OKAY;
}

/** print the best solution found */
static
SCIP_RETCODE SCIPcertificatePrintSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_SOL*             sol                 /**< solution to be printed */
   )
{
   SCIP_VAR** vars;
   SCIP_Rational** vals;
   int nvars;
   int nnonz;
   int i;
   int j;

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

   SCIP_CALL( RatCreateBufferArray(SCIPbuffer(scip), &vals, nvars) );

   /* get number of non-zero coefficient in the solution */
   nnonz = 0;
   for( i = 0; i < nvars; i++)
   {
      SCIPsolGetValExact(vals[i], sol, scip->set, scip->stat, vars[i]);
      if( !RatIsZero(vals[i]) )
         nnonz++;
   }

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "SOL 1\nbest %d", nnonz);

   for( i = 0; i < nvars; i++ )
   {
      if( !RatIsZero(vals[i]) )
      {
         /* print the solution into certificate */
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, " %d ", SCIPvarGetCertificateIndex(vars[i]));
         SCIPcertificatePrintProblemRational(certificate, isorigfile, vals[i], 10);
      }
   }
   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "\n");

   RatFreeBufferArray(SCIPbuffer(scip), &vals, nvars);

   return SCIP_OKAY;
}

/** set the node to have its own bound proof */
SCIP_RETCODE SCIPcertificateSetInheritanceData(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_NODE*            node,               /**< node data structure */
   SCIP_Longint          fileindex,          /**< index of new bound */
   SCIP_Rational*        newbound            /**< the inherited bound */
   )
{
   SCIP_CERTNODEDATA* nodedata;

   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);

   /* do nothing if the newbound is worse than the inherited bound */
   if( RatIsLT(newbound, nodedata->derbound_inherit) )
      return SCIP_OKAY;

   nodedata->inheritedbound = FALSE;
   nodedata->derindex_inherit = fileindex;
   RatSet(nodedata->derbound_inherit, newbound);

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
   (*certificate)->varboundtable = NULL;
   (*certificate)->workbound = NULL;
   (*certificate)->blkmem = NULL;
   (*certificate)->indexcounter = 0;
   (*certificate)->conscounter = 0;
   (*certificate)->origfile = NULL;
   (*certificate)->transfile = NULL;
   (*certificate)->origfilename = NULL;
   (*certificate)->derivationfile = NULL;
   (*certificate)->derivationfilename = NULL;
   (*certificate)->filesize = 0.0;
   (*certificate)->rowdatahash = NULL;
   (*certificate)->boundvals = NULL;
   (*certificate)->boundvalsize = 0;
   (*certificate)->nboundvals = 0;
   (*certificate)->nodedatahash = NULL;
   (*certificate)->rootbound = NULL;
   (*certificate)->finalbound = NULL;
   (*certificate)->derindex_root = -1;
   (*certificate)->rootinfeas = FALSE;
   (*certificate)->objintegral = FALSE;
   (*certificate)->vals = NULL;
   (*certificate)->valssize = 0;

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
   int j;
   char* name = NULL;
   char* compression = NULL;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_Rational* lb;
   SCIP_Rational* ub;

   assert(certificate != NULL);
   assert(set != NULL);
   assert(set->certificate_filename != NULL);
   assert(certificate->derivationfile == NULL);
   assert(certificate->nodedatahash == NULL);
   assert(certificate->rowdatahash == NULL);

   if( !(set->exact_enabled) || (set->certificate_filename[0] == '-' && set->certificate_filename[1] == '\0') )
      return SCIP_OKAY;

   filenamelen = strlen(set->certificate_filename);
   SCIP_CALL( SCIPsetAllocBufferArray(set, &name, filenamelen + 1) );
   BMScopyMemoryArray(name, set->certificate_filename, filenamelen);
   name[filenamelen] = '\0';

   /** @todo exip: this currently strips the .vipr from the filename if there is no compression...
    * the whole compression code currently makes no sense, since the certificate has to be unpacked to be read by vipr
    * again anyway, so it is just disabled.
    */
   //SCIPsplitFilename(name, NULL, NULL, NULL, &compression);

   SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      "storing certificate information in file <%s>\n", set->certificate_filename);

   if( NULL != compression && 0 == strncmp(compression, "gz", 2) )
      certificate->transfile = SCIPfopen(set->certificate_filename, "wb");
   else
      certificate->transfile = SCIPfopen(set->certificate_filename, "wT");

   bufferlen = strlen(name);
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
   if( NULL != compression && 0 == strncmp(compression, "gz", 2) )
   {
      certificate->derivationfilename[bufferlen+4] = '.';
      certificate->derivationfilename[bufferlen+5] = 'g';
      certificate->derivationfilename[bufferlen+6] = 'z';
      certificate->derivationfilename[bufferlen+7] = '\0';

      certificate->origfilename[bufferlen+4] = '.';
      certificate->origfilename[bufferlen+5] = 'g';
      certificate->origfilename[bufferlen+6] = 'z';
      certificate->origfilename[bufferlen+7] = '\0';

      certificate->derivationfile = SCIPfopen(certificate->derivationfilename, "wb");
      certificate->origfile = SCIPfopen(certificate->origfilename, "wb");
   }
   else
   {
      certificate->derivationfilename[bufferlen+4] = '\0';
      certificate->origfilename[bufferlen+4] = '\0';
      certificate->derivationfile = SCIPfopen(certificate->derivationfilename, "wT");
      certificate->origfile = SCIPfopen(certificate->origfilename, "wT");
   }

   if( certificate->origfile == NULL || certificate->derivationfile == NULL )
   {
      SCIPerrorMessage("error creating file <%s> and derivation file\n", set->certificate_filename);
      SCIPprintSysError(set->certificate_filename);
      return SCIP_FILECREATEERROR;
   }

   /* initialisation of hashmaps and hashtables */
   SCIP_CALL( SCIPhashmapCreate(&certificate->nodedatahash, blkmem, SCIP_HASHSIZE_CERTIFICATE) );
   SCIP_CALL( SCIPhashmapCreate(&certificate->rowdatahash, blkmem, SCIP_HASHSIZE_CERTIFICATE) );

   certificate->blkmem = blkmem;
   SCIPsetFreeBufferArray(set, &name);

   SCIP_CALL( SCIPgetOrigVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   nboundconss = 0;
   for ( j = 0 ; j < nvars ; j++ )
   {
      lb = SCIPvarGetLbGlobalExact(vars[j]);
      ub = SCIPvarGetUbGlobalExact(vars[j]);
      if( !RatIsAbsInfinity(lb) )
         nboundconss++;
      if( !RatIsAbsInfinity(ub) )
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
      SCIPvarSetCertificateIndex(SCIPvarGetTransVar(vars[j]), j);
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
         SCIPcertificatePrintProblemMessage(certificate, TRUE, "%d \n", SCIPvarGetIndex(vars[j]));
      }
   }

   {
      SCIP_Rational** objcoefs;
      SCIP_CALL( RatCreateBufferArray(SCIPbuffer(scip), &objcoefs, nvars) );

      for( j = 0; j < nvars; j++)
         RatSet(objcoefs[j], SCIPvarGetObjExact(vars[j]));

      /* print the objective function into certificate header */
      SCIPcertificateSetAndPrintObjective(certificate, TRUE, blkmem, objcoefs, nvars);

      RatFreeBufferArray(SCIPbuffer(scip), &objcoefs, nvars);
   }

   conss = SCIPgetOrigConss(scip);
   ncertcons = 0;
   for( j = 0; j < SCIPgetNOrigConss(scip); j++ )
   {
      SCIP_CONS* cons;
      SCIP_CONSHDLR* conshdlr;

      cons = conss[j];
      conshdlr = SCIPconsGetHdlr(cons);

      if( strcmp(SCIPconshdlrGetName(conshdlr), "linear-exact") == 0 )
      {
         lb = SCIPgetLhsExactLinear(scip, cons);
         ub = SCIPgetRhsExactLinear(scip, cons);
         /* if the constraint a singleton we count it as a bound constraint in the certificate */
         if( SCIPgetNVarsExactLinear(scip, cons) == 1 )
         {
            if( !RatIsAbsInfinity(lb) )
               nboundconss++;
            if( !RatIsAbsInfinity(ub) )
               nboundconss++;
         }
         /* else we check if it is a ranged row with 2 sides (then we count it twice). otherwise one line gets printed */
         else if( !RatIsEqual(lb, ub) && !RatIsAbsInfinity(lb) && !RatIsAbsInfinity(ub) )
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

   for( j = 0; j < nvars; j++ )
   {
      if( !RatIsAbsInfinity(SCIPvarGetLbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, TRUE, NULL, SCIPvarGetCertificateIndex(vars[j]), SCIPvarGetLbGlobalExact(vars[j]), FALSE) );
      }
      if( !RatIsAbsInfinity(SCIPvarGetUbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, TRUE, NULL, SCIPvarGetCertificateIndex(vars[j]), SCIPvarGetUbGlobalExact(vars[j]), TRUE) );
      }
   }

   SCIP_CALL( RatCreateBlock(blkmem, &certificate->rootbound) );
   SCIP_CALL( RatCreateBlock(blkmem, &certificate->finalbound) );
   certificate->valssize = SCIPgetNVars(scip) + SCIPgetNConss(scip);
   SCIP_CALL( RatCreateBlockArray(SCIPblkmem(scip), &(certificate->vals), certificate->valssize) );

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
   SCIP_Rational* lb;
   SCIP_Rational* ub;
   SCIP_CERTIFICATE* certificate;
   BMS_BLKMEM* blkmem;

   assert(scip != NULL);
   assert(scip->set->certificate_filename != NULL);

   certificate = SCIPgetCertificate(scip);
   blkmem = SCIPblkmem(scip);

   assert(certificate != NULL);

   if( !(scip->set->exact_enabled) || (scip->set->certificate_filename[0] == '-' && scip->set->certificate_filename[1] == '\0') )
      return SCIP_OKAY;

   if( certificate->origfile == NULL || certificate->derivationfile == NULL )
   {
      SCIPerrorMessage("error creating file <%s> and derivation file\n", scip->set->certificate_filename);
      SCIPprintSysError(scip->set->certificate_filename);
      return SCIP_FILECREATEERROR;
   }

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   nboundconss = 0;
   for ( j = 0 ; j < nvars ; j++ )
   {
      lb = SCIPvarGetLbGlobalExact(vars[j]);
      ub = SCIPvarGetUbGlobalExact(vars[j]);
      if( !RatIsAbsInfinity(lb) )
         nboundconss++;
      if( !RatIsAbsInfinity(ub) )
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
      SCIP_Rational** objcoefs;
      SCIP_CALL( RatCreateBufferArray(SCIPbuffer(scip), &objcoefs, nvars) );

      for( j = 0; j < nvars; j++)
         RatSet(objcoefs[j], SCIPvarGetObjExact(vars[j]));

      /* print the objective function into certificate header */
      SCIPcertificateSetAndPrintObjective(certificate, FALSE, blkmem, objcoefs, nvars);

      RatFreeBufferArray(SCIPbuffer(scip), &objcoefs, nvars);
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

      if( strcmp(SCIPconshdlrGetName(conshdlr), "linear-exact") == 0 )
      {
         lb = SCIPgetLhsExactLinear(scip, cons);
         ub = SCIPgetRhsExactLinear(scip, cons);
         /* if the constraint a singleton we count it as a bound constraint in the certificate */
         if( SCIPgetNVarsExactLinear(scip, cons) == 1 )
         {
            SCIPdebugMessage("constraint counts as bound constraint \n");

            if( !RatIsAbsInfinity(lb) )
               nboundconss++;
            if( !RatIsAbsInfinity(ub) )
               nboundconss++;
         }
         /* else we check if it is a ranged row with 2 sides (then we count it twice). otherwise one line gets printed */
         else if( !RatIsEqual(lb, ub) && !RatIsAbsInfinity(lb) && !RatIsAbsInfinity(ub) )
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

   for( j = 0; j < nvars; j++ )
   {
      if( !RatIsAbsInfinity(SCIPvarGetLbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, FALSE, NULL, SCIPvarGetCertificateIndex(vars[j]), SCIPvarGetLbGlobalExact(vars[j]), FALSE) );
      }
      if( !RatIsAbsInfinity(SCIPvarGetUbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, FALSE, NULL, SCIPvarGetCertificateIndex(vars[j]), SCIPvarGetUbGlobalExact(vars[j]), TRUE) );
      }
   }

   return SCIP_OKAY;
}

/** Concatenate the certificate and the _der file and delete the _der file  */
static
void concatCert(
   SCIP_CERTIFICATE*     certificate,        /**< The certificate pointer */
   const char*           certificatefilename /**< The name of the certificate file */
   )
{
   SCIP_FILE* derivationfile;
   char buffer[SCIP_MAXSTRLEN];
   size_t size;

   derivationfile = SCIPfopen(certificate->derivationfilename, "r");

   /* append the derivation file to the problem file */
   while( 0 != (size = SCIPfread(buffer, sizeof(char), SCIP_MAXSTRLEN, derivationfile)) )
      SCIPfwrite(buffer, sizeof(char), size, certificate->transfile);

   SCIPfclose(derivationfile);
   /* delete the derivation file */
   remove(certificate->derivationfilename);
}

/** closes the certificate output files */
void SCIPcertificateExit(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   int i;

   assert(certificate != NULL);
   assert(set != NULL);

   if( certificate->origfile != NULL )
   {
      SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
         "closing CERTIFICATE information file (wrote %.1f MB)\n", certificate->filesize);

      if( certificate->derivationfile != NULL )
      {
         /* CERT TODO: DER line with counter and append two files */
         SCIPfclose(certificate->derivationfile);
         certificate->derivationfile = NULL;
         concatCert(certificate, set->certificate_filename);
      }
      SCIPfclose(certificate->origfile);
      SCIPfclose(certificate->transfile);
      certificate->origfile = NULL;
      certificate->transfile = NULL;

      BMSfreeMemoryArray(&certificate->derivationfilename);
      BMSfreeMemoryArray(&certificate->origfilename);
      if( certificate->varboundtable != NULL )
      {
         for( i = 0; i < certificate->nboundvals; i++ )
         {
            RatFreeBlock(certificate->blkmem, &certificate->boundvals[i]->boundval);
            BMSfreeBlockMemory(certificate->blkmem, &certificate->boundvals[i]);
         }
         /**@todo fix memory leak: mpq_clear and free all elements */
         SCIPhashtableRemoveAll(certificate->varboundtable);
         SCIPhashtableFree(&certificate->varboundtable);
         BMSfreeBlockMemoryArray(certificate->blkmem, &certificate->boundvals, certificate->boundvalsize);
      }
      if( certificate->workbound != NULL )
      {
         RatFreeBlock(certificate->blkmem, &certificate->workbound->boundval);
         BMSfreeBlockMemory(certificate->blkmem, &certificate->workbound);
      }

      if( certificate->rowdatahash)
         SCIPhashmapFree(&certificate->rowdatahash);

   if( certificate->nodedatahash )
      {
         assert(SCIPhashmapIsEmpty(certificate->nodedatahash));
         SCIPhashmapFree(&certificate->nodedatahash);
      }

   RatFreeBlock(certificate->blkmem, &certificate->rootbound);
   RatFreeBlock(certificate->blkmem, &certificate->finalbound);
   RatFreeBlockArray(certificate->blkmem, &certificate->vals, certificate->valssize);
   }
}

/** returns whether the certificate output is activated? */
SCIP_Bool SCIPcertificateIsActive(
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   return (set->exact_enabled && certificate != NULL && (certificate->transfile != NULL || strcmp(set->certificate_filename, "-")) );
}

/** returns current certificate file size in MB */
SCIP_Real SCIPcertificateGetFilesize(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   if( certificate == NULL || certificate->transfile == NULL)
      return 0.0;
   else
      return certificate->filesize;
}

/** sets the objective function used when printing dual bounds */
SCIP_RETCODE SCIPcertificateSetAndPrintObjective(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Rational**       coefs,              /**< objective function coefficients */
   int                   nvars               /**< number of variables */
   )
{
   char obj[SCIP_MAXSTRLEN - 2];
   int nnonz;
   int i;

   assert(coefs != NULL);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   /* create a hash table for the variable bounds (we work on the tranformed problem) */
   if( isorigfile )
   {
      SCIP_CALL( SCIPhashtableCreate(&certificate->varboundtable, blkmem, nvars,
            hashGetKeyVarbound, hashKeyEqVarbound, hashKeyValVarbound, NULL) );

      /* create working memory for bound struct */
      SCIP_ALLOC( BMSallocBlockMemory(blkmem, &certificate->workbound) );
      SCIP_CALL( RatCreateBlock(blkmem, &certificate->workbound->boundval) );
   }

   nnonz = 0;
   for( i = 0; i < nvars; i++ )
   {
      if( !RatIsZero(coefs[i]) )
         nnonz++;
   }

   SCIPsnprintf(obj, SCIP_MAXSTRLEN - 2, "OBJ min\n %d ", nnonz);

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s ", obj);

   for( i = 0; i < nvars; i++ )
   {
      if( !RatIsZero(coefs[i]) )
      {
         SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s%d ", (i > 0 ? " " : ""), i );
         SCIPcertificatePrintProblemRational(certificate, isorigfile, coefs[i], 10);
      }
   }

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "\n");

   return SCIP_OKAY;
}

/** prints the last part of the certificate header (RTP range/sol, ...) */
SCIP_RETCODE SCIPcertificatePrintResult(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   SCIP_Rational* primalbound;
   SCIP_Rational* dualbound;
   SCIP_SOL* bestsol;

   assert(scip != NULL);

   if( certificate->transfile == NULL)
      return SCIP_OKAY;

   SCIP_CALL( RatCreateBuffer(set->buffer, &primalbound) );
   SCIP_CALL( RatCreateBuffer(set->buffer, &dualbound) );

   /* check status first: OPTIMAL / INFEAS / NOT OPTIMAL */
   if( SCIPisInRestart(scip) )
      return SCIP_ERROR;

   if( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL )
   {
      bestsol = SCIPgetBestSol(scip);

      if( !SCIPisExactSol(scip, bestsol) )
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

      assert(!RatIsAbsInfinity(primalbound));

      /* print RTP range (same when optimal solution found) */
      SCIPcertificatePrintRtpRange(certificate, isorigfile, primalbound, primalbound);

      /* print optimal solution into certificate */
      SCIP_CALL( SCIPcertificatePrintSol(scip, isorigfile, certificate, bestsol) );
   }
   else if( SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE )
   {
      SCIPcertificatePrintRtpInfeas(certificate, isorigfile);
      SCIP_CALL(SCIPcertificatePrintSol(scip, isorigfile, certificate, NULL) );
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
         RatSetString(primalbound, "inf");
      }

      SCIPcertificatePrintRtpRange(certificate, isorigfile, certificate->finalbound, primalbound);
      SCIP_CALL( SCIPcertificatePrintSol(scip, isorigfile, certificate, bestsol) );
   }
   if( !isorigfile )
      SCIPcertificatePrintDerHeader(certificate);

   RatFreeBuffer(set->buffer, &dualbound);
   RatFreeBuffer(set->buffer, &primalbound);

   return SCIP_OKAY;
}

/** prints the last part of the certificate header (RTP range/sol, ...) */
SCIP_RETCODE SCIPcertificateSaveFinalbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   assert(scip != NULL);

   if( certificate->transfile == NULL)
      return SCIP_OKAY;

   SCIPgetLowerboundExact(scip, certificate->finalbound);

   return SCIP_OKAY;
}

/** prints a string to the problem section of the certificate file */
void SCIPcertificatePrintProblemMessage(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[3 * SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   va_start(ap, formatstr);
   vsnprintf(buffer, 3 * SCIP_MAXSTRLEN, formatstr, ap);

   if( isorigfile )
      SCIPfprintf(certificate->origfile, "%s", buffer);
   else
      SCIPfprintf(certificate->transfile, "%s", buffer);

   va_end(ap);
   updateFilesize(certificate, strlen(buffer));
}

/** prints a string to the proof section of the certificate file */
void SCIPcertificatePrintProofMessage(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->derivationfile == NULL )
      return;
   va_start(ap, formatstr);
   vsprintf(buffer, formatstr, ap);

   SCIPfprintf(certificate->derivationfile, "%s", buffer); // todo: is this correct?
   va_end(ap);
   updateFilesize(certificate, strlen(buffer));
}


/** prints a rational number to the problem section of the certificate file */
void SCIPcertificatePrintProblemRational(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   SCIP_Rational*        val,                /**< Rational to print to the problem*/
   int                   base                /**< The base representation*/
   )
{
   SCIP_Longint len = RatStrlen(val) + 1;
   char* formatstr;

   assert(len <= INT_MAX);

   /* check if certificate output should be created */
   if( certificate->derivationfile == NULL )
     return;

   BMSallocMemoryArray(&formatstr, len);
   RatToString(val, formatstr, len);
   if( isorigfile )
      SCIPfprintf(certificate->origfile, "%s", formatstr);
   else
      SCIPfprintf(certificate->transfile, "%s", formatstr);

   BMSfreeMemoryArray(&formatstr);
}


/** prints a rational number to the proof section of the certificate file */
void SCIPcertificatePrintProofRational(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Rational*        val,                /**< Rational to print to the problem*/
   int                   base                /**< The base representation*/
   )
{
   SCIP_Longint len = RatStrlen(val) + 1;
   char* formatstr;

   assert(len <= INT_MAX);

   /* check if certificate output should be created */
   if( certificate->derivationfile == NULL )
     return;

   BMSallocMemoryArray(&formatstr, len);
   RatToString(val, formatstr, len);
   SCIPfprintf(certificate->derivationfile, "%s", formatstr);
   BMSfreeMemoryArray(&formatstr);
}

/** prints a comment to the problem section of the certificate file */
void SCIPcertificatePrintProblemComment(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPfprintf(certificate->origfile, "# ");

   va_start(ap, formatstr);
   vsprintf(buffer, formatstr, ap);

   if( isorigfile )
      SCIPfprintf(certificate->origfile, "%s", formatstr); // todo: is this correct?
   else
      SCIPfprintf(certificate->transfile, "%s", formatstr); // todo: is this correct?

   va_end(ap);
   updateFilesize(certificate, 2 + strlen(formatstr));
}

/** prints a comment to the proof section of the certificate file */
void SCIPcertificatePrintProofComment(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->derivationfile == NULL )
      return;

   SCIPfprintf(certificate->derivationfile, "# ");

   va_start(ap, formatstr);
   vsprintf(buffer, formatstr, ap);

   SCIPfprintf(certificate->derivationfile, "%s", formatstr); // todo: is this correct?
   va_end(ap);
   updateFilesize(certificate, 2 + strlen(formatstr));
}

/** prints version header */
void SCIPcertificatePrintVersionHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile          /**< shoud the line be printed to the origfile or the transfile */
   )
{
   assert(certificate != NULL);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "VER 1.0 \n");
}

/** prints variable section header */
void SCIPcertificatePrintVarHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   int                   nvars               /**< number of variables */
   )
{
   assert(certificate != NULL);
   assert(nvars >= 0);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "VAR %d \n", nvars);
}

/** prints integer section header */
void SCIPcertificatePrintIntHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   int                   nints               /**< number of integer variables */
   )
{
   assert(certificate != NULL);
   assert(nints >= 0);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "INT %d\n", nints);
}

/** prints constraint section header */
void SCIPcertificatePrintConsHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   int                   nconss,             /**< number of all constraints */
   int                   nboundconss         /**< number of bound constraints */
   )
{
   assert(certificate != NULL);
   assert(nconss >= 0);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "CON %d %d\n", nconss + nboundconss, nboundconss);
}

/** prints derivation section header */
void SCIPcertificatePrintDerHeader(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   int nders;
   assert(certificate != NULL);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   nders = certificate->indexcounter - certificate->conscounter;
   SCIPcertificatePrintProblemMessage(certificate, FALSE, "DER %d\n", nders);
}

/** prints constraint */
void SCIPcertificatePrintCons(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   const char*           consname,           /**< name of the constraint */
   const char            sense,              /**< sense of the constraint, i.e., G, L, or E */
   SCIP_Rational*        side,               /**< left/right-hand side */
   int                   len,                /**< number of nonzeros */
   int*                  ind,                /**< index array */
   SCIP_Rational**       val                 /**< coefficient array */
   )
{
   int i;

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   if( consname == NULL )
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "C%d %c ", certificate->indexcounter, sense);
   else
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s %c ", consname, sense);

   SCIPcertificatePrintProblemRational(certificate, isorigfile, side, 10);

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, " %d", len);

   for( i = 0; i < len; i++ )
   {
      /** @todo exip: perform line breaking before exceeding maximum line length */
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, " %d ", ind[i]);
      SCIPcertificatePrintProblemRational(certificate, isorigfile, val[i], 10);
   }
   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "\n");

   certificate->indexcounter++;

   if( !isorigfile )
   {
      certificate->conscounter++;
   }
}

/** @todo exip: refactor this so duplicates and redundant bounds do not get printed */
/** prints a variable bound to the problem section of the certificate file and returns line index */
SCIP_RETCODE SCIPcertificatePrintBoundCons(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Bool             isorigfile,         /**< shoud the line be printed to the origfile or the transfile */
   const char*           boundname,          /**< name of the bound constraint */
   int                   varindex,           /**< index of the variable */
   SCIP_Rational*        boundval,           /**< value of the bound */
   SCIP_Bool             isupper             /**< is it the upper bound? */
   )
{
   void* image;

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   if( !isorigfile )
   {
      /* install bound information in working struct */
      certificate->workbound->fileindex = certificate->indexcounter;
      certificate->workbound->varindex = varindex;
      RatSet(certificate->workbound->boundval, boundval);
      certificate->workbound->isupper = isupper;

      SCIPdebugMessage("Printing bound at line %" SCIP_LONGINT_FORMAT ": <variable %d> %s approx. %g\n",
         certificate->indexcounter, varindex, (isupper ? "<=" : ">="), RatApproxReal(boundval));

      /* bounds in the problem should be created only once, but row singletons get handled as bounds, as well */
      image = SCIPhashtableRetrieve(certificate->varboundtable, (void*)certificate->workbound);

      if( image != NULL )
      {
         SCIPdebugMessage("Duplicate bound in certificate hashtable. Ignoring \n");
         assert(RatIsEqual(boundval, ((SCIP_CERTIFICATEBOUND*) image)->boundval));

         //return SCIP_OKAY;
      }
      /* add the new bound */
      {
         SCIP_CERTIFICATEBOUND* insertbound;

         SCIP_ALLOC( BMSduplicateBlockMemory(certificate->blkmem, &insertbound, certificate->workbound) );
         SCIP_CALL( RatCopy(certificate->blkmem, &insertbound->boundval, boundval) );
         /* ensure size and insert boundval in array to be able to free it at the end */
         if( certificate->nboundvals >= certificate->boundvalsize )
         {
            BMSreallocBlockMemoryArray(certificate->blkmem, &certificate->boundvals, certificate->boundvalsize, certificate->boundvalsize + 100);
            certificate->boundvalsize += 100;
         }
         certificate->boundvals[certificate->nboundvals] = insertbound;
         certificate->nboundvals++;

         SCIP_CALL( SCIPhashtableInsert(certificate->varboundtable, (void*)insertbound) );
         certificate->indexcounter++;
         certificate->conscounter++;
      }
   }

   if( boundname == NULL )
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "B%d %c ", certificate->indexcounter - 1, (isupper ? 'L' : 'G'));
   else
     SCIPcertificatePrintProblemMessage(certificate, isorigfile, "%s %c ", boundname, (isupper ? 'L' : 'G'));

   SCIPcertificatePrintProblemRational(certificate, isorigfile, boundval, 10);
   SCIPcertificatePrintProblemMessage(certificate, isorigfile, " 1 %d 1\n", varindex);

   return SCIP_OKAY;
}

/** checks whether variable bound assumption is present; prints it if not; returns index */
SCIP_Longint SCIPcertificatePrintBoundAssumption(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           assumptionname,     /**< name of the bound constraint */
   int                   varindex,           /**< index of the variable */
   SCIP_Rational*        boundval,           /**< value of the bound */
   SCIP_Bool             isupper             /**< is it the upper bound? */
   )
{
   void* image;

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return 0;

   /* install bound information in working struct */
   certificate->workbound->fileindex = certificate->indexcounter;
   certificate->workbound->varindex = varindex;
   RatSet(certificate->workbound->boundval, boundval);
   certificate->workbound->isupper = isupper;

   image = SCIPhashtableRetrieve(certificate->varboundtable, (void*)certificate->workbound);
   if( image != NULL )
   {
      SCIP_CERTIFICATEBOUND* foundbound;

      foundbound = (SCIP_CERTIFICATEBOUND*)image;

      SCIPdebugMessage("Found bound assumption at line %" SCIP_LONGINT_FORMAT ": <variable %d> %s approx. %g\n",
         foundbound->fileindex, varindex, (isupper ? "<=" : ">="), RatApproxReal(boundval));

      assert(foundbound->fileindex >= 0);
      assert(foundbound->varindex == varindex);
      assert(RatIsEqual(foundbound->boundval, boundval));
      assert(foundbound->isupper == isupper);
      return foundbound->fileindex;
   }
   else
   {
      SCIP_CERTIFICATEBOUND* insertbound;

      SCIPdebugMessage("Print bound assumption at line %" SCIP_LONGINT_FORMAT ": <variable %d> %s approx. %g\n",
         certificate->workbound->fileindex, varindex, (isupper ? "<=" : ">="), RatApproxReal(boundval));

      SCIP_ALLOC( BMSduplicateBlockMemory(certificate->blkmem, &insertbound, certificate->workbound) );
      SCIP_CALL( RatCopy(certificate->blkmem, &insertbound->boundval, boundval) );

      /* ensure size and insert boundval in array to be able to free it at the end */
      if( certificate->nboundvals >= certificate->boundvalsize )
      {
         BMSreallocBlockMemoryArray(certificate->blkmem, &certificate->boundvals, certificate->boundvalsize, certificate->boundvalsize + 100);
         certificate->boundvalsize += 100;
      }
      certificate->boundvals[certificate->nboundvals] = insertbound;
      certificate->nboundvals++;

      SCIP_CALL( SCIPhashtableInsert(certificate->varboundtable, (void*)insertbound) );

      certificate->indexcounter++;

      if( assumptionname == NULL )
         SCIPcertificatePrintProofMessage(certificate, "A%d %c ", insertbound->fileindex, (isupper ? 'L' : 'G'));
      else
        SCIPcertificatePrintProofMessage(certificate,  "%s %c ", assumptionname, (isupper ? 'L' : 'G'));

      SCIPcertificatePrintProofRational(certificate, boundval, 10);
      SCIPcertificatePrintProofMessage(certificate,  " 1 %d 1 { asm } -1\n", varindex);

      return insertbound->fileindex;
   }
}

/** installs updated node data in parent node */
SCIP_RETCODE SCIPcertificateUpdateParentData(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_NODE*            node,               /**< node data structure */
   SCIP_Longint          fileindex,          /**< index of new bound */
   SCIP_Rational*        newbound            /**< pointer to value of new bound, NULL if infeasible */
   )
{
   SCIP_CERTNODEDATA* nodedataparent;
   SCIP_CERTNODEDATA* nodedata;

   assert(node != NULL);
   assert(fileindex >= 0);

   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   /* Retrieve node data */
   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);

   if( RatIsLT(newbound, nodedata->derbound_inherit) )
      return SCIP_OKAY;

   /* If the node is the root node, then when only update the index and bound */
   if( SCIPnodeGetParent(node) == NULL )
   {
      certificate->derindex_root = fileindex;
      if( newbound == NULL )
         certificate->rootinfeas = TRUE;
      else
         RatSet(certificate->rootbound, newbound);
      return SCIP_OKAY;
   }

   /* Retrieve parent node data */
   assert(SCIPhashmapExists(certificate->nodedatahash, SCIPnodeGetParent(node)));
   nodedataparent = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, SCIPnodeGetParent(node));

   /* First ensure whether the node is left/right child of node parent, then left/rightfilled tells us if,
    * a bound has already been derived for this node, if yes then depending on the value of the bound we update it.
    * Else it is the first time this node is processed.Therefore, we also have to set left/rightfilled to TRUE
    */
   if( certificateIsLeftNode(certificate, node) )
   {
      nodedataparent->leftfilled = TRUE;
      if( RatIsGT(newbound, nodedataparent->derbound_left) )
      {
         nodedataparent->derindex_left = fileindex;
         RatSet(nodedataparent->derbound_left, newbound);
      }
      if( RatIsNegInfinity(newbound) )
      {
         nodedataparent->derindex_left = fileindex;
         nodedataparent->leftinfeas = TRUE;
      }
   }
   else
   {
      nodedataparent->rightfilled = TRUE;
      if( RatIsGT(newbound, nodedataparent->derbound_right) )
      {
         nodedataparent->derindex_right = fileindex;
         RatSet(nodedataparent->derbound_right, newbound);
      }
      if( RatIsNegInfinity(newbound) )
      {
         nodedataparent->rightinfeas = TRUE;
         nodedataparent->derindex_right = fileindex;
      }
   }

   return SCIP_OKAY;
}

/** Print a dual bound from an exact lp solution */
SCIP_RETCODE SCIPcertificatePrintDualboundExactLP(
   SCIP_CERTIFICATE*     certificate,        /**< scip certificate struct */
   SCIP_LPEXACT*         lpexact,            /**< the exact lp */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_NODE*            node,               /**< the current node */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Bool             usefarkas           /**< should an infeasibility proof be printed? */
   )
{
   SCIP_Rational** vals;
   SCIP_Rational* val;
   SCIP_Rational* lowerbound;
   SCIP_Rational* farkasrhs;
   SCIP_Rational* tmp;
   SCIP_Longint* ind;
   int len;
   int i;
   unsigned long key;
   void* image;

   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   SCIP_CALL( RatCreateBuffer(set->buffer, &tmp) );
   SCIPlpExactGetObjval(lpexact, set, prob, tmp);

   /* only print line if bound improved */
   if( !usefarkas && RatIsLT(tmp, SCIPnodeGetLowerboundExact(node)) )
   {
      RatFreeBuffer(set->buffer, &tmp);
      return SCIP_OKAY;
   }

   /* if at root node set, objintegral flag */
   if( SCIPnodeGetParent(node) == NULL )
      certificate->objintegral = SCIPprobIsObjIntegral(prob);

   assert(lpexact!= NULL);
   assert(certificate->transfile != NULL);

   vals = certificate->vals;
   /* if needed extend vals array */
   if( lpexact->ncols + lpexact->nrows > certificate->valssize )
   {
      SCIP_ALLOC( BMSreallocBlockMemoryArray(certificate->blkmem, &certificate->vals,
         certificate->valssize, lpexact->ncols + lpexact->nrows) );
      for( i = certificate->valssize; i < lpexact->ncols + lpexact->nrows; i++ )
      {
         SCIP_CALL( RatCreateBlock(certificate->blkmem, &(certificate->vals[i])) );
      }
      certificate->valssize =  lpexact->ncols + lpexact->nrows;
   }

   SCIP_CALL( RatCreateBuffer(set->buffer, &farkasrhs) );

   SCIPsetAllocBufferArray(set, &ind, lpexact->nrows + lpexact->ncols);

   len = 0;
   for( i = 0; i < lpexact->ncols; ++i )
   {
      SCIP_COLEXACT* col = lpexact->cols[i];
      if( usefarkas )
         val = col->farkascoef;
      else
         val = col->redcost;

      assert(!RatIsAbsInfinity(val));

      if( !RatIsZero(val) )
      {
         if( usefarkas )
            RatNegate(vals[len], val);
         else
            RatSet(vals[len], val);

         ind[len] = getVarBoundFileIndex(certificate, col, vals[len]);

         RatDebugMessage("Column %d for var %s has index %l and farkas coef %q \n", col->index, col->var->name, ind[len], val);

         /* update farkasrhs */
         if( usefarkas )
         {
            val = certificate->workbound->boundval;
            RatAddProd(farkasrhs, vals[len], val);
         }
         len++;
      }
   }

   for( i = 0; i < lpexact->nrows; ++i )
   {
      SCIP_ROWEXACT* row;
      row = lpexact->rows[i];
      val = usefarkas ? row->dualfarkas : row->dualsol;
      assert(!RatIsAbsInfinity(val));

      if( !RatIsZero(val) )
      {
         RatSet(vals[len], val);
         key = (size_t)SCIPhashmapGetImage(certificate->rowdatahash, (void*) row);
         if( key >= INT_MAX - 1 )
         {
            SCIPerrorMessage("row not in rowdata-hash \n");
            SCIPABORT();
            return SCIP_ERROR;
         }

         if( key == 0 && SCIProwExactGetNNonz(row) == 1 )
         {
            key = getVarBoundFileIndex(certificate, SCIProwExactGetCols(row)[0], val);
         }

         assert(key != 0);

         ind[len] = key;
         /* if we have a ranged row, and the dual corresponds to the upper bound,
          * the index for the rhs-constraint is one larger in the certificate */
         if( !RatIsEqual(row->lhs, row->rhs) && !RatIsAbsInfinity(row->lhs) && RatIsNegative(val) )
             ind[len] += 1;

         RatDebugMessage("Row (index %d, %s has index %l and farkas coef %q ", row->index, row->fprow->name, ind[len], val);
         SCIPdebug(SCIProwExactPrint(row, set->scip->messagehdlr, NULL) );

         /* update farkasrhs */
         if( usefarkas )
         {
            val = RatIsPositive(vals[len]) ? row->lhs : row->rhs;
            RatAddProd(farkasrhs, vals[len], val);
            RatDiffProd(farkasrhs, vals[len], row->constant);
         }

         len++;
      }
   }

   SCIP_CALL( RatCreateBuffer(set->buffer, &lowerbound) );
   if( usefarkas )
   {
      RatSetString(lowerbound, "inf");
      /* Scale proof to have RHS = 1 */
      // for( i = 0; i < len; ++i)
      //    RatDiv(vals[i], vals[i], farkasrhs);
   }
   else
   {
      /* vipr does not accept infinity, so in case of objlimit, get the objval from the lpi */
      if( !RatIsInfinity(lpexact->lpobjval) )
         RatSet(lowerbound, lpexact->lpobjval);
      else
      {
         SCIPlpiExactGetObjval(lpexact->lpiexact, lowerbound);
      }
   }

   SCIPcertificatePrintDualbound(certificate, NULL, lowerbound, len, ind, vals);
   SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1, lowerbound);
   SCIPcertificateSetInheritanceData(certificate, node, certificate->indexcounter - 1, lowerbound);

   RatFreeBuffer(set->buffer, &lowerbound);
   SCIPsetFreeBufferArray(set, &ind);
   RatFreeBuffer(set->buffer, &farkasrhs);
   RatFreeBuffer(set->buffer, &tmp);

   return SCIP_OKAY;
}

/** Print a dual bound from the pseudo solution */
SCIP_RETCODE  SCIPcertificatePrintDualboundPseudo(
   SCIP_CERTIFICATE*     certificate,        /**< scip certificate struct */
   SCIP_LPEXACT*         lpexact,            /**< the exact lp */
   SCIP_NODE*            node,               /**< current node */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             psval               /**< the pseudo obj value (or inf to use exact lp value) */
   )
{
   SCIP_VAR** vars;
   SCIP_Rational* pseudoobjval;
   SCIP_Rational** bounds;
   SCIP_Longint* dualind;
   int nvars;
   int duallen;
   int i;
   int nnonzeros;
   void* image;

   assert(certificate != NULL);
   assert(prob != NULL);
   assert(set != NULL);

   /* only print if not -infinity and certificate is active */
   if( !set->exact_enabled || !SCIPcertificateIsActive(set, certificate) || SCIPsetIsInfinity(set, -psval) )
      return SCIP_OKAY;

   if( psval < SCIPnodeGetLowerbound(node) )
      return SCIP_OKAY;

   vars = SCIPprobGetVars(prob);
   nvars = SCIPprobGetNVars(prob);
   SCIP_CALL( RatCreateBuffer(set->buffer, &pseudoobjval) );

   /* infinity means we use the exact lp value */
   if( SCIPsetIsInfinity(set, psval) )
      SCIPlpExactGetPseudoObjval(lpexact, set, prob, pseudoobjval);
   else
      RatSetReal(pseudoobjval, psval);
   duallen = SCIPprobGetNObjVars(prob, set);
   SCIP_CALL( RatCreateBufferArray(set->buffer, &bounds, duallen) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualind, duallen) );
   /* computes pseudo objective value (with bound change if necessary) */
   /*lint --e{838}*/

   nnonzeros = 0;
   for( i = 0; i < nvars; i++ )
   {
      SCIP_Rational* obj = SCIPvarGetObjExact(vars[i]);
      if( !RatIsZero(obj) )
      {
         RatSet(bounds[nnonzeros], obj);

         assert(!RatIsAbsInfinity(bounds[nnonzeros]));

         /* retrieve the line in the certificate of the bound */
         RatSet(certificate->workbound->boundval, SCIPvarGetBestBoundLocalExact(vars[i]));
         certificate->workbound->varindex = SCIPvarGetCertificateIndex(vars[i]);
         certificate->workbound->isupper = RatIsNegative(obj);

         image = SCIPhashtableRetrieve(certificate->varboundtable, (void*)certificate->workbound);

         if( image != NULL )
            dualind[nnonzeros] = ((SCIP_CERTIFICATEBOUND*)image)->fileindex;
         else
         {
            SCIPerrorMessage("Bound should be present in certificate \n");
            SCIPABORT();
            return SCIP_ERROR;
         }
         nnonzeros++;
      }
   }
   assert(nnonzeros == duallen);
   /* print pseudo solution into certificate file */
   SCIPcertificatePrintDualbound(certificate, NULL, pseudoobjval, duallen,
      dualind, bounds);

   SCIP_CALL( SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1,
      pseudoobjval) );

   SCIPcertificateSetInheritanceData(certificate, node, certificate->indexcounter - 1, pseudoobjval);

   SCIPsetFreeBufferArray(set, &dualind);
   RatFreeBufferArray(set->buffer, &bounds, nnonzeros);
   RatFreeBuffer(set->buffer, &pseudoobjval);

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPcertificatePrintInheritedBound(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_NODE*            node                /**< node data */
   )
{

   SCIP_CERTNODEDATA* nodedata;
   SCIP_Longint unsplit_index;
   SCIP_Rational* lowerbound;
   SCIP_Bool infeas;

   assert(node != NULL);

   /* check whether certificate output should be created */
   if( certificate->transfile == NULL || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* get the current node data */
   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*) SCIPhashmapGetImage(certificate->nodedatahash, node);
   infeas = FALSE;
   if( nodedata->inheritedbound && nodedata->assumptionindex_self != - 1 )
   {
      SCIP_Longint ind[1];
      ind[0] = nodedata->derindex_inherit;
      SCIP_Rational* lowerbound;
      SCIP_Rational* val;

      RatCreateBuffer(set->buffer, &lowerbound);
      RatCreateBuffer(set->buffer, &val);

      RatSet(lowerbound, nodedata->derbound_inherit);
      RatSetInt(val, 1, 1);

      SCIPcertificatePrintDualbound(certificate, NULL, lowerbound, 1, ind, &val);
      SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1, lowerbound);

      RatFreeBuffer(set->buffer, &lowerbound);
      RatFreeBuffer(set->buffer, &val);
   }

   return SCIP_OKAY;
}

/** prints dual bound to proof section */
SCIP_Longint SCIPcertificatePrintDualbound(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   const char*           linename,           /**< name of the unsplitting line */
   SCIP_Rational*        lowerbound,         /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   int                   len,                /**< number of dual multipiers */
   SCIP_Longint*         ind,                /**< index array */
   SCIP_Rational**       val                 /**< array of dual multipliers */
   )
{
   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return 0;

   certificate->indexcounter++;

   if( linename == NULL )
   {
      SCIPcertificatePrintProofMessage(certificate, "L%d ", certificate->indexcounter - 1);
   }
   else
   {
      SCIPcertificatePrintProofMessage(certificate, "%s ", linename);
   }

   if( RatIsInfinity(lowerbound) )
   {
      SCIPcertificatePrintProofMessage(certificate, "G 1 0");
   }
   else
   {
      SCIPcertificatePrintProofMessage(certificate, "G ");
      SCIPcertificatePrintProofRational(certificate, lowerbound, 10);
      SCIPcertificatePrintProofMessage(certificate, " ");
      SCIPcertificatePrintProofMessage(certificate, "OBJ");
   }

   if( val == NULL )
   {
      SCIPcertificatePrintProofMessage(certificate, " { lin ... } -1\n");
   }
   else
   {
      SCIPcertificatePrintProofMessage(certificate, " { lin %d", len);
      int i;
      for( i = 0; i < len; i++ )
      {
         /** @todo exip: perform line breaking before exceeding maximum line length */
         assert(!RatIsAbsInfinity(val[i]));
         SCIPcertificatePrintProofMessage(certificate, " %d ", ind[i]);
         SCIPcertificatePrintProofRational(certificate, val[i], 10);
      }
      SCIPcertificatePrintProofMessage(certificate, " } -1\n");
   }

   /* print rounding derivation */
   if( !RatIsNegInfinity(lowerbound) && certificate->objintegral && !RatIsIntegral(lowerbound) )
   {
      long int ceilint;

      certificate->indexcounter++;

      SCIPcertificatePrintProofMessage(certificate, "R%d G ", certificate->indexcounter - 1);
      updateFilesize(certificate, 4.0 + ceil(log10(certificate->indexcounter - 1 + 1)));
      RatRound(lowerbound, lowerbound, SCIP_ROUND_UPWARDS);

      SCIPcertificatePrintProofRational(certificate, lowerbound, 10);

      SCIPcertificatePrintProofMessage(certificate, " ");
      SCIPcertificatePrintProofMessage(certificate, "OBJ");
      SCIPcertificatePrintProofMessage(certificate, " { rnd 1 %d 1 } -1\n", certificate->indexcounter - 2);
   }

   return (certificate->indexcounter - 1);
}

/** update the parent certificate node data when branching, print branching into certificate if not already present */
SCIP_RETCODE SCIPcertificatePrintBranching(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP informations */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_NODE*            node,               /**< node data */
   SCIP_VAR*             branchvar,          /**< the variable that gets branched on */
   SCIP_BOUNDTYPE        boundtype,          /**< the bounding type */
   SCIP_Real             newbound            /**< the new bound */
   )
{
   SCIP_CERTNODEDATA* nodedataparent;
   SCIP_CERTNODEDATA* nodedata;
   SCIP_Rational* branchbound;

   assert(node != NULL);
   assert(stat != NULL);

   /* certificate is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE || lp->diving )
      return SCIP_OKAY;

   /* check whether output should be created */
   if ( certificate->transfile == NULL )
      return SCIP_OKAY;

   SCIP_CALL( RatCreateBuffer(set->buffer, &branchbound) );
   RatSetReal(branchbound, newbound);

   assert(RatIsIntegral(branchbound));

   nodedata = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, node);

   if( branchvar != NULL )
   {
      nodedata->assumptionindex_self = printBoundAssumption(set, certificate, SCIPvarGetCertificateIndex(branchvar),
         branchbound, boundtype);
   }

   if( SCIPnodeGetParent(node) != NULL && certificate->transfile != NULL )
   {
      nodedataparent = (SCIP_CERTNODEDATA*)SCIPhashmapGetImage(certificate->nodedatahash, SCIPnodeGetParent(node));

      if( branchvar != NULL )
      {
         if( boundtype == SCIP_BOUNDTYPE_LOWER )
            nodedataparent->assumptionindex_right = nodedata->assumptionindex_self;
         if( boundtype == SCIP_BOUNDTYPE_UPPER )
            nodedataparent->assumptionindex_left = nodedata->assumptionindex_self;
      }
   }

   RatFreeBuffer(set->buffer, &branchbound);

   return SCIP_OKAY;
}

/** create a new node data structure for the current node */
SCIP_RETCODE SCIPcertificateNewNodeData(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< SCIP certificate */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_NODE*            node                /**< new node, that was created */
   )
{
   SCIP_CERTNODEDATA* nodedata;

   assert(stat != NULL );
   assert(node != NULL );

   /* certificate is disabled on probing nodes */
   if( SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* check whether output should be created */
   if( certificate->transfile == NULL )
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocBlockMemory(certificate->blkmem, &nodedata) );

   nodedata->derindex_left = -1;
   nodedata->derindex_right = -1;
   nodedata->assumptionindex_left = -1;
   nodedata->assumptionindex_right = -1;
   SCIP_CALL( RatCreateString(certificate->blkmem, &nodedata->derbound_left, "-inf") );
   SCIP_CALL( RatCreateString(certificate->blkmem, &nodedata->derbound_right, "-inf") );
   SCIP_CALL( RatCreateString(certificate->blkmem, &nodedata->derbound_inherit, "-inf") );
   nodedata->assumptionindex_self = -1;
   nodedata->leftinfeas = FALSE;
   nodedata->leftfilled = FALSE;
   nodedata->rightinfeas = FALSE;
   nodedata->rightfilled = FALSE;
   nodedata->inheritedbound = TRUE;
   nodedata->derindex_inherit = -1;
   if( SCIPnodeGetParent(node) != NULL )
   {
      SCIP_NODE* parent = SCIPnodeGetParent(node);
      SCIP_CERTNODEDATA* parentdata;
      assert(SCIPhashmapExists(certificate->nodedatahash, parent));
      parentdata = SCIPhashmapGetImage(certificate->nodedatahash, (void*) parent);
      assert(parentdata != NULL);

      nodedata->derindex_inherit = parentdata->derindex_inherit;
      RatSet(nodedata->derbound_inherit, parentdata->derbound_inherit);
   }

   /* link the node to its nodedata in the corresponding hashmap */
   SCIP_CALL( SCIPhashmapSetImage(certificate->nodedatahash, node, (void*)nodedata) );

   return SCIP_OKAY;
}

/** prints unsplitting information to proof section */
int SCIPcertificatePrintUnsplitting(
   SCIP_SET*             set,                /**< general SCIP settings */
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_NODE*            node                /**< node data */
   )
{
   SCIP_CERTNODEDATA* nodedata;
   SCIP_Longint unsplit_index;
   SCIP_Rational* lowerbound;
   SCIP_Bool infeas;

   assert(node != NULL);

   /* check whether certificate output should be created */
   if( certificate->transfile == NULL || SCIPnodeGetType(node) == SCIP_NODETYPE_PROBINGNODE )
      return SCIP_OKAY;

   /* get the current node data */
   assert(SCIPhashmapExists(certificate->nodedatahash, node));
   nodedata = (SCIP_CERTNODEDATA*) SCIPhashmapGetImage(certificate->nodedatahash, node);
   SCIP_CALL( RatCreateBuffer(set->buffer, &lowerbound) );
   infeas = FALSE;

   assert(nodedata != NULL);

   if( nodedata->leftfilled && nodedata->rightfilled )
   {
      if( nodedata->leftinfeas && nodedata->rightinfeas )
      {
         infeas = TRUE;
         RatSetString(lowerbound, "inf");
      }
      else if( nodedata->leftinfeas )
         RatSet(lowerbound, nodedata->derbound_right);
      else if( nodedata->rightinfeas )
         RatSet(lowerbound, nodedata->derbound_left);
      else
         RatMIN(lowerbound, nodedata->derbound_left, nodedata->derbound_right);

      if( RatIsInfinity(nodedata->derbound_left) && RatIsInfinity(nodedata->derbound_right) )
         infeas = TRUE;

      certificate->indexcounter++;

      SCIPcertificatePrintProofMessage(certificate, "U%d ", certificate->indexcounter - 1);

      if( infeas )
      {
         SCIPcertificatePrintProofMessage(certificate, "G 1 0");
      }
      else
      {
         SCIPcertificatePrintProofMessage(certificate, "G ");
         SCIPcertificatePrintProofRational(certificate, lowerbound, 10);
         SCIPcertificatePrintProofMessage(certificate, " ");
         SCIPcertificatePrintProofMessage(certificate, "OBJ");
      }

      assert(nodedata->derindex_right < certificate->indexcounter - 1);
      assert(nodedata->derindex_left < certificate->indexcounter - 1);

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
         ind[0] = nodedata->derindex_inherit;
         SCIP_Rational* lowerbound;
         SCIP_Rational* val;

         RatCreateBuffer(set->buffer, &lowerbound);
         RatCreateBuffer(set->buffer, &val);

         RatSet(lowerbound, nodedata->derbound_inherit);
         RatSetInt(val, 1, 1);

         SCIPcertificatePrintDualbound(certificate, NULL, lowerbound, 1, ind, &val);
         SCIPcertificateUpdateParentData(certificate, node, certificate->indexcounter - 1, lowerbound);

         RatFreeBuffer(set->buffer, &lowerbound);
         RatFreeBuffer(set->buffer, &val);
      }
   }

   certificateFreeNodeData(certificate, node);

   RatFreeBuffer(set->buffer, &lowerbound);
   return (certificate->indexcounter - 1);
}

/** prints RTP section with lowerbound and upperbound range */
void SCIPcertificatePrintRtpRange(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   SCIP_Rational*        lowerbound,         /**< pointer to lower bound on the objective */
   SCIP_Rational*        upperbound          /**< pointer to upper bound on the objective */
   )
{
   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "RTP range ");
   if( RatIsNegInfinity(lowerbound) )
   {
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "-inf");
   }
   else
   {
      SCIPcertificatePrintProblemRational(certificate, isorigfile, lowerbound, 10);
   }

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, " ");
   if( RatIsInfinity(upperbound) )
   {
      SCIPcertificatePrintProblemMessage(certificate, isorigfile, "inf");
   }
   else
   {
      SCIPcertificatePrintProblemRational(certificate, isorigfile, upperbound, 10);
   }
   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "\n");
 }

/** prints RTP section for infeasibility */
void SCIPcertificatePrintRtpInfeas(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_Bool             isorigfile          /**< should the original solution be printed or in transformed space */
   )
{
   /* check if certificate output should be created */
   if( certificate->transfile == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, isorigfile, "RTP infeas\n");
 }
