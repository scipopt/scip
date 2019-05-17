/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
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
#include "scip/scip.h"
#include "scip/set.h"
#include "scip/lpex.h"
#include "scip/pub_misc.h"
#include "scip/prob.h"
#include "scip/certificate.h"
#include "scip/struct_certificate.h"
#include "scip/solex.h"
#include "scip/struct_scip.h"
#include "scip/pub_varex.h"

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

   if( !RisEqual(((SCIP_CERTIFICATEBOUND*)key1)->boundval, ((SCIP_CERTIFICATEBOUND*)key2)->boundval) )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValVarbound)
{
   SCIP_CERTIFICATEBOUND* bound;

   bound = (SCIP_CERTIFICATEBOUND*)key;

   return ((((unsigned int)(10*RgetRealApprox(bound->boundval))) << 22) + (bound->varindex << 2) + (unsigned int)bound->isupper);
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
   (*certificate)->file = NULL;
   (*certificate)->derivationfile = NULL;
   (*certificate)->derivationfilename = NULL;
   (*certificate)->objstring = NULL;
   (*certificate)->filesize = 0.0;
   (*certificate)->rowdatahash = NULL;
   (*certificate)->boundvals = NULL;
   (*certificate)->boundvalsize = 0;

   return SCIP_OKAY;
}

/** frees certificate data structure */
void SCIPcertificateFree(
   SCIP_CERTIFICATE**    certificate         /**< pointer to store the certificate information */
   )
{
   assert(certificate != NULL);
   assert(*certificate != NULL);
   assert((*certificate)->file == NULL);
   assert((*certificate)->derivationfile == NULL);
   assert((*certificate)->objstring == NULL);

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
   int j;
   char* name = NULL;
   char* compression = NULL;
   SCIP_VAR** vars;
   SCIP_Rational* lb;
   SCIP_Rational* ub;

   assert(certificate != NULL);
   assert(set != NULL);
   assert(set->certificate_filename != NULL);
   assert(certificate->derivationfile == NULL);

   if( set->certificate_filename[0] == '-' && set->certificate_filename[1] == '\0' )
      return SCIP_OKAY;

   filenamelen = strlen(set->certificate_filename);

   SCIP_CALL( SCIPsetAllocBufferArray(set, &name, filenamelen + 1) );
   BMScopyMemoryArray(name, set->certificate_filename, filenamelen);

   SCIPsplitFilename(name, NULL, NULL, NULL, &compression);

   SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      "storing certificate information in file <%s>\n", set->certificate_filename);

   if( NULL != compression && 0 == strncmp(compression, "gz", 2) )
      certificate->file = SCIPfopen(set->certificate_filename, "wb");
   else
      certificate->file = SCIPfopen(set->certificate_filename, "wT");

   bufferlen = strlen(name);
   SCIP_ALLOC( BMSallocMemoryArray(&certificate->derivationfilename, filenamelen+5) );
   BMScopyMemoryArray(certificate->derivationfilename, name, bufferlen);
   certificate->derivationfilename[bufferlen] = '_';
   certificate->derivationfilename[bufferlen+1] = 'd';
   certificate->derivationfilename[bufferlen+2] = 'e';
   certificate->derivationfilename[bufferlen+3] = 'r';
   if( NULL != compression && 0 == strncmp(compression, "gz", 2) )
   {
      certificate->derivationfilename[bufferlen+4] = '.';
      certificate->derivationfilename[bufferlen+5] = 'g';
      certificate->derivationfilename[bufferlen+6] = 'z';
      certificate->derivationfilename[bufferlen+7] = '\0';

      certificate->derivationfile = SCIPfopen(certificate->derivationfilename, "wb");
   }
   else
   {
      certificate->derivationfilename[bufferlen+4] = '\0';
      certificate->derivationfile = SCIPfopen(certificate->derivationfilename, "wT");
   }

   if( certificate->file == NULL || certificate->derivationfile == NULL )
   {
      SCIPerrorMessage("error creating file <%s> and derivation file\n", set->certificate_filename);
      SCIPprintSysError(set->certificate_filename);
      return SCIP_FILECREATEERROR;
   }

   assert(certificate->rowdatahash == NULL);
   SCIP_CALL( SCIPhashmapCreate(&certificate->rowdatahash, blkmem, SCIP_HASHSIZE_CERTIFICATE) );

   certificate->blkmem = blkmem;
   SCIPsetFreeBufferArray(set, &name);

   SCIP_CALL( SCIPgetOrigVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   nboundconss = 0;
   for ( j = 0 ; j < nvars ; j++ )
   {
      lb = SCIPvarGetLbGlobalExact(vars[j]);
      ub = SCIPvarGetUbGlobalExact(vars[j]);
      if( !RisAbsInfinity(lb) )
         nboundconss++;
      if( !RisAbsInfinity(ub) )
         nboundconss++;
   }

   /* print the Version Header into certificate */
   SCIPcertificatePrintVersionHeader(certificate);

   /* print the Variable Header into certificate */
   SCIPcertificatePrintVarHeader(certificate, nvars);
   for( j = 0; j < nvars; j++ )
   {
      const char* varname;

      varname = SCIPvarGetName(vars[j]);
      if( strstr(varname, " ") != NULL || strstr(varname, "\t") != NULL || strstr(varname, "\n") != NULL
         || strstr(varname, "\v") != NULL || strstr(varname, "\f") != NULL || strstr(varname, "\r") != NULL )
      {
         SCIPerrorMessage("Variable name <%s> cannot be printed to certificate file because it contains whitespace.\n",
            varname);
         return SCIP_ERROR;
      }

      SCIPcertificatePrintProblemMessage(certificate, "%s\n", varname);
   }

   /* print the Integer Variable Header into certificate */
   SCIPcertificatePrintIntHeader(certificate, nintvars + nbinvars);
   for( j = 0; j < nvars; j++ )
   {
      if( SCIPvarGetType(vars[j]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(vars[j]) == SCIP_VARTYPE_INTEGER )
      {
         SCIPcertificatePrintProblemMessage(certificate, "%d \n", SCIPvarGetIndex(vars[j]));
      }
   }

   {
      SCIP_Rational** objcoefs = RcreateArrayTemp(SCIPbuffer(scip), nvars);
      for( j = 0; j < nvars; j++)
         Rset(objcoefs[j], SCIPvarGetObjExact(vars[j]));

      /* print the objective function into certificate header */
      SCIPcertificateSetAndPrintObjective(certificate, blkmem, objcoefs, nvars);

      RdeleteArrayTemp(SCIPbuffer(scip), &objcoefs, nvars);
   }

   SCIPcertificatePrintConsHeader(certificate, SCIPgetNConss(scip), nboundconss);

   for( j = 0; j < nvars; j++ )
   {
      if( !RisAbsInfinity(SCIPvarGetLbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, NULL, SCIPvarGetIndex(vars[j]), SCIPvarGetLbGlobalExact(vars[j]), FALSE) );
      }
      if( RisAbsInfinity(SCIPvarGetUbGlobalExact(vars[j])) )
      {
         SCIP_CALL( SCIPcertificatePrintBoundCons(certificate, NULL, SCIPvarGetIndex(vars[j]), SCIPvarGetUbGlobalExact(vars[j]), TRUE) );
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
      SCIPfwrite(buffer, sizeof(char), size, certificate->file);

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

   if( certificate->file != NULL )
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
      SCIPfclose(certificate->file);
      certificate->file = NULL;

      BMSfreeMemoryArray(&certificate->derivationfilename);
      if( certificate->varboundtable != NULL )
      {
         for( i = 0; i < SCIPhashtableGetNElements(certificate->varboundtable); i++ )
         {
            Rdelete(certificate->blkmem, &certificate->boundvals[i]->boundval);
            BMSfreeBlockMemory(certificate->blkmem, &certificate->boundvals[i]);
         }
         /**@todo fix memory leak: mpq_clear and free all elements */
         SCIPhashtableRemoveAll(certificate->varboundtable);
         SCIPhashtableFree(&certificate->varboundtable);
         BMSfreeBlockMemoryArray(certificate->blkmem, &certificate->boundvals, certificate->boundvalsize);
      }
      if( certificate->workbound != NULL )
      {
         Rdelete(certificate->blkmem, &certificate->workbound->boundval);
         BMSfreeBlockMemory(certificate->blkmem, &certificate->workbound);
      }
      BMSfreeMemoryArrayNull(&certificate->objstring);

      if( certificate->rowdatahash)
         SCIPhashmapFree(&certificate->rowdatahash);
   }
}

/** returns whether the certificate output is activated? */
SCIP_Bool SCIPcertificateIsActive(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   return (certificate != NULL && certificate->file != NULL);
}

/** returns current certificate file size in MB */
SCIP_Real SCIPcertificateGetFilesize(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   if( certificate == NULL || certificate->file == NULL)
      return 0.0;
   else
      return certificate->filesize;
}

/** sets the objective function used when printing dual bounds */
SCIP_RETCODE SCIPcertificateSetAndPrintObjective(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Rational**       coefs,              /**< objective function coefficients */
   int                   nvars               /**< number of variables */
   )
{
   char obj[SCIP_MAXSTRLEN - 2];
   char* objstring;
   int printlen;
   int nnonz;
   int buffpos;
   int i;
   int buflength;
   int allocsize;
   int leftsize;

   assert(coefs != NULL);

   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return SCIP_OKAY;

   /* create a hash table for the variable bounds */
   SCIP_CALL( SCIPhashtableCreate(&certificate->varboundtable, blkmem, nvars,
         hashGetKeyVarbound, hashKeyEqVarbound, hashKeyValVarbound, NULL) );

   /* create working memory for bound struct */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &certificate->workbound) );
   certificate->workbound->boundval = Rcreate(blkmem);

   nnonz = 0;
   buflength = SCIP_MAXSTRLEN - 2;
   for( i = 0; i < nvars; i++ )
   {
      if( !RisZero(coefs[i]) )
         nnonz++;
   }

   allocsize = nnonz * 1024;
   if( certificate->objstring == NULL )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&certificate->objstring, allocsize) );
   }
   else
   {
      SCIP_ALLOC( BMSreallocMemoryArray(&certificate->objstring, allocsize) );
   }

   certificate->objstring[0] = '\0';
   objstring = certificate->objstring;
   leftsize = allocsize;

   buffpos = 0;
   printlen = SCIPsnprintf(obj, SCIP_MAXSTRLEN - 2, "OBJ min\n %d ", nnonz);
   if( printlen >= leftsize )
   {
      SCIPerrorMessage("Objective function string exceeds %d characters: currently not implemented.\n",
         allocsize);
      return SCIP_ERROR;
   }
   else
   {
      buffpos += printlen;
   }

   for( i = 0; i < nvars; i++ )
   {
      if( !RisZero(coefs[i]) )
      {
         printlen = gmp_snprintf(&obj[buffpos], SCIP_MAXSTRLEN - buffpos, "%s%d %Qd", (i > 0 ? " " : ""), i, *RgetGMP(coefs[i]));
         if( printlen >= SCIP_MAXSTRLEN - buffpos )
         {
            obj[buffpos] = '\0';
            SCIPcertificatePrintProblemMessage(certificate, "%s \n", obj);
            buffpos = 0;
            i--;
         }
         else
         {
            objstring += printlen;
            buffpos += printlen;
            leftsize -= printlen;
         }
      }
   }

   SCIPcertificatePrintProblemMessage(certificate, "%s \n", obj); 

   return SCIP_OKAY;
}

/** prints a string to the problem section of the certificate file */
void SCIPcertificatePrintProblemMessage(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[3 * SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return;

   va_start(ap, formatstr);
   vsnprintf(buffer, 3 * SCIP_MAXSTRLEN, formatstr, ap);

   SCIPfprintf(certificate->file, "%s", buffer);
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
   SCIP_Rational*        val,                /**< Rational to print to the problem*/
   int                   base                /**< The base representation*/
   )
{
   char formatstr[SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->file == NULL )
     return;
   RtoString(val, formatstr);
   SCIPcertificatePrintProblemMessage(certificate, "%s", formatstr);
}


/** prints a rational number to the proof section of the certificate file */
void SCIPcertificatePrintProofRational(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   SCIP_Rational*        val,                /**< Rational to print to the problem*/
   int                   base                /**< The base representation*/
   )
{
   char formatstr[SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->derivationfile == NULL )
     return;
   RtoString(val, formatstr);
   SCIPcertificatePrintProofMessage(certificate, "%s", formatstr);
}

/** prints a comment to the problem section of the certificate file */
void SCIPcertificatePrintProblemComment(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;
   char buffer[SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return;
   
   SCIPfprintf(certificate->file, "# ");

   va_start(ap, formatstr);
   vsprintf(buffer, formatstr, ap);

   SCIPfprintf(certificate->file, "%s", formatstr); // todo: is this correct?
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
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   assert(certificate != NULL);

   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, "VER 1.0 \n");
}

/** prints variable section header */
void SCIPcertificatePrintVarHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   int                   nvars               /**< number of variables */
   )
{
   assert(certificate != NULL);
   assert(nvars >= 0);

   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, "VAR %d \n", nvars);
}

/** prints integer section header */
void SCIPcertificatePrintIntHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   int                   nints               /**< number of integer variables */
   )
{
   assert(certificate != NULL);
   assert(nints >= 0);

   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, "INT %d\n", nints);
}

/** prints constraint section header */
void SCIPcertificatePrintConsHeader(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   int                   nconss,             /**< number of all constraints */
   int                   nboundconss         /**< number of bound constraints */
   )
{
   assert(certificate != NULL);
   assert(nconss >= 0);

   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, "CON %d %d\n", nconss, nboundconss);
}

/** prints derivation section header */
void SCIPcertificatePrintDerHeader(
   SCIP_CERTIFICATE*     certificate         /**< certificate information */
   )
{
   int nders;
   assert(certificate != NULL);

   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return;

   nders = certificate->indexcounter - certificate->conscounter;
   SCIPcertificatePrintProblemMessage(certificate, "DER %d\n", nders);
}

/** prints constraint */
void SCIPcertificatePrintCons(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
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
   if( certificate->file == NULL )
      return;

   if( consname == NULL )
      SCIPcertificatePrintProblemMessage(certificate, "C%d %c ", certificate->indexcounter, sense);
   else
      SCIPcertificatePrintProblemMessage( certificate, "%s %c ", consname, sense);

   SCIPcertificatePrintProblemRational(certificate, side, 10);

   SCIPcertificatePrintProblemMessage(certificate, " %d", len);

   for( i = 0; i < len; i++ )
   {
      /* TODO: perform line breaking before exceeding maximum line length */
      SCIPcertificatePrintProblemMessage(certificate, " %d ", ind[i]);
      SCIPcertificatePrintProblemRational(certificate, val[i], 10);
   }
   SCIPcertificatePrintProblemMessage(certificate, "\n");

   certificate->indexcounter++;
   certificate->conscounter++;
}

/** prints a variable bound to the problem section of the certificate file and returns line index */
SCIP_RETCODE SCIPcertificatePrintBoundCons(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const char*           boundname,          /**< name of the bound constraint */
   int                   varindex,           /**< index of the variable */
   SCIP_Rational*        boundval,           /**< value of the bound */
   SCIP_Bool             isupper             /**< is it the upper bound? */
   )
{
   void* image;

   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return SCIP_OKAY;

   /* install bound information in working struct */
   certificate->workbound->fileindex = certificate->indexcounter;
   certificate->workbound->varindex = varindex;
   Rset(certificate->workbound->boundval, boundval);
   certificate->workbound->isupper = isupper;

   SCIPdebugMessage("Printing bound at line %" SCIP_LONGINT_FORMAT ": <variable %d> %s approx. %g\n",
      certificate->indexcounter, varindex, (isupper ? "<=" : ">="), RgetRealApprox(boundval));

   /* bounds in the problem should be created only once */
   image = SCIPhashtableRetrieve(certificate->varboundtable, (void*)certificate->workbound);
   if( image == NULL )
   {
      SCIP_CERTIFICATEBOUND* insertbound;

      SCIP_ALLOC( BMSduplicateBlockMemory(certificate->blkmem, &insertbound, certificate->workbound) );
      insertbound->boundval = Rcopy(certificate->blkmem, boundval);
      /* ensure size and insert boundval in array to be able to free it at the end */
      if( SCIPhashtableGetNElements(certificate->varboundtable) >= certificate->boundvalsize )
      {
         BMSreallocBlockMemoryArray(certificate->blkmem, &certificate->boundvals, certificate->boundvalsize, certificate->boundvalsize + 100);
         certificate->boundvalsize += 100;
      }
      certificate->boundvals[SCIPhashtableGetNElements(certificate->varboundtable)] = insertbound;

      SCIP_CALL( SCIPhashtableInsert(certificate->varboundtable, (void*)insertbound) );
      certificate->indexcounter++;
      certificate->conscounter++;
   }
   else
   {
      SCIPerrorMessage("Duplicate bound in certificate hashtable.\n");
      return SCIP_ERROR;
   }

   if( boundname == NULL )
      SCIPcertificatePrintProblemMessage(certificate, "B%d %c ", certificate->workbound->fileindex, (isupper ? 'L' : 'G'));
   else
     SCIPcertificatePrintProblemMessage(certificate, "%s %c ", boundname, (isupper ? 'L' : 'G'));

   SCIPcertificatePrintProblemRational(certificate, boundval, 10);
   SCIPcertificatePrintProblemMessage(certificate, " 1 %d 1\n", varindex);

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
   if( certificate->file == NULL )
      return 0;

   /* install bound information in working struct */
   certificate->workbound->fileindex = certificate->indexcounter;
   certificate->workbound->varindex = varindex;
   Rset(certificate->workbound->boundval, boundval);
   certificate->workbound->isupper = isupper;

   image = SCIPhashtableRetrieve(certificate->varboundtable, (void*)certificate->workbound);
   if( image != NULL )
   {
      SCIP_CERTIFICATEBOUND* foundbound;

      foundbound = (SCIP_CERTIFICATEBOUND*)image;

      SCIPdebugMessage("Found bound assumption at line %" SCIP_LONGINT_FORMAT ": <variable %d> %s approx. %g\n",
         foundbound->fileindex, varindex, (isupper ? "<=" : ">="), RgetRealApprox(boundval));

      assert(foundbound->fileindex >= 0);
      assert(foundbound->varindex == varindex);
      assert(RisEqual(foundbound->boundval, boundval));
      assert(foundbound->isupper == isupper);
      return foundbound->fileindex;
   }
   else
   {
      SCIP_CERTIFICATEBOUND* insertbound;

      SCIPdebugMessage("Print bound assumption at line %" SCIP_LONGINT_FORMAT ": <variable %d> %s approx. %g\n",
         certificate->workbound->fileindex, varindex, (isupper ? "<=" : ">="), RgetRealApprox(boundval));

      SCIP_ALLOC( BMSduplicateBlockMemory(certificate->blkmem, &insertbound, certificate->workbound) );
      insertbound->boundval = Rcopy(certificate->blkmem, boundval);
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

/** Print a dual bound from an exact lp solution */
SCIP_RETCODE SCIPcertificatePrintDualboundExactLP(
   SCIP_CERTIFICATE*     certificate,        /**< scip certificate struct */
   SCIP_LPEX*            lpex,               /**< the exact lp */
   SCIP_SET*             set,                /**< scip settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   SCIP_Rational** vals;
   int* ind;
   int len;
   int i;
   unsigned long key;
   SCIP_Rational* lowerbound;

   assert(lpex!= NULL);
   assert(certificate->file != NULL);

   vals = RcreateArrayTemp(set->buffer, lpex->nrows);
   SCIPsetAllocBufferArray(set, &ind, lpex->nrows);

   len = 0;
   for( i = 0; i < lpex->nrows; ++i)
   {
      SCIP_ROWEX* row;
      row = lpex->rows[i];

      if( !RisZero(row->dualsol) )
      {
         Rset(vals[i], row->dualsol);
         key = (size_t)SCIPhashmapGetImage(certificate->rowdatahash, (void*) row);
         if( key >= INT_MAX - 1 )
         {
            SCIPerrorMessage("row not in rowdata-hash \n");
            SCIPABORT();
            return SCIP_ERROR;
         }

         ind[len] = key;
         /* if we have a ranged row, and the dual corresponds to the upper bound,
          * the index for the rhs-constraint is one larger in the certificate */
         if( !RisEqual(row->lhs, row->rhs) && !RisAbsInfinity(row->lhs) && RisNegative(row->dualsol) )
             ind[len] += 1;
         
         len++;
      }
   }

   lowerbound = lpex->lpobjval;

   SCIPcertificatePrintDualbound(certificate, prob, NULL, lowerbound, len, ind, vals);

   SCIPsetFreeBufferArray(set, &ind);
   RdeleteArrayTemp(set->buffer, &vals, lpex->nrows);

   return SCIP_OKAY;
}

/** prints dual bound to proof section */
SCIP_Longint SCIPcertificatePrintDualbound(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_PROB*            prob,               /**< problem data */
   const char*           linename,           /**< name of the unsplitting line */
   SCIP_Rational*        lowerbound,         /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   int                   len,                /**< number of dual multipiers */
   int*                  ind,                /**< index array */
   SCIP_Rational**       val                 /**< array of dual multipliers */
   )
{
   /* check if certificate output should be created */
   if( certificate->file == NULL )
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

   if( lowerbound == NULL )
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
         /* TODO: perform line breaking before exceeding maximum line length */
         SCIPcertificatePrintProofMessage(certificate, " %d ", ind[i]);
         SCIPcertificatePrintProofRational(certificate, val[i], 10);
      }
      SCIPcertificatePrintProofMessage(certificate, " } -1\n");
   }

   /* print rounding derivation */
   if( lowerbound != NULL && SCIPprobIsObjIntegral(prob) && !RisIntegral(lowerbound) )
   {
      long int ceilint;

      certificate->indexcounter++;

      SCIPcertificatePrintProofMessage(certificate, "R%d G ", certificate->indexcounter - 1);
      updateFilesize(certificate, 4.0 + ceil(log10(certificate->indexcounter - 1 + 1)));

      if( !RroundInteger(&ceilint, lowerbound, SCIP_ROUND_UPWARDS) )
      {
         SCIPerrorMessage("too large to be represented as long long \n");
         SCIPABORT();
      }
      SCIPcertificatePrintProofMessage(certificate, "%ld", ceilint);

      SCIPcertificatePrintProofMessage(certificate, " ");
      SCIPcertificatePrintProofMessage(certificate, "OBJ");
      SCIPcertificatePrintProofMessage(certificate, " { rnd 1 %d 1 } -1\n", certificate->indexcounter - 2);
   }

   return (certificate->indexcounter - 1);
}

/** prints unsplitting information to proof section */
int SCIPcertificatePrintUnsplitting(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   const char*           linename,           /**< name of the unsplitting line */
   SCIP_Rational*        lowerbound,         /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   int                   derindex_left,      /**< index of the first derivation */
   int                   assumptionindex_left,/**< index of the first unsplitting assumption */
   int                   derindex_right,     /**< index of the second derivation */
   int                   assumptionindex_right/**< index of the second unsplitting assumption */
   )
{
   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return 0;

   certificate->indexcounter++;

   if( linename == NULL )
   {
      SCIPcertificatePrintProofMessage(certificate, "U%d ", certificate->indexcounter - 1);
   }
   else
   {
      SCIPcertificatePrintProofMessage(certificate, "%s ", linename);
   }

   if( lowerbound == NULL )
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

   assert(derindex_right < certificate->indexcounter - 1);
   assert(derindex_left < certificate->indexcounter - 1);

   SCIPcertificatePrintProofMessage(certificate, " { uns %d %d  %d %d  } -1\n", derindex_left, assumptionindex_left, derindex_right, assumptionindex_right);

   return (certificate->indexcounter - 1);
}

/** prints RTP section with lowerbound and upperbound range */
void SCIPcertificatePrintRtpRange(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP_Rational*        lowerbound,         /**< pointer to lower bound on the objective */
   SCIP_Rational*        upperbound          /**< pointer to upper bound on the objective */
   )
{
   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, "RTP range ");
   if( RisNegInfinity(lowerbound) )
   {
      SCIPcertificatePrintProblemMessage(certificate, "-inf");
   }
   else
   {
      SCIPcertificatePrintProblemRational(certificate, lowerbound, 10);
   }

   SCIPcertificatePrintProblemMessage(certificate, " ");
   if( RisInfinity(upperbound) )
   {
      SCIPcertificatePrintProblemMessage(certificate, "inf");
   }
   else
   {
      SCIPcertificatePrintProblemRational(certificate, upperbound, 10);
   }
   SCIPcertificatePrintProblemMessage(certificate, "\n");
 }

/** prints RTP section for infeasibility */
void SCIPcertificatePrintRtpInfeas(
   SCIP_CERTIFICATE*     certificate         /**< certificate data structure */
   )
{
   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, "RTP infeas\n");
 }

/** prints SOL header and exact solution to certificate file */
void SCIPcertificatePrintSolex(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOLEX*           sol                 /**< primal CIP solution, may be NULL */
   )
{
   SCIP_VAR** vars;
   SCIP_Rational* solval;
   int nvars;
   int nnonz;
   int i;

   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return;

   assert(scip != NULL);

   if( sol == NULL )
   {
      SCIPcertificatePrintProblemMessage(certificate, "SOL 0\n");
      updateFilesize(certificate, 6.0);
      return;
   }

   SCIP_CALL_ABORT( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   solval = RcreateTemp(SCIPbuffer(scip));
   nnonz = 0;
   for( i = 0; i < nvars; i ++)
   {
      SCIPsolexGetVal(solval, sol, scip->set, scip->stat, vars[i]);
      if( !RisZero(solval) )
         nnonz++;
   }

   SCIPcertificatePrintProblemMessage(certificate, "SOL 1\nbest %d", nnonz);

   for( i = 0; i < nvars; i ++)
   {
      SCIPsolexGetVal(solval, sol, scip->set, scip->stat, vars[i]);
      if( !RisZero(solval) )
      {
         SCIPcertificatePrintProblemMessage(certificate, " %d ", i);
         SCIPcertificatePrintProblemRational(certificate, solval, 10);
      }
   }

   SCIPcertificatePrintProblemMessage(certificate, "\n");
}
