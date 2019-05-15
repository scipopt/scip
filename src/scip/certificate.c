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
#include "scip/pub_misc.h"
#include "scip/certificate.h"
#include "scip/struct_certificate.h"
#include "scip/solex.h"


/** checks, if value is integral */
static
SCIP_Bool mpqIsIntegral(
   const mpq_t           val                  /**< value to process */
   )
{
   return (mpz_cmp_ui(mpq_denref(val), 1) == 0);
}

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

   if( mpq_cmp(((SCIP_CERTIFICATEBOUND*)key1)->boundval, ((SCIP_CERTIFICATEBOUND*)key2)->boundval) != 0 )
      return FALSE;

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValVarbound)
{
   SCIP_CERTIFICATEBOUND* bound;

   bound = (SCIP_CERTIFICATEBOUND*)key;

   return ((((unsigned int)(10*mpq_get_d(bound->boundval))) << 22) + (bound->varindex << 2) + (unsigned int)bound->isupper);
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
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   int filenamelen;
   int bufferlen;
   char* name;
   char* extension;

   assert(certificate != NULL);
   assert(set != NULL);
   assert(set->certificate_filename != NULL);
   assert(certificate->derivationfile == NULL);

   if( set->certificate_filename[0] == '-' && set->certificate_filename[1] == '\0' )
      return SCIP_OKAY;

   filenamelen = strlen(set->certificate_filename);

   SCIP_ALLOC( BMSallocMemoryArray(&name, filenamelen) );
   BMScopyMemoryArray(name, set->certificate_filename, filenamelen);

   SCIPsplitFilename(name, NULL, NULL, NULL, &extension);

   SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      "storing certificate information in file <%s>\n", set->certificate_filename);

   if( NULL != extension && 0 == strncmp(extension, "gz", 2) )
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
   certificate->derivationfilename[bufferlen+4] = '\0';
   if( NULL != extension && 0 == strncmp(extension, "gz", 2) )
   {
      certificate->derivationfilename[bufferlen+5] = '.';
      certificate->derivationfilename[bufferlen+6] = 'g';
      certificate->derivationfilename[bufferlen+7] = 'z';

      certificate->derivationfile = SCIPfopen(certificate->derivationfilename, "wb");
   }
   else
      certificate->derivationfile = SCIPfopen(certificate->derivationfilename, "wT");

   if( certificate->file == NULL || certificate->derivationfile == NULL )
   {
      SCIPerrorMessage("error creating file <%s> and derivation file\n", set->certificate_filename);
      SCIPprintSysError(set->certificate_filename);
      return SCIP_FILECREATEERROR;
   }

   certificate->blkmem = blkmem;
   BMSfreeMemoryArray(&name);


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

   derivationfile = SCIPfopen(certificate->derivationfilename, "r");

   /* append the derivation file to the problem file */


   while( 0 != SCIPfread(buffer, sizeof(char), SCIP_MAXSTRLEN, derivationfile) )
      SCIPfwrite(buffer, sizeof(char), SCIP_MAXSTRLEN, certificate->file);

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
         /**@todo fix memory leak: mpq_clear and free all elements */
         SCIPhashtableFree(&certificate->varboundtable);
      }
      if( certificate->workbound != NULL )
      {
         mpq_clear(certificate->workbound->boundval);
         BMSfreeBlockMemory(certificate->blkmem, &certificate->workbound);
      }
      BMSfreeMemoryArrayNull(&certificate->objstring);
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
   const mpq_t*          coefs,              /**< objective function coefficients */
   int                   nvars               /**< number of variables */
   )
{
   char obj[SCIP_MAXSTRLEN - 2];
   int hashtablesize;
   int printlen;
   int nnonz;
   int buffpos;
   int i;
   int buflength;

   assert(coefs != NULL);

   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return SCIP_OKAY;

   /* create a hash table for the variable bounds */
   hashtablesize = SCIPcalcHashtableSize(1024 * nvars);
   SCIP_CALL( SCIPhashtableCreate(&certificate->varboundtable, blkmem, hashtablesize,
         hashGetKeyVarbound, hashKeyEqVarbound, hashKeyValVarbound, NULL) );

   /* create working memory for bound struct */
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, &certificate->workbound) );
   mpq_init(certificate->workbound->boundval);

   nnonz = 0;
   buflength = SCIP_MAXSTRLEN - 2;
   for( i = 0; i < nvars; i++ )
   {
      if( mpq_sgn(coefs[i]) != 0 )
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

   buffpos = 0;
   printlen = gmp_snprintf(obj, SCIP_MAXSTRLEN - 2, "%d ", nnonz);
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
      if( mpq_sgn(coefs[i]) != 0 )
      {
         printlen = gmp_snprintf(&obj[buffpos], SCIP_MAXSTRLEN - buffpos, "%s%d %Qd", (i > 0 ? " " : ""), i, coefs[i]);
         if( printlen >= SCIP_MAXSTRLEN - buffpos )
         {
            obj[buffpos] = '\0';
            SCIPcertificatePrintProblemMessage(certificate, "%s \n", obj);
         }
         else
         {
            objstring += printlen;
            leftsize -= printlen;
         }
      }
   }

   SCIPcertificatePrintProblemMessage(certificate, "OBJ min\n%s\n", certificate->objstring); 

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
   const mpq_t           val,                /**< Rational to print to the problem*/
   int                   base                /**< The base representation*/
   )
{
   char formatstr[SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->file == NULL )
     return;
   mpq_get_str(formatstr, base, val);
   SCIPcertificatePrintProblemMessage(certificate, "%s", formatstr);
}


/** prints a rational number to the proof section of the certificate file */
void SCIPcertificatePrintProofRational(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const mpq_t           val,                /**< Rational to print to the problem*/
   int                   base                /**< The base representation*/
   )
{
   char formatstr[SCIP_MAXSTRLEN];
   /* check if certificate output should be created */
   if( certificate->derivationfile == NULL )
     return;
   mpq_get_str(formatstr, base, val);
   SCIPcertificatePrintProofMessage(certificate, "%s", formatstr);
}

/** prints an integer to the problem section of the certificate file */
void SCIPcertificatePrintProblemInteger(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const mpz_t           val,                /**< Integer to print to the problem*/
   int                   base                /**< The base representation*/
   )
{
   char* formatstr;
   /* check if certificate output should be created */
   if( certificate->file == NULL )
     return;
   formatstr = mpz_get_str(NULL, base, val);
   SCIPcertificatePrintProblemMessage(certificate, "%s", formatstr);
}

/** prints an integer to the proof section of the certificate file */
void SCIPcertificatePrintProofInteger(
   SCIP_CERTIFICATE*     certificate,        /**< certificate information */
   const mpz_t           val,                /**< Integer to print to the problem*/
   int                   base                /**< The base representation*/
   )
{
  char* formatstr;
   /* check if certificate output should be created */
   if( certificate->derivationfile == NULL )
     return;
   formatstr = mpz_get_str(NULL, base, val);
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
   const mpq_t           side,               /**< left/right-hand side */
   int                   len,                /**< number of nonzeros */
   int*                  ind,                /**< index array */
   mpq_t*                val                 /**< coefficient array */
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
   const mpq_t           boundval,           /**< value of the bound */
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
   mpq_set(certificate->workbound->boundval, boundval);
   certificate->workbound->isupper = isupper;

   SCIPdebugMessage("Printing bound at line %" SCIP_LONGINT_FORMAT ": <variable %d> %s approx. %g\n",
      certificate->indexcounter, varindex, (isupper ? "<=" : ">="), mpq_get_d(boundval));

   /* bounds in the problem should be created only once */
   image = SCIPhashtableRetrieve(certificate->varboundtable, (void*)certificate->workbound);
   if( image == NULL )
   {
      SCIP_CERTIFICATEBOUND* insertbound;

      SCIP_ALLOC( BMSduplicateBlockMemory(certificate->blkmem, &insertbound, certificate->workbound) );
      mpq_init(insertbound->boundval);
      mpq_set(insertbound->boundval, boundval);
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
   const mpq_t           boundval,           /**< value of the bound */
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
   mpq_set(certificate->workbound->boundval, boundval);
   certificate->workbound->isupper = isupper;

   image = SCIPhashtableRetrieve(certificate->varboundtable, (void*)certificate->workbound);
   if( image != NULL )
   {
      SCIP_CERTIFICATEBOUND* foundbound;

      foundbound = (SCIP_CERTIFICATEBOUND*)image;

      SCIPdebugMessage("Found bound assumption at line %" SCIP_LONGINT_FORMAT ": <variable %d> %s approx. %g\n",
         foundbound->fileindex, varindex, (isupper ? "<=" : ">="), mpq_get_d(boundval));

      assert(foundbound->fileindex >= 0);
      assert(foundbound->varindex == varindex);
      assert(mpq_cmp(foundbound->boundval, boundval) == 0);
      assert(foundbound->isupper == isupper);
      return foundbound->fileindex;
   }
   else
   {
      SCIP_CERTIFICATEBOUND* insertbound;

      SCIPdebugMessage("Print bound assumption at line %" SCIP_LONGINT_FORMAT ": <variable %d> %s approx. %g\n",
         certificate->workbound->fileindex, varindex, (isupper ? "<=" : ">="), mpq_get_d(boundval));

      SCIP_ALLOC( BMSduplicateBlockMemory(certificate->blkmem, &insertbound, certificate->workbound) );
      mpq_init(insertbound->boundval);
      mpq_set(insertbound->boundval, boundval);
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

/** prints dual bound to proof section */
SCIP_Longint SCIPcertificatePrintDualbound(
   SCIP_CERTIFICATE*     certificate,        /**< certificate data structure */
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           linename,           /**< name of the unsplitting line */
   const mpq_t*          lowerbound,         /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   int                   len,                /**< number of dual multipiers */
   int*                  ind,                /**< index array */
   const mpq_t*          val                 /**< array of dual multipliers */
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
      SCIPcertificatePrintProofRational(certificate, *lowerbound, 10);
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
   if( lowerbound != NULL && SCIPgetStage(scip) == SCIP_STAGE_SOLVING && SCIPisObjIntegral(scip) && !mpqIsIntegral(*lowerbound) )
   {
      mpz_t ceilint;

      certificate->indexcounter++;

      SCIPcertificatePrintProofMessage(certificate, "R%d G ", certificate->indexcounter - 1);
      updateFilesize(certificate, 4.0 + ceil(log10(certificate->indexcounter - 1 + 1)));

      mpz_init(ceilint);
      mpz_cdiv_q(ceilint, mpq_numref(*lowerbound), mpq_denref(*lowerbound));
      SCIPcertificatePrintProofInteger(certificate, ceilint, 10);
      mpz_clear(ceilint);

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
   const mpq_t*          lowerbound,         /**< pointer to lower bound on the objective, NULL indicating infeasibility */
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
      SCIPcertificatePrintProofRational(certificate, *lowerbound, 10);
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
   const mpq_t*          lowerbound,         /**< pointer to lower bound on the objective, NULL if negative infinity */
   const mpq_t*          upperbound          /**< pointer to upper bound on the objective, NULL if positive infinity */
   )
{
   /* check if certificate output should be created */
   if( certificate->file == NULL )
      return;

   SCIPcertificatePrintProblemMessage(certificate, "RTP range ");
   if( lowerbound == NULL )
   {
      SCIPcertificatePrintProblemMessage(certificate, "-inf");
   }
   else
   {
      SCIPcertificatePrintProblemRational(certificate, *lowerbound, 10);
   }

   SCIPcertificatePrintProblemMessage(certificate, " ");
   if( upperbound == NULL )
   {
      SCIPcertificatePrintProblemMessage(certificate, "inf");
    }
   else
   {
      SCIPcertificatePrintProblemRational(certificate, *upperbound, 10);
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
   mpq_t solval;
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

   mpq_init(solval);
   nnonz = 0;
   for( i = 0; i < nvars; i ++)
   {
      SCIPsolexGetVal(sol, vars[i], solval);
      if( mpq_sgn(solval) != 0 )
         nnonz++;
   }

   SCIPcertificatePrintProblemMessage(certificate, "SOL 1\nbest %d", nnonz);

   for( i = 0; i < nvars; i ++)
   {
      SCIPsolexGetVal(sol, vars[i], solval);
      if( mpq_sgn(solval) != 0 )
      {
         SCIPcertificatePrintProblemMessage(certificate, " %d ", i);
         SCIPcertificatePrintProblemRational(certificate, solval, 10);
      }
   }

   SCIPcertificatePrintProblemMessage(certificate, "\n");
}
