/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_fileio.h,v 1.1 2005/07/15 17:20:15 bzfpfend Exp $"

/**@file   pub_fileio.h
 * @brief  wrapper header to map file i/o to standard or zlib file i/o
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_FILEIO_H__
#define __SCIP_PUB_FILEIO_H__


#ifdef WITH_ZLIB

/* file i/o using zlib */
#include <zlib.h>

#define SCIPFILE                             gzFile
#define SCIPfopen                            gzopen
#define SCIPfdopen                           gzdopen
#define SCIPfread(ptr,size,nmemb,stream)     gzread(stream,ptr,(size)*(nmemb))
#define SCIPfwrite(ptr,size,nmemb,stream)    gzwrite(stream,ptr,(size)*(nmemb))
#define SCIPfprintf                          gzprintf
#define SCIPfputs(s,stream)                  gzputs(stream,s)
#define SCIPfgets(s,size,stream)             gzgets(stream,s,size)
#define SCIPfputc(c,stream)                  gzputc(stream,c)
#define SCIPfgetc                            gzgetc
#define SCIPfflush(stream)                   gzflush(stream,Z_SYNC_FLUSH)
#define SCIPfseek                            gzseek
#define SCIPrewind                           gzrewind
#define SCIPftell                            gztell
#define SCIPfeof                             gzeof
#define SCIPfclose                           gzclose

#else

/* file i/o using standard i/o */
#include <stdio.h>

#define SCIPFILE      FILE
#define SCIPfopen     fopen
#define SCIPfdopen    fdopen
#define SCIPfread     fread
#define SCIPfwrite    fwrite
#define SCIPfprintf   fprintf
#define SCIPfputs     fputs
#define SCIPfgets     fgets
#define SCIPfputc     fputc
#define SCIPfgetc     fgetc
#define SCIPfflush    fflush
#define SCIPfseek     fseek
#define SCIPrewind    rewind
#define SCIPftell     ftell
#define SCIPfeof      feof
#define SCIPfclose    fclose

#endif


#endif
