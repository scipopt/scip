/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_fileio.h
 * @ingroup PUBLICMETHODS
 * @brief  wrapper functions to map file i/o to standard or zlib file i/o
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_FILEIO_H__
#define __SCIP_PUB_FILEIO_H__

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_File SCIP_FILE;          /**< file data structure */

extern SCIP_FILE* SCIPfopen(const char *path, const char *mode);
extern SCIP_FILE* SCIPfdopen(int fildes, const char *mode);
extern size_t SCIPfread(void *ptr, size_t size, size_t nmemb, SCIP_FILE *stream);
extern size_t SCIPfwrite(const void *ptr, size_t size, size_t nmemb, SCIP_FILE *stream);
extern int SCIPfprintf(SCIP_FILE *stream, const char *format, ...);
extern int SCIPfputc(int c, SCIP_FILE *stream);
extern int SCIPfputs(const char *s, SCIP_FILE *stream);
extern int SCIPfgetc(SCIP_FILE *stream);
extern char* SCIPfgets(char *s, int size, SCIP_FILE *stream);
extern int SCIPfflush(SCIP_FILE *stream);
extern int SCIPfseek(SCIP_FILE *stream, long offset, int whence);
extern void SCIPrewind(SCIP_FILE *stream);
extern long SCIPftell(SCIP_FILE *stream);
extern int SCIPfeof(SCIP_FILE *stream);
extern int SCIPfclose(SCIP_FILE *fp);

#ifdef __cplusplus
}
#endif

#endif
