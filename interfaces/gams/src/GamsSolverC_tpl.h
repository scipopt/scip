// Copyright (C) GAMS Development and others 2011
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Stefan Vigerske

/* check that three-letter id of solver is defined */
#ifndef GAMSSOLVERC_ID
#error You need to define GAMSSOLVERC_ID
#endif

#if defined(_WIN32)
# if ! defined(STDCALL)
#  define STDCALL   __stdcall
# endif
# if ! defined(DllExport)
#  define DllExport __declspec( dllexport )
# endif
#else
# if ! defined(STDCALL)
#  define STDCALL
# endif
# if ! defined(DllExport)
#  define DllExport
# endif
#endif

/* we need to redirect the CONCAT macro to another one to ensure that macro arguments are substituted too */

#ifndef GAMSSOLVERC_CONCAT_
#define GAMSSOLVERC_CONCAT_(a, b) a ## b
#endif

#ifndef GAMSSOLVERC_CONCAT
#define GAMSSOLVERC_CONCAT(a, b) GAMSSOLVERC_CONCAT_(a, b)
#endif

#ifndef GAMSSOLVERC_CONCAT3_
#define GAMSSOLVERC_CONCAT3_(a, b, c) a ## b ## c
#endif

#ifndef GAMSSOLVERC_CONCAT3
#define GAMSSOLVERC_CONCAT3(a, b, c) GAMSSOLVERC_CONCAT3_(a, b, c)
#endif

#if defined(__cplusplus)
extern "C" {
#endif

#ifdef GAMS_BUILD

typedef void GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,Rec_t);

DllExport int  STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,Create)(void** Cptr, char* msgBuf, int msgBufLen);

DllExport void STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,Free)(void** Cptr);

DllExport int  STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,CallSolver)(void* Cptr);

DllExport int  STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,HaveModifyProblem)(void* Cptr);

DllExport int  STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,ModifyProblem)(void* Cptr);

DllExport int  STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,ReadyAPI)(void* Cptr, struct gmoRec* Gptr, struct optRec* Optr);

#else

DllExport void STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,XCreate)(void** Cptr);

DllExport int  STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,create)(void** Cptr, char* msgBuf, int msgBufLen);

DllExport void STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,XFree)(void** Cptr);

DllExport int  STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,free)(void** Cptr);

DllExport int  STDCALL GAMSSOLVERC_CONCAT3(C__,GAMSSOLVERC_ID,XAPIVersion)(int api, char* Msg, int* comp);

DllExport int  STDCALL GAMSSOLVERC_CONCAT3(D__,GAMSSOLVERC_ID,XAPIVersion)(int api, char* Msg, int* comp);

DllExport int  STDCALL GAMSSOLVERC_CONCAT3(C__,GAMSSOLVERC_ID,XCheck)(const char* funcn, int ClNrArg, int Clsign[], char* Msg);

DllExport int  STDCALL GAMSSOLVERC_CONCAT3(D__,GAMSSOLVERC_ID,XCheck)(const char* funcn, int ClNrArg, int Clsign[], char* Msg);

DllExport int  STDCALL GAMSSOLVERC_CONCAT3(C__,GAMSSOLVERC_ID,CallSolver)(void* Cptr);

DllExport int  STDCALL GAMSSOLVERC_CONCAT3(C__,GAMSSOLVERC_ID,HaveModifyProblem)(void* Cptr);

DllExport int  STDCALL GAMSSOLVERC_CONCAT3(C__,GAMSSOLVERC_ID,ModifyProblem)(void* Cptr);

DllExport int  STDCALL GAMSSOLVERC_CONCAT3(C__,GAMSSOLVERC_ID,ReadyAPI)(void* Cptr, struct gmoRec* Gptr, struct optRec* Optr);

#endif

#if defined(__cplusplus)
}
#endif
