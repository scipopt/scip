// Copyright (C) GAMS Development and others 2011
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Stefan Vigerske

#include <cstdlib>
#include <cassert>

#include "GamsSolverC_tpl.h"

/* check that corresponding solver class is defined */
#ifndef GAMSSOLVERC_CLASS
#error You need to define GAMSSOLVERC_CLASS
#endif

/* check whether constructor arguments are given */
#ifndef GAMSSOLVERC_CONSTRARGS
#define GAMSSOLVERC_CONSTRARGS
#endif

#ifdef GAMS_BUILD

DllExport int STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,Create)(void** Cptr, char* msgBuf, int msgBufLen)
{
   assert(Cptr != NULL);
   assert(msgBufLen > 0);
   assert(msgBuf != NULL);

   *Cptr = (void*) new GAMSSOLVERC_CLASS(GAMSSOLVERC_CONSTRARGS);
   msgBuf[0] = 0;

   return 1;
}

DllExport void STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,Free)(void** Cptr)
{
   assert(Cptr != NULL);
   delete (GAMSSOLVERC_CLASS*)*Cptr;
   *Cptr = NULL;
}

DllExport int STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,CallSolver)(void* Cptr)
{
   assert(Cptr != NULL);
   return ((GAMSSOLVERC_CLASS*)Cptr)->callSolver();
}

DllExport int STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,HaveModifyProblem)(void* Cptr)
{
   return GAMSSOLVERC_CLASS::haveModifyProblem();
}

DllExport int STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,ModifyProblem)(void* Cptr)
{
   assert(Cptr != NULL);
   return ((GAMSSOLVERC_CLASS*)Cptr)->modifyProblem();
}

DllExport int STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,ReadyAPI)(void* Cptr, gmoHandle_t Gptr, optHandle_t Optr)
{
   assert(Cptr != NULL);
   assert(Gptr != NULL);
   char msg[256];
   if( !gmoGetReady(msg, sizeof(msg)) )
      return 1;
   if( !gevGetReady(msg, sizeof(msg)) )
      return 1;
   return ((GAMSSOLVERC_CLASS*)Cptr)->readyAPI(Gptr, Optr);
}

#else

DllExport void STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,XCreate)(void** Cptr)
{
   assert(Cptr != NULL);
   *Cptr = (void*) new GAMSSOLVERC_CLASS(GAMSSOLVERC_CONSTRARGS);
}

DllExport int STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,create)(void** Cptr, char* msgBuf, int msgBufLen)
{
   assert(Cptr != NULL);
   assert(msgBufLen > 0);
   assert(msgBuf != NULL);

   *Cptr = (void*) new GAMSSOLVERC_CLASS(GAMSSOLVERC_CONSTRARGS);
   msgBuf[0] = 0;

   return 1;
}

DllExport void STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,XFree)(void** Cptr)
{
   assert(Cptr != NULL);
   delete (GAMSSOLVERC_CLASS*)*Cptr;
   *Cptr = NULL;
   gmoLibraryUnload();
   gevLibraryUnload();
}

DllExport int STDCALL GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,free)(void** Cptr)
{
   GAMSSOLVERC_CONCAT(GAMSSOLVERC_ID,XFree)(Cptr);
   return 1;
}

/* comp returns the compatibility mode:
   0: client is too old for the DLL, no compatibility
   1: client version and DLL version are the same, full compatibility
   2: client is older than DLL, but defined as compatible, backward compatibility
   3: client is newer than DLL, forward compatibility
   FIXME: for now, we just claim full compatibility
*/
DllExport int STDCALL GAMSSOLVERC_CONCAT3(C__,GAMSSOLVERC_ID,XAPIVersion)(int api, char* Msg, int* comp)
{
  *comp = 1;
  return 1;
}

DllExport int STDCALL GAMSSOLVERC_CONCAT3(D__,GAMSSOLVERC_ID,XAPIVersion)(int api, char* Msg, int* comp)
{
  *comp = 1;
  return 1;
}

DllExport int STDCALL GAMSSOLVERC_CONCAT3(C__,GAMSSOLVERC_ID,XCheck)(const char* funcn, int ClNrArg, int Clsign[], char* Msg)
{
   return 1;
}

DllExport int STDCALL GAMSSOLVERC_CONCAT3(D__,GAMSSOLVERC_ID,XCheck)(const char* funcn, int ClNrArg, int Clsign[], char* Msg)
{
   return 1;
}

DllExport int STDCALL GAMSSOLVERC_CONCAT3(C__,GAMSSOLVERC_ID,CallSolver)(void* Cptr)
{
   assert(Cptr != NULL);
   return ((GAMSSOLVERC_CLASS*)Cptr)->callSolver();
}

DllExport int STDCALL GAMSSOLVERC_CONCAT3(C__,GAMSSOLVERC_ID,HaveModifyProblem)(void* Cptr)
{
   return GAMSSOLVERC_CLASS::haveModifyProblem();
}

DllExport int STDCALL GAMSSOLVERC_CONCAT3(C__,GAMSSOLVERC_ID,ModifyProblem)(void* Cptr)
{
   assert(Cptr != NULL);
   return ((GAMSSOLVERC_CLASS*)Cptr)->modifyProblem();
}

DllExport int STDCALL GAMSSOLVERC_CONCAT3(C__,GAMSSOLVERC_ID,ReadyAPI)(void* Cptr, gmoHandle_t Gptr, optHandle_t Optr)
{
   assert(Cptr != NULL);
   assert(Gptr != NULL);
   char msg[256];
   if( !gmoGetReady(msg, sizeof(msg)) )
      return 1;
   if( !gevGetReady(msg, sizeof(msg)) )
      return 1;
   return ((GAMSSOLVERC_CLASS*)Cptr)->readyAPI(Gptr, Optr);
}

#endif

#undef GAMSSOLVERC_ID
#undef GAMSSOLVERC_CLASS
#undef GAMSSOLVERC_CONSTRARGS
