/* Copyright (C) 2008-2011 GAMS Development and others
   All Rights Reserved.
   This code is published under the Eclipse Public License.

   Author: Stefan Vigerske

   WARNING:
   The NL instruction (codes, format, semantics, etc.) are subject to change.
   Thus, there is no warranty that a code based on this file is going to work fine
   after changes in the GAMS I/O libraries as provided by GAMS.
   GAMS is NOT obligated to issue a warning when it changes the NL instructions.
 */

#ifndef GAMSNLINSTR_H_
#define GAMSNLINSTR_H_

#ifdef __cplusplus
extern "C" {
#endif

/** GAMS nonlinear expressions opcodes */
typedef enum
{
   nlNoOp     =  0, /* no operation */
   nlPushV    =  1, /* push variable */
   nlPushI    =  2, /* push immediate (constant) */
   nlStore    =  3, /* store row */
   nlAdd      =  4, /* add */
   nlAddV     =  5, /* add variable */
   nlAddI     =  6, /* add immediate */
   nlSub      =  7, /* minus */
   nlSubV     =  8, /* subtract variable */
   nlSubI     =  9, /* subtract immediate */
   nlMul      = 10, /* multiply */
   nlMulV     = 11, /* multiply variable */
   nlMulI     = 12, /* multiply immediate */
   nlDiv      = 13, /* divide */
   nlDivV     = 14, /* divide variable */
   nlDivI     = 15, /* divide immediate */
   nlUMin     = 16, /* unary minus */
   nlUMinV    = 17, /* unary minus variable */
   nlHeader   = 18, /* header */
   nlEnd      = 19, /* end of instruction list */
   nlCallArg1 = 20,
   nlCallArg2 = 21,
   nlCallArgN = 22,
   nlFuncArgN = 23,
   nlMulIAdd  = 24,
   nlPushZero = 25,
   nlChk      = 26,
   nlAddO     = 27,
   nlPushO    = 28,
   nlInvoc    = 29,
   nlStackIn  = 30,
   MAXINS     = 31
} GamsOpCode;

/** GAMS nonlinear expression opcode names as strings */
static const char* GamsOpCodeName[MAXINS] =
{
   "nlNoOp",
   "nlPushV",
   "nlPushI",
   "nlStore",
   "nlAdd",
   "nlAddV",
   "nlAddI",
   "nlSub",
   "nlSubV",
   "nlSubI",
   "nlMul",
   "nlMulV",
   "nlMulI",
   "nlDiv",
   "nlDivV",
   "nlDivI",
   "nlUMin",
   "nlUMinV",
   "nlHeader",
   "nlEnd",
   "nlCallArg1",
   "nlCallArg2",
   "nlCallArgN",
   "nlFuncArgN",
   "nlMulIAdd",
   "nlPushZero",
   "nlChk",
   "nlAddO",
   "nlPushO",
   "nlInvoc",
   "nlStackIn"
};

/** GAMS nonlinear function codes */
typedef enum
{
   fnmapval=0,fnceil,fnfloor,fnround,
   fnmod,fntrunc,fnsign,fnmin,
   fnmax,fnsqr,fnexp,fnlog,
   fnlog10,fnsqrt,fnabs,fncos,
   fnsin,fnarctan,fnerrf,fndunfm,
   fndnorm,fnpower,fnjdate,fnjtime,
   fnjstart,fnjnow,fnerror,fngyear,
   fngmonth,fngday,fngdow,fngleap,
   fnghour,fngminute,fngsecond,
   fncurseed,fntimest,fntimeco,
   fntimeex,fntimecl,fnfrac,fnerrorl,
   fnheaps,fnfact,fnunfmi,fnpi,
   fnncpf,fnncpcm,fnentropy,fnsigmoid,
   fnlog2,fnboolnot,fnbooland,
   fnboolor,fnboolxor,fnboolimp,
   fnbooleqv,fnrelopeq,fnrelopgt,
   fnrelopge,fnreloplt,fnrelople,
   fnrelopne,fnifthen,fnrpower,
   fnedist,fndiv,fndiv0,fnsllog10,
   fnsqlog10,fnslexp,fnsqexp,fnslrec,
   fnsqrec,fncvpower,fnvcpower,
   fncentropy,fngmillisec,fnmaxerror,
   fntimeel,fngamma,fnloggamma,fnbeta,
   fnlogbeta,fngammareg,fnbetareg,
   fnsinh,fncosh,fntanh,fnmathlastrc,
   fnmathlastec,fnmathoval,fnsignpower,
   fnhandle,fnncpvusin,fnncpvupow,
   fnbinomial,fnrehandle,fngamsver,
   fndelhandle,fntan,fnarccos,
   fnarcsin,fnarctan2,fnsleep,fnheapf,
   fncohandle,fngamsrel,fnpoly,
   fnlicensestatus,fnlicenselevel,fnheaplimit,
   fndummy
} GamsFuncCode;

static const char* GamsFuncCodeName[fndummy+1] =
{
   "mapval","ceil","floor","round",
   "mod","trunc","sign","min",
   "max","sqr","exp","log",
   "log10","sqrt","abs","cos",
   "sin","arctan","errf","dunfm",
   "dnorm","power","jdate","jtime",
   "jstart","jnow","error","gyear",
   "gmonth","gday","gdow","gleap",
   "ghour","gminute","gsecond",
   "curseed","timest","timeco",
   "timeex","timecl","frac","errorl",
   "heaps","fact","unfmi","pi",
   "ncpf","ncpcm","entropy","sigmoid",
   "log2","boolnot","booland",
   "boolor","boolxor","boolimp",
   "booleqv","relopeq","relopgt",
   "relopge","reloplt","relople",
   "relopne","ifthen","rpower",
   "edist","div","div0","sllog10",
   "sqlog10","slexp","sqexp","slrec",
   "sqrec","cvpower","vcpower",
   "centropy","gmillisec","maxerror",
   "timeel","gamma","loggamma","beta",
   "logbeta","gammareg","betareg",
   "sinh","cosh","tanh","mathlastrc",
   "mathlastec","mathoval","signpower",
   "handle","ncpvusin","ncpvupow",
   "binomial","rehandle","gamsver",
   "delhandle","tan","arccos",
   "arcsin","arctan2","sleep","heapf",
   "cohandle","gamsrel","poly",
   "licensestatus","licenselevel","heaplimit",
   "dummy"
};


#ifdef __cplusplus
}
#endif

#endif /*GAMSNLINSTR_H_*/
