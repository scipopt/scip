/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*        This file is part of the program PolySCIP                          */
/*                                                                           */
/*    Copyright (C) 2012-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  PolySCIP is distributed under the terms of the ZIB Academic License.     */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with PolySCIP; see the file LICENCE.                               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   ReaderMOP.h
 * @brief  mop reader
 * @author Sebastian Schenker, Timo Strunk
 *
 * The mop-reader adapts the mps-reader of SCIP by the functionality to read multiple objectives.
 *
 * The input file has to follow some simple conventions
 * - It has to contain a problem in 
 * <a href="http://en.wikipedia.org/wiki/MPS_%28format%29">MPS</a> format
 * - The file extension must be <code>.mop</code>
 * - Every row marked <code>N</code> is treated as an objective
 *
 */

#ifndef _GUARD_READER_MOP_H_
#define _GUARD_READER_MOP_H_

#include "objscip/objscip.h"

class ReaderMOP : public scip::ObjReader {
 public:

 ReaderMOP(SCIP* scip)
   : scip::ObjReader(scip, "MOP Reader", "file reader for MOP file", "mop")
    {}  
  
  virtual ~ReaderMOP() {}
  
  virtual SCIP_DECL_READERFREE(scip_free);
  
  virtual SCIP_DECL_READERREAD(scip_read);

};

#endif
