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

/**
 * @brief .mop file format reader
 * @author Sebastian Schenker, Timo Strunk
 *
 * Adaption of SCIP MPS reader towards MOP format with multiple objectives.
 * The input file has to follow some simple conventions
 * - It has to contain a problem in 
 * <a href="http://en.wikipedia.org/wiki/MPS_%28format%29">MPS</a> format
 * - The file extension must be <code>.mop</code>
 * - Every row marked <code>N</code> is treated as an objective
 */

#ifndef POLYSCIP_SRC_READER_MOP_H_INCLUDED
#define POLYSCIP_SRC_READER_MOP_H_INCLUDED

#include "objscip/objscip.h"

class ReaderMOP : public scip::ObjReader {
public:
    ReaderMOP(SCIP *scip)
            : scip::ObjReader(scip, "MOP Reader", "file reader for MOP file", "mop") { };

    virtual ~ReaderMOP() { };

    virtual SCIP_DECL_READERFREE(scip_free);

    virtual SCIP_DECL_READERREAD(scip_read);

};

#endif //POLYSCIP_SRC_READER_MOP_H_INCLUDED
