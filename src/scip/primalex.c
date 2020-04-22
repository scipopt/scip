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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   primalex.c
 * @brief  methods for collecting exact primal CIP solutions and primal informations
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/visual.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/lpex.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/sol.h"
#include "scip/solex.h"
#include "scip/primal.h"
#include "scip/tree.h"
#include "scip/reopt.h"
#include "scip/disp.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/rational.h"
#include "scip/struct_sol.h"

#ifdef SCIP_WITH_GMP
#include "gmp.h"
#endif

/*
 * memory growing methods for dynamically allocated arrays
 */


