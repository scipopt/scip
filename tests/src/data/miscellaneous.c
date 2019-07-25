/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  this file is part of the program and library             */
/*         scip --- solving constraint integer programs                      */
/*                                                                           */
/*    copyright (c) 2002-2019 konrad-zuse-zentrum                            */
/*                            fuer informationstechnik berlin                */
/*                                                                           */
/*  scip is distributed under the terms of the zib academic license.         */
/*                                                                           */
/*  you should have received a copy of the zib academic license              */
/*  along with scip; see the file copying. if not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   arrays.c
 * @brief  unittest for arrays in scip_datastructures.c
 * @author Merlin Viernickel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"
#include "scip/misc.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_QUEUE* queue;


/* creates scip and problem */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* create queue */
   SCIP_CALL( SCIPqueueCreate(&queue, 8, 2) );
}

static
void teardown(void)
{
   /* test queue */
   SCIP_CALL( SCIPqueueInsert(queue, (void*) scip) );
   SCIP_CALL( SCIPqueueInsertUInt(queue, 2) );
   assert((void*) scip == SCIPqueueFirst(queue));
   assert(2 == SCIPqueueFirstUInt(queue));
   assert(2 == SCIPqueueNElems(queue));

   assert(2 == SCIPqueueRemoveUInt(queue));
   SCIPqueueRemove(queue);
   assert(SCIPqueueIsEmpty(queue));

   SCIPqueueFree(&queue);

   SCIP_CALL( SCIPfree(&scip) );
}
