/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Function                                                      */
/*   File....: sdtest.c                                                      */
/*   Name....: Special Distance Test                                         */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "grph.h"
#include "portab.h"
#include "scip/scip.h"
#include "probdata_stp.h"

double compute_node_lb(
   double*  radius,
   double*  closetermsdist,
   int*     closetermshops,
   int*     closeterms,
   int*     radiushops,
   int      termcount,
   int      nodegrad,
   int      source,
   int      probtype,
   int*     hopsbound
   )
{

   double lowerbound = 0;

   return lowerbound;
}
