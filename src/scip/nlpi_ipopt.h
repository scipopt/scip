#pragma ident "@(#) $Id: nlpi_ipopt.h,v 1.1 2009/07/27 20:04:13 bzfviger Exp $"
/**@file    nlpi_ipopt.h
 * @brief   Ipopt NLP interface
 * @ingroup NLPINTERFACES
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_IPOPT_H__
#define __SCIP_NLPI_IPOPT_H__

#include "type_nlpi.h"
#include "scip/scip.h"

/** create solver interface for Ipopt solver
 */
extern
SCIP_RETCODE SCIPcreateNlpSolverIpopt( 
   SCIP*       scip,  /**< central scip datastructure */
   SCIP_NLPI** nlpi   /**< pointer to buffer for nlpi address */
);

#endif /* __SCIP_NLPI_IPOPT_H__ */
