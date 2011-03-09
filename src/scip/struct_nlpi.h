/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_nlpi.h
 * @brief  data definitions for an NLP solver interface
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_NLPI_H__
#define __SCIP_STRUCT_NLPI_H__

#include "scip/scip.h"
#include "scip/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** NLP interface data */
struct SCIP_Nlpi
{
   char*                           name;                        /**< name of NLP solver */
   char*                           description;                 /**< description of NLP solver */
   int                             priority;                    /**< priority of NLP interface */
   SCIP_DECL_NLPIINIT              ((*nlpiinit));               /**< initialize NLPI user data */
   SCIP_DECL_NLPIFREE              ((*nlpifree));               /**< free NLPI user data */
   SCIP_DECL_NLPIADDVARS           ((*nlpiaddvars));            /**< add variables */
   SCIP_DECL_NLPIADDCONSTRAINTS    ((*nlpiaddconstraints));     /**< add constraints */
   SCIP_DECL_NLPISETOBJECTIVE      ((*nlpisetobjective));       /**< set objective */
   SCIP_DECL_NLPICHGVARBOUNDS      ((*nlpichgvarbounds));       /**< change variable bounds */
   SCIP_DECL_NLPICHGCONSBOUNDS     ((*nlpichgconsbounds));      /**< change constraint bounds */
   SCIP_DECL_NLPIDELVARSET         ((*nlpidelvarset));          /**< delete a set of variables */
   SCIP_DECL_NLPIDELCONSSET        ((*nlpidelconsset));         /**< delete a set of constraints */
   SCIP_DECL_NLPICHGLINEARCOEFS    ((*nlpichglinearcoefs));     /**< change one coefficient in linear part */
   SCIP_DECL_NLPICHGQUADCOEFS      ((*nlpichgquadcoefs));       /**< change one coefficient in quadratic part */
   SCIP_DECL_NLPICHGNONLINCOEF     ((*nlpichgnonlincoef));      /**< change one parameter in nonlinear expressions */
   SCIP_DECL_NLPISETINITIALGUESS   ((*nlpisetinitialguess));    /**< set initial guess for primal variables */
   SCIP_DECL_NLPISOLVE             ((*nlpisolve));              /**< solve NLP */
   SCIP_DECL_NLPIGETSOLSTAT        ((*nlpigetsolstat));         /**< get solution status */
   SCIP_DECL_NLPIGETTERMSTAT       ((*nlpigettermstat));        /**< get termination status */
   SCIP_DECL_NLPIGETSOLUTION       ((*nlpigetsolution));        /**< get solution */
   SCIP_DECL_NLPIGETSTATISTICS     ((*nlpigetstatistics));      /**< get solve statistics */
   SCIP_DECL_NLPIGETWARMSTARTSIZE  ((*nlpigetwarmstartsize));   /**< get size for warmstart object buffer */
   SCIP_DECL_NLPIGETWARMSTARTMEMO  ((*nlpigetwarmstartmemo));   /**< get warmstart object */
   SCIP_DECL_NLPISETWARMSTARTMEMO  ((*nlpisetwarmstartmemo));   /**< set warmstart object */
   SCIP_DECL_NLPIGETSOLVERPOINTER  ((*nlpigetsolverpointer));   /**< get solver pointer */
   SCIP_DECL_NLPIGETINTPAR         ((*nlpigetintpar));          /**< get value of integer parameter */
   SCIP_DECL_NLPISETINTPAR         ((*nlpisetintpar));          /**< set value of integer parameter */
   SCIP_DECL_NLPIGETREALPAR        ((*nlpigetrealpar));         /**< get value of floating point parameter */
   SCIP_DECL_NLPISETREALPAR        ((*nlpisetrealpar));         /**< set value of floating point parameter */
   SCIP_NLPIDATA*                  nlpidata;                    /**< NLP interface local data */
};

/** Statistics from an NLP solve */
struct SCIP_NlpStatistics
{
   int       niterations;   /**< number of iterations the NLP solver spend in the last solve command */
   SCIP_Real totaltime;     /**< total time in CPU sections the NLP solver spend in the last solve command */ 
};

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_STRUCT_NLPI_H__ */
