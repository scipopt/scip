/* Copyright (C) GAMS Development and others 2011
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Author: Stefan Vigerske
 */

#ifndef GAMSCOMPATIBILITY_H_
#define GAMSCOMPATIBILITY_H_

#ifndef GMOAPIVERSION
#error "gmomcc.h need to be included before this file to have GMOAPIVERSION defined"
#endif

#if GMOAPIVERSION < 13
#define gmoModelStat_Feasible gmoModelStat_NonOptimalIntermed 
#define gmoNameOutput         gmoNameInput
#endif

#endif /* GAMSCOMPATIBILITY_H_ */
