/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_cons.h,v 1.29 2005/02/25 14:27:09 bzfpfend Exp $"

/**@file   struct_cons.h
 * @brief  datastructures for constraints and constraint handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_CONS_H__
#define __STRUCT_CONS_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_cons.h"



/** constraint data structure */
struct Cons
{
   Real             age;                /**< age of constraint: number of successive times, the constraint was irrelevant */
   char*            name;               /**< name of the constraint */
   CONSHDLR*        conshdlr;           /**< constraint handler for this constraint */
   CONSDATA*        consdata;           /**< data for this specific constraint */
   CONS*            transorigcons;      /**< for original constraints: associated transformed constraint or NULL,
                                         *   for transformed constraints: associated original constraint or NULL */
   CONSSETCHG*      addconssetchg;      /**< constraint change that added constraint to current subproblem, or NULL if
                                         *   constraint is from global problem */
   int              addarraypos;        /**< position of constraint in the conssetchg's/prob's addedconss/conss array */
   int              consspos;           /**< position of constraint in the handler's conss array */
   int              sepaconsspos;       /**< position of constraint in the handler's sepaconss array */
   int              enfoconsspos;       /**< position of constraint in the handler's enfoconss array */
   int              checkconsspos;      /**< position of constraint in the handler's checkconss array */
   int              propconsspos;       /**< position of constraint in the handler's propconss array */
   int              nuses;              /**< number of times, this constraint is referenced */
   int              nlockspos;          /**< number of times, the constraint locked rounding of its variables */
   int              nlocksneg;          /**< number of times, the constraint locked vars for the constraint's negation */
   int              activedepth;        /**< depth level of constraint activation (-2: inactive, -1: problem constraint) */
   unsigned int     initial:1;          /**< TRUE iff LP relaxation of constraint should be in initial LP, if possible */
   unsigned int     separate:1;         /**< TRUE iff constraint should be separated during LP processing */
   unsigned int     enforce:1;          /**< TRUE iff constraint should be enforced during node processing */
   unsigned int     check:1;            /**< TRUE iff constraint should be checked for feasibility */
   unsigned int     propagate:1;        /**< TRUE iff constraint should be propagated during node processing */
   unsigned int     propenabled:1;      /**< TRUE iff constraint should be propagated in the next propagation call */
   unsigned int     local:1;            /**< TRUE iff constraint is only valid locally */
   unsigned int     modifiable:1;       /**< TRUE iff constraint is modifiable (subject to column generation) */
   unsigned int     removeable:1;       /**< TRUE iff constraint should be removed from the LP due to aging or cleanup */
   unsigned int     original:1;         /**< TRUE iff constraint belongs to original problem */
   unsigned int     active:1;           /**< TRUE iff constraint is active in the current node */
   unsigned int     enabled:1;          /**< TRUE iff constraint is enforced, separated, and propagated in current node */
   unsigned int     obsolete:1;         /**< TRUE iff constraint is too seldomly used and therefore obsolete */
   unsigned int     deleted:1;          /**< TRUE iff constraint was globally deleted */
   unsigned int     update:1;           /**< TRUE iff constraint has to be updated in update phase */
   unsigned int     updateactivate:1;   /**< TRUE iff constraint has to be activated in update phase */
   unsigned int     updatedeactivate:1; /**< TRUE iff constraint has to be deactivated in update phase */
   unsigned int     updateenable:1;     /**< TRUE iff constraint has to be enabled in update phase */
   unsigned int     updatedisable:1;    /**< TRUE iff constraint has to be disabled in update phase */
   unsigned int     updatepropenable:1; /**< TRUE iff constraint's propagation has to be enabled in update phase */
   unsigned int     updatepropdisable:1;/**< TRUE iff constraint's propagation has to be disabled in update phase */
   unsigned int     updatedelete:1;     /**< TRUE iff constraint has to be deleted in update phase */
   unsigned int     updateobsolete:1;   /**< TRUE iff obsolete status of constraint has to be updated in update phase */
};

/** tracks additions and removals of the set of active constraints */
struct ConsSetChg
{
   CONS**           addedconss;         /**< constraints added to the set of active constraints */
   CONS**           disabledconss;      /**< constraints disabled in the set of active constraints */
   int              addedconsssize;     /**< size of added constraints array */
   int              naddedconss;        /**< number of added constraints */
   int              disabledconsssize;  /**< size of disabled constraints array */
   int              ndisabledconss;     /**< number of disabled constraints */
};

/** constraint handler */
struct Conshdlr
{
   Longint          nsepacalls;         /**< number of times, the separator was called */
   Longint          nenfolpcalls;       /**< number of times, the LP enforcer was called */
   Longint          nenfopscalls;       /**< number of times, the pseudo enforcer was called */
   Longint          npropcalls;         /**< number of times, the propagator was called */
   Longint          ncutoffs;           /**< number of cutoffs found so far by this constraint handler */
   Longint          ncutsfound;         /**< number of cuts found by this constraint handler */
   Longint          nconssfound;        /**< number of additional constraints added by this constraint handler */
   Longint          ndomredsfound;      /**< number of domain reductions found so far by this constraint handler */
   Longint          nchildren;          /**< number of children the constraint handler created during branching */
   Longint          lastpropdomchgcount;/**< last bound change number, where the domain propagation was called */
   Longint          lastenfolpdomchgcount;/**< last bound change number, where the LP enforcement was called */
   Longint          lastenfopsdomchgcount;/**< last bound change number, where the pseudo enforcement was called */
   char*            name;               /**< name of constraint handler */
   char*            desc;               /**< description of constraint handler */
   DECL_CONSFREE    ((*consfree));      /**< destructor of constraint handler */
   DECL_CONSINIT    ((*consinit));      /**< initialize constraint handler */
   DECL_CONSEXIT    ((*consexit));      /**< deinitialize constraint handler */
   DECL_CONSINITPRE ((*consinitpre));   /**< presolving initialization method of constraint handler */
   DECL_CONSEXITPRE ((*consexitpre));   /**< presolving deinitialization method of constraint handler */
   DECL_CONSINITSOL ((*consinitsol));   /**< solving process initialization method of constraint handler */
   DECL_CONSEXITSOL ((*consexitsol));   /**< solving process deinitialization method of constraint handler */
   DECL_CONSDELETE  ((*consdelete));    /**< free specific constraint data */
   DECL_CONSTRANS   ((*constrans));     /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSINITLP  ((*consinitlp));    /**< initialize LP with relaxations of "initial" constraints */
   DECL_CONSSEPA    ((*conssepa));      /**< separate cutting planes */
   DECL_CONSENFOLP  ((*consenfolp));    /**< enforcing constraints for LP solutions */
   DECL_CONSENFOPS  ((*consenfops));    /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHECK   ((*conscheck));     /**< check feasibility of primal solution */
   DECL_CONSPROP    ((*consprop));      /**< propagate variable domains */
   DECL_CONSPRESOL  ((*conspresol));    /**< presolving method */
   DECL_CONSRESPROP ((*consresprop));   /**< propagation conflict resolving method */
   DECL_CONSLOCK    ((*conslock));      /**< variable rounding lock method */
   DECL_CONSACTIVE  ((*consactive));    /**< activation notification method */
   DECL_CONSDEACTIVE((*consdeactive));  /**< deactivation notification method */
   DECL_CONSENABLE  ((*consenable));    /**< enabling notification method */
   DECL_CONSDISABLE ((*consdisable));   /**< disabling notification method */
   DECL_CONSPRINT   ((*consprint));     /**< constraint display method */
   CONSHDLRDATA*    conshdlrdata;       /**< constraint handler data */
   CONS**           conss;              /**< array with all transformed constraints, active ones preceed incative ones */
   CONS**           sepaconss;          /**< array with active constraints that must be separated during LP processing */
   CONS**           enfoconss;          /**< array with active constraints that must be enforced during node processing */
   CONS**           checkconss;         /**< array with active constraints that must be checked for feasibility */
   CONS**           propconss;          /**< array with active constraints that must be propagated during node processing */
   CONS**           updateconss;        /**< array with constraints that changed and have to be update in the handler */
   CLOCK*           presoltime;         /**< time used for presolving of this constraint handler */
   CLOCK*           sepatime;           /**< time used for separation of this constraint handler */
   CLOCK*           enfolptime;         /**< time used for LP enforcement of this constraint handler */
   CLOCK*           enfopstime;         /**< time used for pseudo enforcement of this constraint handler */
   CLOCK*           proptime;           /**< time used for propagation of this constraint handler */
   int              sepapriority;       /**< priority of the constraint handler for separation */
   int              enfopriority;       /**< priority of the constraint handler for constraint enforcing */
   int              checkpriority;      /**< priority of the constraint handler for checking infeasibility */
   int              sepafreq;           /**< frequency for separating cuts; zero means to separate only in the root node */
   int              propfreq;           /**< frequency for propagating domains; zero means only preprocessing propagation */
   int              eagerfreq;          /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
   int              maxprerounds;       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   int              consssize;          /**< size of conss array */
   int              nconss;             /**< total number of constraints */
   int              nactiveconss;       /**< total number of active constraints */
   int              maxnactiveconss;    /**< maximal number of active constraints existing at the same time */
   int              startnactiveconss;  /**< number of active constraints existing when problem solving started */
   int              sepaconsssize;      /**< size of sepaconss array */
   int              nsepaconss;         /**< number of active constraints that may be separated during LP processing */
   int              nusefulsepaconss;   /**< number of non-obsolete active constraints that should be separated */
   int              enfoconsssize;      /**< size of enfoconss array */
   int              nenfoconss;         /**< number of active constraints that must be enforced during node processing */
   int              nusefulenfoconss;   /**< number of non-obsolete active constraints that must be enforced */
   int              checkconsssize;     /**< size of checkconss array */
   int              ncheckconss;        /**< number of active constraints that must be checked for feasibility */
   int              nusefulcheckconss;  /**< number of non-obsolete active constraints that must be checked */
   int              propconsssize;      /**< size of propconss array */
   int              npropconss;         /**< number of active constraints that may be propagated during node processing */
   int              nusefulpropconss;   /**< number of non-obsolete active constraints that should be propagated */
   int              updateconsssize;    /**< size of updateconss array */
   int              nupdateconss;       /**< number of update constraints */
   int              nenabledconss;      /**< total number of enabled constraints of the handler */
   int              lastsepalpcount;    /**< last LP number, where the separations was called */
   int              lastenfolplpcount;  /**< last LP number, where the LP enforcement was called */
   int              lastnusefulpropconss;/**< number of already propagated useful constraints on current domains */
   int              lastnusefulsepaconss;/**< number of already separated useful constraints on current solution */
   int              lastnusefulenfoconss;/**< number of already enforced useful constraints on current solution */
   int              lastnfixedvars;     /**< number of variables fixed before the last call to the presolver */
   int              lastnaggrvars;      /**< number of variables aggregated before the last call to the presolver */
   int              lastnchgvartypes;   /**< number of variable type changes before the last call to the presolver */
   int              lastnchgbds;        /**< number of variable bounds tightend before the last call to the presolver */
   int              lastnaddholes;      /**< number of domain holes added before the last call to the presolver */
   int              lastndelconss;      /**< number of deleted constraints before the last call to the presolver */
   int              lastnupgdconss;     /**< number of upgraded constraints before the last call to the presolver */
   int              lastnchgcoefs;      /**< number of changed coefficients before the last call to the presolver */
   int              lastnchgsides;      /**< number of changed left or right hand sides before the last call */
   int              nfixedvars;         /**< total number of variables fixed by this presolver */
   int              naggrvars;          /**< total number of variables aggregated by this presolver */
   int              nchgvartypes;       /**< total number of variable type changes by this presolver */
   int              nchgbds;            /**< total number of variable bounds tightend by this presolver */
   int              naddholes;          /**< total number of domain holes added by this presolver */
   int              ndelconss;          /**< total number of deleted constraints by this presolver */
   int              nupgdconss;         /**< total number of upgraded constraints by this presolver */
   int              nchgcoefs;          /**< total number of changed coefficients by this presolver */
   int              nchgsides;          /**< total number of changed left or right hand sides by this presolver */
   Bool             delaysepa;          /**< should separation method be delayed, if other separators found cuts? */
   Bool             delayprop;          /**< should propagation method be delayed, if other propagators found reductions? */
   Bool             delaypresol;        /**< should presolving method be delayed, if other presolvers found reductions? */
   Bool             needscons;          /**< should the constraint handler be skipped, if no constraints are available? */
   Bool             sepawasdelayed;     /**< was the separation method delayed at the last call? */
   Bool             propwasdelayed;     /**< was the propagation method delayed at the last call? */
   Bool             presolwasdelayed;   /**< was the presolving method delayed at the last call? */
   Bool             initialized;        /**< is constraint handler initialized? */
   Bool             delayupdates;       /**< must the updates of the constraint arrays be delayed until processUpdates()? */
};


#endif
