/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal.c
 * @brief  main document page
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage SCIP (Solving Constraint Integer Programs)
 * @version  0.90k
 * @author   Tobias Achterberg
 * @author   Timo Berthold
 * @author   Thorsten Koch
 * @author   Alexander Martin
 * @author   Kati Wolter
 *
 * SCIP is a program and library to solve constraint integer programs (CIPs).
 *
 * - \ref CODE    "Coding style guidelines"
 * - \ref DOC     "How to search the documentation for interface methods"
 * - \ref CONS    "How to add constraint handlers"
 * - \ref PRICER  "How to add variable pricers"
 * - \ref PRESOL  "How to add presolvers"
 * - \ref SEPA    "How to add separators"
 * - \ref PROP    "How to add propagators"
 * - \ref BRANCH  "How to add branching rules"
 * - \ref NODESEL "How to add node selectors"
 * - \ref HEUR    "How to add heuristics"
 * - \ref RELAX   "How to add additional relaxations"
 * - \ref READER  "How to add file readers"
 * - \ref DIALOG  "How to add dialog options"
 * - \ref DISP    "How to add display columns"
 * - \ref OBJ     "Creating, capturing, releasing, and adding data objects"
 * - \ref PARAM   "Adding additional user parameters"
*/

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CODE Coding style guidelines
 *
 * Here we explain how your code should look like, so it is easy for us to include it into the SCIP distribution.
 *
 * - Indentation is 3 spaces. No tabs anywhere in the code.
 * - Always only one declaration in a line.
 * - Braces are on a new line and not indented.
 * - Spaces around all operators.
 * - Use assert() to show preconditions for the parameters, invariants and postconditions.
 * - All global functions start with "SCIP". In the usual naming scheme this is followed by the object and a method name
 *   like in SCIPlpAddRow(). Functions return TRUE or FALSE should be named like SCIPlpiIsOptimal().
 * - Make all functions that are not used outside the module 'static'. Naming should start with a lower case letter.
 * - Variable names should be all lower case.
 * - For each structure there is a typedef with the name in all upper case.
 * - Defines should be named all upper case.
 * - Document functions, parameters and variables doxygen conform.
 *
 * As an example have a look at tree.c .
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DOC How to search the documentation for interface methods
 *
 * If you are looking for a method in order to perform a specific task, there are usually two places to look at:
 * - The file "scip.h" in the file list.
 *   In this main header file, you find all methods that perform "complex" operations that affect or need data from
 *   different components of SCIP.
 *   For these methods, you always have to provide the SCIP pointer that is created by SCIPcreate().
 *   The documentation of "scip.h" is grouped into several blocks, each dealing with methods for a specific kind of
 *   object.
 *   For example, all methods operating on variables are grouped together.
 * - The files "pub_<...>.h" contains methods that perform "easy" operations that only affect the corresponding
 *   objects.
 *   Usually, with these methods you can access the data of the object.
 *   For example, in "pub_var.h" you find methods to get information about a variable.
 *
 * The file "pub_misc.h" contains methods for data structures like priority queues, hash tables, and hash maps,
 * as well as methods for sorting, numerics, random numbers, string operations, and file operations.
 *
 * If you are looking for a description of a callback method of a plugin that you want to implement, you have to
 * look at the corresponding "type_<...>.h".
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CONS How to add constraint handlers
 *
 * A constraint handler defines the semantics and the algorithms to process constraints of a certain class.
 * A single constraint handler is responsible for all the constraints belonging to his constraint class.
 * For example, there is one knapsack constraint handler that ensures that only solutions are accepted that
 * satisfy all the knapsack constraints in the model.
 *
 * In the following, we explain how the user can add an own constraint handler.
 * For example, look into the subtour constraint handler (examples/TSP/src/ConshdlrSubtour.cpp) of the
 * TSP example project.
 * The example is written in C++ and uses the C++ wrapper classes.
 * However, we will explain the implementation of a constraint handler using the C interface.
 * It is very easy to transfer the C explanation to C++: whenever a method should be implemented using the
 * SCIP_DECL_CONS... notion, reimplement the corresponding virtual member function of the abstract ObjConshdlr
 * base class.
 *
 * Additional documentation for the callback methods of a constraint handler can be found in the file
 * "type_cons.h".
 *
 * Here is what you have to do (assuming your constraint handler should be named "nosubtour"):
 * -# Copy the template files "src/scip/cons_xxx.c" and "src/scip/cons_xxx.h" into files "cons_nosubtour.c"
 *    and "cons_nosubtour.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "nosubtour".
 * -# Adjust the properties of the constraint handler (see \ref CONS_PROPERTIES).
 * -# Define the constraint data and the constraint handler data (see \ref CONS_DATA).
 * -# Implement the interface methods (see \ref CONS_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref CONS_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref CONS_ADDITIONALCALLBACKS).
 *
 * 
 * @section CONS_PROPERTIES Properties of a Constraint Handler
 *
 * At the top of the new file "cons_subtour.c" you can find the constraint handler properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the constraint handler properties by calling the constructor
 * of the abstract base class ObjConshdlr from within your constructor (see the TSP example).
 * The properties you have to set have the following meaning:
 *
 * \par CONSHDLR_NAME: the name of the constraint handler.
 * This name is used in the interactive shell to address the constraint handler.
 * Additionally, if you are searching a constraint handler with SCIPfindConshdlr(), this name is looked up.
 * Names have to be unique: no two constraint handlers may have the same name.
 *
 * \par CONSHDLR_DESC: the description of the constraint handler.
 * This string is printed as description of the constraint handler in the interactive shell.
 *
 * \par CONSHDLR_SEPAPRIORITY: the priority of the constraint handler for separation.
 * In each separation round during the price-and-cut loop of the subproblem processing or the separation loop
 * of the primal solution separation, the separators and separation methods of the constraint handlers are called in
 * a predefined order, which is given by the priorities of the separators and the separation priorities of the
 * constraint handlers.
 * First, the separators with non-negative priority are called in the order of decreasing priority.
 * Next, the separation methods of the different constraint handlers are called in the order of decreasing separation
 * priority.
 * Finally, the separators with negative priority are called in the order of decreasing priority.
 * \n
 * The separation priority of the constraint handler should be set according to the complexity of the cut separation
 * algorithm and the impact of the resulting cuts:
 * Constraint handlers that provide fast algorithms that usually have a high impact (i.e., cut off a large portion of
 * the LP relaxation) should have a high priority.
 * See \ref CONSSEPALP and \ref CONSSEPASOL for further details of the separation callbacks.
 *
 * \par CONSHDLR_ENFOPRIORITY: the priority of the constraint handler for constraint enforcing.
 * Like the separation priority, the enforcement priorities define the order in which the different constraint handlers
 * are called in the constraint enforcement step of the subproblem processing.
 * The constraint enforcement is called after the price-and-cut loop was executed (in the case that the LP is solved
 * at the current subproblem).
 * \n
 * The integrality constraint handler has an enforcement priority of 0.
 * That means, if a constraint handler has negative enforcement priority, it only has to deal with integral solutions
 * in its enforcement methods, because for fractional solutions, the integrality constraint handler would have
 * created a branching, thereby aborting the enforcement step.
 * If you want to implement a constraint-depending branching rule (for example, SOS branching on special ordered
 * set constraints), you have to assign a positive enforcement priority to your constraint handler.
 * In this case, you have to be able to deal with fractional solutions.
 * \n
 * See \ref CONSENFOLP and \ref CONSENFOPS for further details of the separation callback.
 *
 * \par CONSHDLR_CHECKPRIORITY: the priority of the constraint handler for checking feasibility.
 * Like the separation priority, the checking priorities define the order in which the different constraint handlers
 * are called to check the feasibility of a given primal solution candidate.
 * The integrality constraint handler has a checking priority of 0.
 * That means, constraint handlers with negative checking priorities only have to deal with integral solutions.
 * 
 * \par CONSHDLR_SEPAFREQ: the default frequency for separating cuts.
 * The separation frequency define the depth levels at which the constraint handler's separation methods \ref CONSSEPALP
 * and \ref CONSSEPASOL are called.
 * For example, a separation frequency of 7 means, that the separation callback is executed for subproblems that are
 * in depth 0, 7, 14, ... of the branching tree.
 * A separation frequency of 0 means, that the separation method is only called at the root node.
 * \n
 * The separation frequency can be adjusted by the user.
 * The property of the constraint handler only defines the default value of the frequency.
 * If you want to have a more flexible control of when to execute the separation algorithm, you have to assign
 * a separation frequency of 1 and implement a check at the beginning of your separation algorithm whether you really 
 * want to execute the separator or not.
 * If you do not want to execute the method, set the result code to SCIP_DIDNOTRUN.
 *
 * \par CONSHDLR_PROPFREQ: the default frequency for propagating domains.
 * This default frequency has the same meaning as the CONSHDLR_SEPAFREQ with respect to the domain propagation
 * callback of the consrtaint handler.
 * A propagation frequency of 0 means that propagation is only applied in preprocessing and at the root node.
 *
 * \par CONSHDLR_EAGERFREQ: the default frequency for using all instead of only the useful constraints in separation, propagation and enforcement.
 * If \em constraint \em aging is activated, some constraints that were not useful in the past for propagation or
 * separation are marked to be \em obsolete.
 * Usually, the obsolete constraints are not presented to the separation and propagation methods of the constraint
 * handlers, such that the constraint handlers only process the non-obsolete constraints.
 * However, every n'th call, with n being the EAGERFREQ of the constraint handler, all constraints are presented to the
 * separation and propagation methods of the constraint handler.
 * This gives obsolete constraints the chance of becoming non-obsolete again.
 * \n
 * If the eager evaluation frequency is set to -1, obsolete constraints are never presented to the separation and
 * propagation methods.
 * A frequency of 0 means, that obsolete constraints are only used in the first call of each method.
 *
 * \par CONSHDLR_MAXPREROUNDS: the default maximal number of presolving rounds the constraint handler participates in.
 * The preprocessing is executed in rounds.
 * If enough changes have been applied to the model, an additional preprocessing round is performed.
 * The MAXPREROUNDS parameter of a constraint handler denotes the maximal number of preprocessing rounds, the constraint
 * handler participates in.
 * A value of -1 means, that there is no limit on the number of rounds.
 * A value of 0 means, the preprocessing callback of the constraint handler is disabled.
 *
 * \par CONSHDLR_DELAYSEPA: the default for whether the separation method should be delayed, if other separators found cuts.
 * If the constraint handler's separation method is marked to be delayed, it is only executed after no other separator
 * or constraint handler found a cut during the price-and-cut loop.
 * If the separation method of the constraint handler is very expensive, you may want to mark it to be delayed after all
 * cheap separation methods have been executed.
 *
 * \par CONSHDLR_DELAYPROP: the default for whether the propagation method should be delayed, if other propagators found reductions.
 * This property is analoguos to the DELAYSEPA flag, but deals with the propagation method of the constraint handler.
 *
 * \par CONSHDLR_DELAYPRESOL: the default for whether the presolving method should be delayed, if other presolvers found reductions.
 * This property is analoguos to the DELAYSEPA flag, but deals with the preprocessing method of the constraint handler.
 *
 * \par CONSHDLR_NEEDSCONS: indicates whether the constraint handler be should skipped, if no constraints are available.
 * Usually, a constraint handler is only executed if there are constraints of its corresponding class in the model.
 * For those constraint handlers, the NEEDSCONS flag should be set to TRUE.
 * However, some constraint handlers must be called without having a constraint of the class in the model, because
 * the constraint is only implicitly available.
 * For example, the integrality constraint handler has the NEEDSCONS flag set to FALSE, because there is no explicit
 * integrality constraint in the model.
 * The integrality conditions are attached to the variables, and the integrality constraint handler has to check
 * all variables that are marked to be integer for integral values.
 *
 *
 * @section CONS_DATA Constraint Data and Constraint Handler Data
 *
 * Below the header "Data Structures" you can find two structs called "struct SCIP_ConsData" and
 * "struct SCIP_ConshdlrData".
 * If you are using C++, you only need to define the "struct SCIP_ConsData".
 * The constraint handler data must be implemented as member variables of your constraint handler class.
 * \n
 * The constraint data are the information that is needed to define a single constraint of the constraint handler's
 * constraint class.
 * For example, the data of a knapsack constraint would consist of a list of variables, a list of weights, and
 * the capacity of the knapsack.
 * The data of a nosubtour constraint consists of the graph on which the problem is defined.
 * In the graph, each edge should be linked to the corresponding binary problem variable.
 * \n
 * The constraint handler data are additional variables, that belong to the constraint handler itself and which are
 * not specific to a single constraint.
 * For example, you can use these data to store parameters of the constraint handler or statistical information.
 * The constraint handler data are optional.
 * You can leave the struct empty.
 *
 *
 * @section CONS_INTERFACE Interface Methods
 *
 * At the bottom of "cons_subtour.c" you can find two interface methods, that also appear in "cons_subtour.h".
 * These are SCIPincludeConshdlrSubtour() and SCIPcreateConsSubtour().
 * \n
 * The method \b SCIPincludeConshdlrSubtour() has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the constraint handler by calling the method
 * SCIPincludeConshdlr().
 * It is called by the user, if he wants to include the constraint handler, i.e. if he wants to make
 * the constraint handler available to the model.
 *
 * If you are using constraint handler data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_ConshdlrData afterwards.
 *
 * If the constraint handler is a specialization of a linear constraint, you may want to include an automatic
 * upgrade mechanism by calling the interface method
 * \code
 * SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdSubtour, LINCONSUPGD_PRIORITY) );
 * \endcode
 * of the linear constraint handler (see cons_linear.h).
 *
 * You may also add user parameters for your constraint handler.
 * An example for this can be found in src/scip/cons_knapsack.c.
 *
 * The method \b SCIPcreateConsSubtour() is called to create a single constraint of the constraint handler's constraint
 * class.
 * It should allocate and fill the constraint data, and call SCIPcreateCons().
 * Take a look at the following example from the logicor constraint handler:
 * \code
 * SCIP_RETCODE SCIPcreateConsLogicor(
 *    SCIP*                 scip,
 *    SCIP_CONS**           cons,
 *    const char*           name,
 *    int                   nvars,
 *    SCIP_VAR**            vars,
 *    SCIP_Bool             initial,
 *    SCIP_Bool             separate,
 *    SCIP_Bool             enforce,
 *    SCIP_Bool             check,
 *    SCIP_Bool             propagate,
 *    SCIP_Bool             local,
 *    SCIP_Bool             modifiable,
 *    SCIP_Bool             dynamic,
 *    SCIP_Bool             removeable
 *    )
 * {
 *    SCIP_CONSHDLR* conshdlr;
 *    SCIP_CONSDATA* consdata;
 * 
 *    assert(scip != NULL);
 * 
 *    conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
 *    if( conshdlr == NULL )
 *    {
 *       SCIPerrorMessage("logic or constraint handler not found\n");
 *       return SCIP_INVALIDCALL;
 *    }
 * 
 *    SCIP_CALL( consdataCreate(scip, &consdata, nvars, vars) );
 * 
 *    SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
 *          local, modifiable, dynamic, removeable) );
 * 
 *    return SCIP_OKAY;
 * }
 * \endcode
 * In this example, consdataCreate() is a local method that allocates memory for the given consdata and fills the
 * data with the given vars array.
 *
 * 
 * @section CONS_FUNDAMENTALCALLBACKS Fundamental Callback Methods
 *
 * By implementing the fundamental callbacks, you define the semantics of the constraint class the constraint handler
 * deals with.
 * If these methods are implemented, the resulting code is already correct and finds the optimal solution to the
 * given problem instance.
 * However, it might be very slow, because the additional features like cut separation and domain propagation are
 * missing.
 * In the C++ wrapper class ObjConshdlr, the fundamental callback methods are virtual abstract member functions.
 * You have to implement them in order to be able to construct an object of your constraint handler class.
 *
 * There are three fundamental callback methods that are all dealing with the feasibility of a given solution.
 * They are called at different places in the algorithm and have slightly different meaning.
 * However, it is usually reasonable to implement a single local method that is called by all of the three callback
 * methods with slightly modified parameters.
 * The fourth method provides dual information that is used for example in preprocessing.
 *
 * Additional documentation to the callback methods can be found in "type_cons.h".
 *
 * @subsection CONSCHECK
 *
 * The CONSCHECK callback gets a primal solution candidate in a SOL* data structure and has to check this solution
 * for global feasibility.
 * It has to return a result SCIP_FEASIBLE, if the solution satisfies all the constraints of the constraint handler,
 * and a result SCIP_INFEASIBLE if there is at least one constraint that is violated.
 * The callback is used by primal heuristics to check a constructed solution for feasibility.
 * That means, the constraint handler has to deal with arbitrary solutions that do not necessarily satisfy the bounds
 * and constraints of the local subproblem.
 *
 * The value of a variable \em var in the given solution \em sol can be accessed by calling
 * \code
 * SCIPgetSolVal(scip, sol, var)
 * \endcode
 *
 * For example, the knapsack constraint handler loops over his constraints and calculates the scalar product
 * \f$w^T x\f$ of weights \f$w\f$ with the solution vector \f$x\f$.
 * This scalar product is compared with the capacity of the knapsack constraint.
 * If it exceeds the capacity, the CONSCHECK method is immediately aborted with the result SCIP_INFEASIBLE.
 * If all knapsack constraints are satisfied, a result SCIP_FEASIBLE is returned.
 *
 * @subsection CONSENFOLP
 *
 * The CONSENFOLP method is called after the price-and-cut loop was finished and an LP solution is available.
 * Like the CHECK method, the ENFOLP method should check the solution (in this case, the LP solution) for
 * feasibility.
 * However, the solution is not given as a SOL* data structure.
 *
 * The value of a variable \em var in the LP solution can be accessed by calling
 * \code
 * SCIPgetVarSol(scip, var)
 * \endcode
 * or by
 * \code
 * SCIPgetSolVal(scip, NULL, var)
 * \endcode
 * By using the latter method, you can have a single local method to check a solution for feasibility by passing
 * the given \em sol in the CHECK call and by passing a NULL pointer as \em sol in the ENFOLP and ENFOPS calls.
 *
 * Like the CHECK call, the ENFOLP method should return a result SCIP_FEASIBLE, if the solution satisfies all the
 * constraints.
 * However, the behaviour should be different, if the solution violates one or more constraints.
 * The constraint handler may return a result SCIP_INFEASIBLE in this situation, but this is not the best what
 * one can do.
 * The ENFOLP method has the possibility of \em resolving the infeasibility by
 * - stating that the current subproblem is infeasible (result SCIP_CUTOFF),
 * - adding an additional constraint that resolves the infeasibility (result SCIP_CONSADDED),
 * - reducing the domain of a variable (result SCIP_REDUCEDDOM),
 * - adding a cutting plane (result SCIP_SEPARATED),
 * - performing a branching (result SCIP_BRANCHED).
 *
 * @subsection CONSENFOPS
 *
 * The CONSENFOPS callback is similar to the CONSENFOLP callback, but deals with \em pseudo \em solutions instead
 * of LP solutions.
 *
 * If the LP was not solved at the current subproblem (either because the user did not want to solve it, or because
 * numerical difficulties in the LP solving process were detected) no LP solution is available.
 * In this situation, the pseudo solution is used instead.
 * In this solution, the variables are set to the local bound which is best with respect to the objective function.
 * You can think of the pseudo solution as solution to the LP relaxation with all constraints except the bounds
 * being removed.
 *
 * Like the ENFOLP callback, the ENFOPS callback has to check whether the pseudo solution satisfies all the constraints
 * of the constraint handler.
 * The pseudo solution can be accessed by the same methods as the LP solution (SCIP knows, if the LP was solved at the
 * current subproblem, and returns either the LP solution or the pseudo solution).
 *
 * Unlike the ENFOLP callback, the ENFOPS callback must not add cuts and cannot return the result SCIP_SEPARATED.
 * It is, however, possible to force the solving of the LP by returning the result SCIP_SOLVELP.
 * For example, the infeasibility of a linear constraint that contains continuous variables cannot be resolved,
 * if all integer variables in the constraint are already fixed.
 * In this case, the LP has to be solved in order to get a solution that satisfies the linear constraint.
 *
 * @subsection CONSLOCK
 *
 * The CONSLOCK callback provides dual information for a single constraint.
 * It has to tell SCIP, which variables are existing in the given constraint, and in which way modifications of these
 * variables may affect the feasibility of the constraint.
 *
 * For each variable that is affected by the constraint, the callback should call SCIPaddVarLocks():
 *  - If the constraint may get violated by decreasing the value of a variable, it should call
 *    SCIPaddVarLocks(scip, var, nlockspos, nlocksneg), saying that rounding down is potentially rendering the
 *    (positive) constraint infeasible and rounding up is potentially rendering the negation of the constraint
 *    infeasible.
 *  - If the constraint may get violated by increasing the value of a variable, it should call
 *    SCIPaddVarLocks(scip, var, nlocksneg, nlockspos), saying that rounding up is potentially rendering the
 *    constraint's negation infeasible and rounding up is potentially rendering the constraint itself
 *    infeasible.
 *  - If the constraint may get violated by changing the variable in any direction, it should call
 *    SCIPaddVarLocks(scip, var, nlockspos + nlocksneg, nlockspos + nlocksneg).
 *
 *  Consider the linear constraint "3x -5y +2z <= 7" as an example. The CONSLOCK callback method of the
 *  linear constraint handler should call SCIPaddVarLocks(scip, x, nlocksneg, nlockspos), 
 *  SCIPaddVarLocks(scip, y, nlockspos, nlocksneg) and SCIPaddVarLocks(scip, z, nlocksneg, nlockspos) to tell SCIP,
 *  that rounding up of x and z and rounding down of y can destroy the feasibility of the constraint, while rounding
 *  down of x and z and rounding up of y can destroy the feasibility of the constraint's negation "3x -5y +2z > 7".
 *  A linear constraint "2 <= 3x -5y +2z <= 7" should call
 *  SCIPaddVarLocks(scip, ..., nlockspos + nlocksneg, nlockspos + nlocksneg) on all variables, since rounding in both
 *  directions of each variable can destroy both the feasibility of the constraint and it's negation
 *  "3x -5y +2z < 2  or  3x -5y +2z > 7".
 * 
 *
 * @section CONS_ADDITIONALCALLBACKS Additional Callback Methods
 *
 * The additional callback methods need not to be implemented in every case.
 * However, some of them have to be implemented for most applications.
 *
 * @subsection CONSFREE
 *
 * If you are using constraint handler data, you have to implement this method in order to free the constraint handler
 * data.
 * This can be done by the following procedure (which is taken from the knapsack constraint handler):
 * \code
 * static
 * SCIP_DECL_CONSFREE(consFreeKnapsack)
 * {
 *    SCIP_CONSHDLRDATA* conshdlrdata;
 *  
 *    conshdlrdata = SCIPconshdlrGetData(conshdlr);
 *    assert(conshdlrdata != NULL);
 *
 *    SCIPfreeMemory(scip, &conshdlrdata);
 *
 *    SCIPconshdlrSetData(conshdlr, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection CONSINIT
 *
 * The CONSINIT callback is executed after the problem was transformed.
 * The constraint handler may, e.g., use this call to replace the original variables in its constraints by transformed
 * variables, or to initialize his statistical constraint handler data.
 *
 * @subsection CONSEXIT
 *
 * The CONSEXIT callback is executed before the transformed problem is freed.
 * In this method, the constraint handler should free all resources that were allocated for the solving process.
 *
 * @subsection CONSINITPRE
 *
 * The CONSINITPRE callback is executed before the preprocessing is started, even if presolving is turned off.
 * The constraint handler may use this call to initialize its presolving data, or to modify its constraints
 * before the presolving process begins.
 * Necessary constraint modifications that have to be performed even if presolving is turned off should be done here
 * or in the presolving deinitialization call.
 *
 * @subsection CONSEXITPRE
 *
 * The CONSEXITPRE callback is executed after the preprocessing has been finished, even if presolving is turned off.
 * The constraint handler may use this call e.g. to clean up its presolving data, or to finally modify its constraints
 * before the branch and bound process begins.
 * Necessary constraint modifications that have to be performed even if presolving is turned off should be done here
 * or in the presolving initialization call.
 * Besides necessary modifications and clean up, no time consuming operations should be done.
 *
 * @subsection CONSINITSOL
 *
 * The CONSINITSOL callback is executed when the presolving was finished and the branch and bound process is about to
 * begin.
 * The constraint handler may use this call to initialize its branch and bound specific data.
 *
 * @subsection CONSEXITSOL
 *
 * The CONSEXITSOL callback is executed before the branch and bound process is freed.
 * The constraint handler should use this call to clean up its branch and bound data, in particular to release
 * all LP rows that he has created or captured.
 *
 * @subsection CONSDELETE
 *
 * The CONSDELETE callback is executed if a constraint should be freed.
 * You can think of it as the destructor of a single constraint.
 * In the callback, you have to free the given constraint data.
 * The CONSDELETE callback is therefore the counterpart of the SCIPcreateCons...() interface method and the CONSTRANS
 * method.
 *
 * @subsection CONSTRANS
 *
 * The CONSTRANS method is called for each constraint of the constraint handler, when the user starts the solving
 * process.
 * It has to copy the original constraint data of the constraint to the memory for the transformed problem.
 * You can think of it as a copy constructor for a single constraint.
 *
 * The original model is copied in order to protect it from transformations that are applied during the solving process,
 * in particular during preprocessing.
 * Preprocessing and solving always operates on the transformed problem.
 * If the solving process data are freed, the original data still exist and the user can, e.g., modify the problem and
 * restart the solving process.
 *
 * If you do not implement the CONSTRANS method, a transformed constraint is created with the same flags and the
 * same constraint data pointer.
 * That means, the transformed constraint points to the original constraint data.
 * This is okay, as long as the constraint data is not changed during the solving process.
 * If you want to implement preprocessing methods or other methods that modify the constraint data, you have to
 * implement the CONSTRANS method and create a copy of the constraint data.
 *
 * Here is an example, which is taken from the logicor constraint handler:
 * \code
 * static
 * SCIP_DECL_CONSTRANS(consTransLogicor)
 * {
 *    SCIP_CONSDATA* sourcedata;
 *    SCIP_CONSDATA* targetdata;
 * 
 *    sourcedata = SCIPconsGetData(sourcecons);
 * 
 *    SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->nvars, sourcedata->vars) );
 * 
 *    SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
 *          SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
 *          SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
 *          SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
 *          SCIPconsIsDynamic(sourcecons), SCIPconsIsRemoveable(sourcecons)) );
 * 
 *    return SCIP_OKAY;
 * }
 * \endcode
 * 
 * @subsection CONSINITLP
 *
 * The CONSINITLP callback is executed at the root node before the first LP relaxation is solved.
 * It should add the LP relaxations of all "initial" constraints to the LP. The method should scan the constraints
 * array for constraints that are marked initial via calls to SCIPconsIsInitial() and put the LP relaxation
 * of all initial constraints to the LP with calls to SCIPaddCut().
 * 
 * @subsection CONSSEPALP
 *
 * The CONSSEPALP callback is executed during the price-and-cut loop of the subproblem processing.
 * It should try to generate cutting planes for the constraints of the constraint handler in order to separate
 * the current LP solution.
 * The method is called in the LP solution loop, which means that a valid LP solution exists.
 *
 * Usually, a separation callback searches and produces cuts, that are added with a call to SCIPaddCut().
 * If the cut should be remembered in the global cut pool, it may also call SCIPaddPoolCut().
 * However, the callback may also produce domain reductions or add other constraints.
 *
 * The CONSSEPALP callback has the following options:
 *  - detecting that the node is infeasible in the variable's bounds and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint (result SCIP_CONSADDED)
 *  - reducing a variable's domain (result SCIP_REDUCEDDOM)
 *  - adding a cutting plane to the LP (result SCIP_SEPARATED)
 *  - stating that the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the separator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the separator was skipped, but should be called again (result SCIP_DELAYED)
 * 
 * @subsection CONSSEPASOL
 *
 * The CONSSEPASOL callback is executed during separation loop on arbitrary primal solutions.
 * It should try to generate cutting planes for the constraints of the constraint handler in order to separate
 * the given primal solution.
 * The method is not called in the LP solution loop, which means that there is no valid LP solution.
 *
 * Usually, a separation callback searches and produces cuts, that are added with a call to SCIPaddCut().
 * If the cut should be remembered in the global cut pool, it may also call SCIPaddPoolCut().
 * However, the callback may also produce domain reductions or add other constraints.
 *
 * The CONSSEPASOL callback has the following options:
 *  - detecting that the node is infeasible in the variable's bounds and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint (result SCIP_CONSADDED)
 *  - reducing a variable's domain (result SCIP_REDUCEDDOM)
 *  - adding a cutting plane to the LP (result SCIP_SEPARATED)
 *  - stating that the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the separator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the separator was skipped, but should be called again (result SCIP_DELAYED)
 * 
 * @subsection CONSPROP
 * 
 * The CONSPROP callback is called during the subproblem processing.
 * It should propagate the constraints, which means that it should infer reductions in the variables' local bounds
 * from the current local bounds.
 * This technique, which is the main workhorse of constraint programming, is called "node preprocessing" in the
 * Integer Programming community.
 *
 * The CONSPROP callback has the following options:
 *  - detecting that the node is infeasible in the variable's bounds and can be cut off (result SCIP_CUTOFF)
 *  - reducing a variable's domain (result SCIP_REDUCEDDOM)
 *  - stating that the propagator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the propagator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the propagator was skipped, but should be called again (result SCIP_DELAYED)
 *
 * @subsection CONSRESPROP
 *
 * If the constraint handler should support conflict analysis, it has to supply a CONSRESPROP method.
 * It also should call SCIPinferVarLbCons() or SCIPinferVarUbCons() in domain propagation instead of SCIPchgVarLb() or
 * SCIPchgVarUb() in order to deduce bound changes on variables.
 * In the SCIPinferVarLbCons() and SCIPinferVarUbCons() calls, the handler provides the constraint, that deduced the
 * variable's bound change, and an integer value "inferinfo" that can be arbitrarily chosen.
 *
 * The propagation conflict resolving method CONSRESPROP must then be implemented, to provide the "reasons" for the bound
 * changes, i.e. the bounds of variables at the time of the propagation, that forced the constraint to set the
 * conflict variable's bound to its current value. It can use the "inferinfo" tag to identify its own propagation
 * rule and thus identify the "reason" bounds. The bounds that form the reason of the assignment must then be provided
 * by calls to SCIPaddConflictLb() and SCIPaddConflictUb() in the propagation conflict resolving method.
 *
 * For example, the logicor constraint c = "x or y or z" fixes variable z to TRUE (i.e. changes the lower bound of z
 * to 1.0), if both, x and y, are assigned to FALSE (i.e. if the upper bounds of these variables are 0.0). It uses
 * SCIPinferVarLbCons(scip, z, 1.0, c, 0) to apply this assignment (an inference information tag is not needed by the
 * constraint handler and is set to 0).
 * In the conflict analysis, the constraint handler may be asked to resolve the lower bound change on z with
 * constraint c, that was applied at a time given by a bound change index "bdchgidx".
 * With a call to SCIPvarGetLbAtIndex(z, bdchgidx), the handler can find out, that the lower bound of variable z was
 * set to 1.0 at the given point of time, and should call SCIPaddConflictUb(scip, x, bdchgidx) and
 * SCIPaddConflictUb(scip, y, bdchgidx) to tell SCIP, that the upper bounds of x and y at this point of time were
 * the reason for the deduction of the lower bound of z.
 *
 * @subsection CONSPRESOL
 *
 * The CONSPRESOL callback is called during preprocessing.
 * It should try to tighten the domains of the variables, tighten the coefficients of the constraints of the constraint
 * handler, delete redundant constraints, aggregate and fix variables if possible, and upgrade constraints to a more
 * specific type.
 *
 * If the CONSPRESOL callback applies changes to the constraint data, you also have to implement the CONSTRANS callback
 * in order to copy the constraint data to the transformed problem space and protect the original problem from the
 * preprocessing changes.
 *
 * @subsection CONSACTIVE
 *
 * The CONSACTIVE callback method is called each time, a constraint of the constraint handler is activated.
 * For example, if a constraint is added locally to a subproblem, the CONSACTIVE callback is called whenever the
 * search enters the subtree where the constraint exists.
 *
 * @subsection CONSDEACTIVE
 *
 * The CONSDEACTIVE callback method is called each time, a constraint of the constraint handler is deactivated.
 * For example, if a constraint is added locally to a subproblem, the CONSDEACTIVE callback is called whenever the
 * search leaves the subtree where the constraint exists.
 *
 * @subsection CONSENABLE
 *
 * The CONSENABLE callback method is called each time, a constraint of the constraint handler is enabled.
 * Constraints might be active without being enabled. In this case, only the feasibility checks are executed,
 * but domain propagation and separation is skipped.
 *
 * @subsection CONSDISABLE
 *
 * The CONSDISABLE callback method is called each time, a constraint of the constraint handler is disabled.
 *
 * @subsection CONSPRINT
 *
 * The CONSPRINT callback method is called, when the user asks SCIP to display the problem to the screen or same
 * the problem into a file.
 * The constraint handler should display the data of the constraint in an appropriate form.
 * The output format that is defined by the CONSPRINT callbacks is called CIP format.
 * In later versions of SCIP, the constraint handlers should also be able to parse (i.e., read) constraints
 * which are given in CIP format.
 *
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page PRICER How to add variable pricers
 *
 * A pricer performs the dynamic generation of new variables in a column generation algorithm.
 * It is an algorithmic representation of a (usually exponential) number of variables.
 * The \ref PRICERREDCOST and \ref PRICERFARKAS methods are called after each LP solve to generate additional
 * variables which may improve the objective value or decrease the LP infeasibility, respectively.
 *
 * If the pricer finds one or more variables with negative reduced costs or negative farkas value, it should
 * call SCIPcreateVar() and SCIPaddPricedVar() to create and add the variable to the problem. Additionally,
 * the pricer has to add the variable to all constraints in which it appears. Therefore, a pricer needs to
 * know the constraints of the model and their meaning. Note that all constraints for which additional variables
 * are generated by a pricer have to be flagged as "modifiable" in the SCIPcreateCons() call.
 *
 * In the following, we explain how the user can add an own pricer.
 * For example, look into the distance pricer for the p-median problem (examples/SamplePricer/src/pricer_distance.h) of the
 * SamplePricer example project.
 * The example is written in C++ and uses the C++ wrapper classes.
 * However, we will explain the implementation of a pricer using the C interface.
 * It is very easy to transfer the C explanation to C++: whenever a method should be implemented using the
 * SCIP_DECL_PRICER... notion, reimplement the corresponding virtual member function of the abstract ObjPricer
 * base class.
 *
 * Additional documentation for the callback methods of a pricer can be found in the file
 * "type_pricer.h".
 *
 * Here is what you have to do to implement a pricer:
 * -# Copy the template files "src/scip/pricer_xxx.c" and "src/scip/pricer_xxx.h" into files "pricer_mypricer.c"
 *    and "pricer_mypricer.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "mypricer".
 * -# Adjust the properties of the pricer (see \ref PRICER_PROPERTIES).
 * -# Define the pricer data (see \ref PRICER_DATA).
 * -# Implement the interface methods (see \ref PRICER_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref PRICER_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref PRICER_ADDITIONALCALLBACKS).
 *
 * 
 * @section PRICER_PROPERTIES Properties of a Pricer
 *
 * At the top of the new file "pricer_mypricer.c" you can find the pricer properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the pricer properties by calling the constructor
 * of the abstract base class ObjPricer from within your constructor (see the SamplePricer example).
 * The properties you have to set have the following meaning:
 *
 * \par PRICER_NAME: the name of the pricer.
 * This name is used in the interactive shell to address the pricer.
 * Additionally, if you are searching a pricer with SCIPfindPricer(), this name is looked up.
 * Names have to be unique: no two pricers may have the same name.
 *
 * \par PRICER_DESC: the description of the pricer.
 * This string is printed as description of the pricer in the interactive shell.
 *
 * \par PRICER_PRIORITY: the priority of the pricer.
 * In each pricing round during the price-and-cut loop of the subproblem processing, the included pricers are
 * called in a predefined order, which is given by the priorities of the pricers.
 * The higher the priority, the earlier the pricer is called.
 * Usually, you will have only one pricer in your application and the priority is therefore irrelevant.
 *
 * \par PRICER_DELAY: the default for whether the pricer should be delayed, if other variables with negative reduced costs have already been found in the current pricing round.
 * Variables may be declared to be "removable" in the SCIPcreateVar() call. This means that SCIP may remove the variable
 * from the LP if it was inactive (i.e., sitting at zero) for a number of LP solves. Nevertheless, after the removal of the
 * column from the LP, the variable still exists, and SCIP can calculate reduced costs and add it to the LP again if
 * necessary.
 * \n
 * If the PRICER_DELAY flag is set to TRUE (which is the common setting), all those existing variables with negative reduced costs
 * are added to the LP, and the LP is resolved before the pricer is called. Thus, the pricer can assume that all existing variables
 * have non-negative reduced costs if the \ref PRICERREDCOST or \ref PRICERFARKAS methods are called.
 * \n
 * In some applications, this inner pricing loop on the already existing variables can significantly slow down the solving process,
 * since it may lead to the addition of only very few variables in each pricing round. If this is an issue in your application,
 * you should consider to set the PRICER_DELAY flag to FALSE. You must, however, be aware of the fact that there may be already
 * existing variables with negative reduced costs. For example, this may lead to the issue that your pricer generates the same
 * variable twice. In some models, this is not critical because an optimal solution would choose only one of the two identical
 * variables anyway, but for other models this can lead to wrong results because the duplication of a variable essentially doubles
 * the upper bound of the variable.
 *
 *
 * @section PRICER_DATA Pricer Data
 *
 * Below the header "Data Structures" you can find a struct which is called "struct SCIP_PricerData".
 * In this data structure, you can store the data of your pricer. For example, it may be convenient to store pointers to the
 * constraints of the problem instance here, because the pricer has to add variables to those constraints.
 * If you are using C++, you can add pricer data as usual as object variables to your class.
 * \n
 * Defining pricer data is optional. You can leave the struct empty.
 *
 *
 * @section PRICER_INTERFACE Interface Methods
 *
 * At the bottom of "pricer_mypricer.c" you can find the interface method SCIPincludePricerMypricer(), which also appears in "pricer_mypricer.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the pricer by calling the method
 * SCIPincludePricer().
 * It is called by the user, if he wants to include the pricer, i.e. if he wants to solve a model for which variables should
 * be generated by this pricer.
 *
 * If you are using pricer data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &pricerdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_PricerData afterwards.
 *
 * You may also add user parameters for your pricer, see the method \b SCIPincludeConshdlrKnapsack() in src/scip/cons_knapsack.c
 * for an example.
 *
 * 
 * @section PRICER_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Pricer
 *
 * The fundamental callback methods have to be implemented in order to obtain an operational algorithm.
 * In the case of a pricer, the fundamental callbacks are the two variable pricing callbacks which search and add
 * new variables to the problem.
 * In the C++ wrapper class ObjConshdlr, the scip_redcost() method (which corresponds to the PRICERREDCOST callback) is a virtual
 * abstract member function.
 * You have to implement it in order to be able to construct an object of your pricer class.
 * The other fundamental method scip_farkas() has an empty default implementation, but if it may happen in your application
 * that infeasible LP relaxations are encountered, you have to overwrite this method as well.
 *
 * Additional documentation to the callback methods can be found in "type_pricer.h".
 *
 * @subsection PRICERREDCOST
 *
 * The PRICERREDCOST callback is called inside the price-and-cut loop of the subproblem solving process if the current LP relaxation
 * is feasible.
 * It should search for additional variables that can contribute to improve the current LP's solution value.
 * In standard branch-and-price, these are variables with negative dual feasibility, that is negative
 * reduced costs for non-negative variables, positive reduced costs for non-positive variables,
 * and non-zero reduced costs for variables that can be negative and positive.
 *
 * Whenever the pricer finds a variable with negative dual feasibility, it should call SCIPcreateVar()
 * and SCIPaddPricedVar() to add the variable to the problem. Furthermore, it should call the appropriate
 * methods of the constraint handlers to add the necessary variable entries to the constraints.
 *
 * Pricers usually need the dual LP solution as input for the pricing algorithm.
 * Since SCIP does not know the semantics of the individual constranits in the problem, the dual solution
 * has to be provided by the constraint handlers.
 * For example, the setppc constraint handler that deals with set partitioning, packing, and covering constraints provides
 * the method SCIPgetDualsolSetppc() (defined in cons_setppc.h) to access the dual solution value for a single constraint.
 * Similarly, the dual solution of a linear constraint can be queried with the method SCIPgetDualsolLinear() of cons_linear.h.
 * The reduced costs of the existing variables can be accessed with the method SCIPgetVarRedcost().
 *
 * @subsection PRICERFARKAS
 *
 * If the current LP relaxation is infeasible, it is the task of the pricer to generate additional variables that can
 * potentially render the LP feasible again. In standard branch-and-price, these are variables with positive farkas values,
 * and the PRICERFARKAS method should identify those variables.
 *
 * If the LP was proven to be infeasible, we have an infeasibility proof by the dual farkas multipliers \f$y\f$.
 * With the values of y, an implicit inequality \f$y^T A x \ge y^T b\f$ is associated, with \f$b\f$ given
 * by the sides of the LP rows and the sign of \f$y\f$:
 *  - if \f$y_i\f$ is positive, \f$b_i\f$ is the left hand side of the row,
 *  - if \f$y_i\f$ is negative, \f$b_i\f$ is the right hand side of the row.
 *
 * \f$y\f$ is chosen in a way, such that the valid inequality  \f$y^T A x \ge y^T b\f$  is violated by all \f$x\f$,
 * especially by the (for this inequality least infeasible solution) \f$x'\f$ defined by 
 *  - \f$x'_i := ub_i\f$, if \f$y^T A_i \ge 0\f$
 *  - \f$x'_i := lb_i\f$, if \f$y^T A_i < 0\f$.
 * Pricing in this case means to add variables \f$i\f$ with positive farkas value, i.e. \f$y^T A_i x'_i > 0\f$.
 *
 * To apply farkas pricing, the pricer needs to know the farkas values of the constraints. Like the dual solution values for
 * feasible LP solutions, the dual farkas values for infeasible solutions can be obtained by constraint handler interface
 * methods like, for example, the SCIPgetDualfarkasLinear() method of the linear constraint handler.
 * The farkas values for the bounds of the variables are just the regular reduced costs and can be accessed with SCIPgetVarRedcost().
 *
 * It is useful to note that farkas pricing is the same as the regular pricing with a zero objective function.
 * Therefore, a typical implementation of a pricer would consist of a generic pricing algorithm that gets a dual solution and an
 * objective function vector as input and generates variables by calling SCIPcreateVar() and SCIPaddPricedVar().
 * The PRICERREDCOST callback would call this function with the regular objective function and the regular dual solution vector,
 * while the PRICERFARKAS callback would call this function with a zero objective function and the farkas vector.
 * From a practical point of view, it is usually the simplest approach to provide just one Boolean flag to the generic pricing
 * algorithm in order to identify whether it is reduced cost or farkas pricing. Then, the algorithm would just call the appropriate
 * methods to access the dual solution or objective function, depending on the Boolean flag.
 *
 *
 * @section PRICER_ADDITIONALCALLBACKS Additional Callback Methods of a Pricer
 *
 * The additional callback methods need not to be implemented in every case.
 * However, some of them have to be implemented for most applications.
 *
 * @subsection PRICERFREE
 *
 * If you are using pricer data, you have to implement this method in order to free the pricer data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_PRICERFREE(pricerFreeMypricer)
 * {
 *    SCIP_PRICERDATA* pricerdata;
 *  
 *    pricerdata = SCIPpricerGetData(pricer);
 *    assert(pricerdata != NULL);
 *
 *    SCIPfreeMemory(scip, &pricerdata);
 *
 *    SCIPpricerSetData(pricer, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection PRICERINIT
 *
 * The PRICERINIT callback is executed after the problem was transformed.
 * The pricer may, e.g., use this call to replace the original constraints stored in its pricer data by transformed
 * constraints, or to initialize other elements of his pricer data.
 *
 * @subsection PRICEREXIT
 *
 * The PRICEREXIT callback is executed before the transformed problem is freed.
 * In this method, the pricer should free all resources that have been allocated for the solving process in PRICERINIT.
 *
 * @subsection PRICERINITSOL
 *
 * The PRICERINITSOL callback is executed when the presolving was finished and the branch and bound process is about to begin.
 * The pricer may use this call to initialize its branch and bound specific data.
 *
 * @subsection PRICEREXITSOL
 *
 * The PRICEREXITSOL callback is executed before the branch and bound process is freed.
 * The pricer should use this call to clean up its branch and bound data, which was allocated in PRICERINITSOL.
 *
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page PRESOL How to add presolvers
 *
 * Presolvers are used to simplify the problem instance and to extract useful information in the presolving step.
 * Constraint based presolving is done in the CONSPRESOL callback methods of the constraint handlers, see \ref CONSPRESOL.
 * The presolver plugins complement the constraint based presolving by additional, usually optimality based, presolving
 * reductions.
 *
 * In the following, we explain how the user can add an own presolver.
 * Take the dual fixing presolver (src/scip/presol_dualfix.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the ObjPresol wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_PRESOL_... callback methods.
 *
 * Additional documentation for the callback methods of a presolver can be found in the file "type_presol.h".
 *
 * Here is what you have to do to implement a presolver:
 * -# Copy the template files "src/scip/presol_xxx.c" and "src/scip/presol_xxx.h" into files names "presol_mypresolver.c"
 *    and "presol_mypresolver.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "mypresolver".
 * -# Adjust the properties of the presolver (see \ref PRESOL_PROPERTIES).
 * -# Define the presolver data (see \ref PRESOL_DATA).
 * -# Implement the interface methods (see \ref PRESOL_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref PRESOL_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref PRESOL_ADDITIONALCALLBACKS).
 *
 * 
 * @section PRESOL_PROPERTIES Properties of a Presolver
 *
 * At the top of the new file "presol_mypresolver.c" you can find the presolver properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the presolver properties by calling the constructor
 * of the abstract base class ObjPresol from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par PRESOL_NAME: the name of the presolver.
 * This name is used in the interactive shell to address the presolver.
 * Additionally, if you are searching a presolver with SCIPfindPresol(), this name is looked up.
 * Names have to be unique: no two presolvers may have the same name.
 *
 * \par PRESOL_DESC: the description of the presolver.
 * This string is printed as description of the presolver in the interactive shell.
 *
 * \par PRESOL_PRIORITY: the priority of the presolver.
 * In each presolving round, the presolvers and presolving methods of the constraint handlers are called in
 * a predefined order, which is given by the priorities of the presolvers and the check priorities of the
 * constraint handlers.
 * First, the presolvers with non-negative priority are called in the order of decreasing priority.
 * Next, the presolving methods of the different constraint handlers are called in the order of decreasing check
 * priority.
 * Finally, the presolvers with negative priority are called in the order of decreasing priority.
 * \n
 * The priority of the presolver should be set according to the complexity of the presolving algorithm and the impact of the reduction:
 * presolvers that provide fast algorithms that usually have a high impact (i.e., remove lots of variables or tighten many bounds)
 * should have a high priority.
 *
 * \par PRESOL_MAXROUNDS: the default maximal number of rounds the presolver participates in.
 * The presolving is conducted in rounds: the presolvers and presolving methods of the constraint handlers
 * are called iteratively until no more reductions have been found or some other abort criterion applies.
 * The "maxrounds" parameter of a presolver imposes a limit on the number of presolving rounds in which the
 * presolver is called. The PRESOL_MAXROUNDS property specifies the default value for this parameter.
 * A value of -1 represents an unlimited number of rounds.
 *
 * \par PRESOL_DELAYPRESOL: the default for whether the presolver should be delayed, if other presolvers found reductions.
 * If the presolver is marked to be delayed, it is only executed if no other presolvers found a reduction during the current
 * presolving round.
 * If the presolver is very expensive, you may want to mark it to be delayed after all cheap presolving methods have been executed.
 * 
 *
 * @section PRESOL_DATA Presolver Data
 *
 * Below the header "Data Structures" you can find a struct which is called "struct SCIP_PresolData".
 * In this data structure, you can store the data of your presolver. For example, you should store the adjustable parameters
 * of the presolver in this data structure.
 * If you are using C++, you can add presolver data as usual as object variables to your class.
 * \n
 * Defining presolver data is optional. You can leave the struct empty.
 *
 *
 * @section PRESOL_INTERFACE Interface Methods
 *
 * At the bottom of "presol_mypresolver.c" you can find the interface method SCIPincludePresolMypresolver(), which also appears in "presol_mypresolver.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the presolver by calling the method
 * SCIPincludePresol().
 * It is called by the user, if he wants to include the presolver, i.e. if he wants to use your presolver in his application.
 *
 * If you are using presolver data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &presoldata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_PresolData afterwards.
 *
 * You may also add user parameters for your presolver, see the method \b SCIPincludeConshdlrKnapsack() in src/scip/cons_knapsack.c
 * for an example.
 *
 * 
 * @section PRESOL_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Presolver
 *
 * Presolver plugins have only one fundamental callback method, namely the PRESOLEXEC method.
 * This method has to be implemented for every presolver; the other callback methods are optional.
 * In the C++ wrapper class ObjConshdlr, the scip_exec() method (which corresponds to the PRESOLEXEC callback) is a virtual
 * abstract member function.
 * You have to implement it in order to be able to construct an object of your presolver class.
 *
 * Additional documentation to the callback methods can be found in "type_presol.h".
 *
 * @subsection PRESOLEXEC
 *
 * The PRESOLEXEC callback is called inside the presolving loop and should perform the actual presolving reductions.
 * It should inspect the problem instance at hand and simplify it by tightening bounds of variables, aggregating or fixing
 * variables, changing the type of variables, modifying the graph that represents the instance of your application, and
 * the like.
 *
 * Typical methods called by a presolver are, for example, SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), SCIPtightenVarLb(),
 * and SCIPtightenVarUb().
 *
 *
 * @section PRESOL_ADDITIONALCALLBACKS Additional Callback Methods of a Presolver
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 *
 * @subsection PRESOLFREE
 *
 * If you are using presolver data, you have to implement this method in order to free the presolver data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_PRESOLFREE(presolFreeMypresolver)
 * {
 *    SCIP_PRESOLDATA* presoldata;
 *  
 *    presoldata = SCIPpresolGetData(presol);
 *    assert(presoldata != NULL);
 *
 *    SCIPfreeMemory(scip, &presoldata);
 *
 *    SCIPpresolSetData(presol, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection PRESOLINIT
 *
 * The PRESOLINIT callback is executed after the problem was transformed.
 * The presolver may, e.g., use this call to initialize his presolver data.
 *
 * @subsection PRESOLEXIT
 *
 * The PRESOLEXIT callback is executed before the transformed problem is freed.
 * In this method, the presolver should free all resources that have been allocated for the solving process in PRESOLINIT.
 *
 * @subsection PRESOLINITPRE
 *
 * The PRESOLINITPRE callback is executed when the presolving is about to begin.
 * The presolver may use this call to initialize its presolving data which only need to exist during the presolving stage.
 *
 * @subsection PRESOLEXITPRE
 *
 * The PRESOLEXITPRE callback is executed after presolving has been finished and before the branch and bound process begins.
 * The presolver should use this call to clean up its presolving data, which was allocated in PRESOLINITPRE.
 *
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page SEPA How to add separators
 *
 * Separators are used to generate general purpose cutting planes. 
 * Constraint based cutting planes, the second type of cutting planes in SCIP, are separated in the CONSSEPALP and 
 * CONSSEPASOL callback methods of the constraint handlers, see \ref CONSSEPALP and \ref CONSSEPASOL. These cuts are 
 * valid inequalities or even facets of the polyhedron described by a single constraint or a subset of the constraints of
 * a single constraint class. In contrast, general purpose cuts do not require or exploit any knowledge about the 
 * underlying problem structure but use only the current LP relaxation and the integrality conditions. 
 *
 * In the following, we explain how the user can add an own separator.
 * Take the separator for the class of Gomory mixed inter inequalities (src/scip/sepa_gomory.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the ObjSepa wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_SEPA... callback methods.
 *
 * Additional documentation for the callback methods of a separator can be found in the file "type_sepa.h".
 *
 * Here is what you have to do to implement a separator:
 * -# Copy the template files "src/scip/sepa_xxx.c" and "src/scip/sepa_xxx.h" into files "sepa_myseparator.c"
 *    and "sepa_myseparator.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "myseparator".
 * -# Adjust the properties of the separator (see \ref SEPA_PROPERTIES).
 * -# Define the separator data (see \ref SEPA_DATA).
 * -# Implement the interface methods (see \ref SEPA_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref SEPA_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref SEPA_ADDITIONALCALLBACKS).
 *
 *
 * @section SEPA_PROPERTIES Properties of a Separator
 *
 * At the top of the new file "sepa_myseparator.c" you can find the separator properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the separator properties by calling the constructor
 * of the abstract base class ObjSepa from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par SEPA_NAME: the name of the separator.
 * This name is used in the interactive shell to address the separator.
 * Additionally, if you are searching a separator with SCIPfindSepa(), this name is looked up.
 * Names have to be unique: no two separators may have the same name.
 *
 * \par SEPA_DESC: the description of the separator.
 * This string is printed as description of the separator in the interactive shell.
 *
 * \par SEPA_PRIORITY: the priority of the separator.
 * In each separation round during the price-and-cut loop of the subproblem processing or the separation loop
 * of the primal solution separation, the separators and separation methods of the constraint handlers are called in
 * a predefined order, which is given by the priorities of the separators and the separation priorities of the
 * constraint handlers.
 * First, the separators with non-negative priority are called in the order of decreasing priority.
 * Next, the separation methods of the different constraint handlers are called in the order of decreasing separation
 * priority.
 * Finally, the separators with negative priority are called in the order of decreasing priority.
 * \n
 * The priority of the separator should be set according to the complexity of the cut separation algorithm and the 
 * impact of the resulting cuts: separators that provide fast algorithms that usually have a high impact (i.e., cut off 
 * a large portion of the LP relaxation) should have a high priority.
 * See \ref SEPAEXECLP and \ref SEPAEXECSOL for further details of the separation callbacks.
 *
 * \par SEPA_FREQ: the default frequency for separating cuts.
 * The frequency defines the depth levels at which the separation methods \ref SEPAEXECLP and \ref SEPAEXECSOL are called.
 * For example, a frequency of 7 means, that the separation callback is executed for subproblems that are in depth 
 * 0, 7, 14, ... of the branching tree. A frequency of 0 means, that the separation method is only called at the root node.
 * \n
 * The frequency can be adjusted by the user. The property of the separator only defines the default value of the frequency.
 * If you want to have a more flexible control of when to execute the separation algorithm, you have to assign
 * a frequency of 1 and implement a check at the beginning of your separation algorithm whether you really 
 * want to execute the separation or not. If you do not want to execute it, set the result code to SCIP_DIDNOTRUN.
 *
 * \par SEPA_MAXBOUNDDIST: the default maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for applying separation.
 * At the current branch-and-bound node, the relative distance from its dual bound (local dual bound) 
 * to the primal bound compared to the best node's dual bound (global dual bound) is considered. The separation method 
 * of the separator will only be applied at the node if this relative distance does not exceed SEPA_MAXBOUNDDIST. 
 * \n
 * For example, if the global dual bound is 0 and the primal bound is 10, SEPA_MAXBOUNDDIST = 0.25 means that separation 
 * is only applied if the current node's dual bound is in the first quater of the interval [0,10], i.e., if it is less 
 * than or equal to 2.5. 
 * \n
 * In particular, the extremal values 0.0 and 1.0 mean that separtion is only applied at the current best node and at all 
 * nodes, respectively. Since separation seems to be most effective when applied at nodes that contribute to the global 
 * dual bound, 0.0 is probably a good choice for SEPA_MAXBOUNDDIST.
 *
 * \par SEPA_DELAY: the default for whether the separation method should be delayed, if other separators or constraint handlers found cuts.
 * If the separator's separation method is marked to be delayed, it is only executed after no other separator
 * or constraint handler found a cut during the price-and-cut loop. 
 * If the separation method of the separator is very expensive, you may want to mark it to be delayed after all cheap 
 * separation methods have been executed.
 *
 * @section SEPA_DATA Separator Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_SepaData".
 * In this data structure, you can store the data of your separator. For example, you should store the adjustable parameters
 * of the separator in this data structure.
 * If you are using C++, you can add separator data as usual as object variables to your class.
 * \n
 * Defining separator data is optional. You can leave the struct empty.
 *
 * @section SEPA_INTERFACE Interface Methods
 *
 * At the bottom of "sepa_myseparator.c" you can find the interface method SCIPincludeSepaMyseparator(), which also 
 * appears in "sepa_myseparator.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the separator by calling the method SCIPincludeSepa().
 * It is called by the user, if he wants to include the separator, i.e. if he wants to use the separator in his application.
 *
 * If you are using separator data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &sepadata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_SepaData afterwards.
 *
 * You may also add user parameters for your separator, see the method \b SCIPincludeConshdlrKnapsack() 
 * in src/scip/cons_knapsack.c for an example.
 *
 * 
 * @section SEPA_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Separator
 *
 * Separator plugins do not have any fundamental callback methods. All callback methods are optional.
 *
 *
 * @section SEPA_ADDITIONALCALLBACKS Additional Callback Methods of a Separator
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 *
 * @subsection SEPAFREE
 *
 * If you are using separator data, you have to implement this method in order to free the separator data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_SEPAFREE(sepaFreeMyseparator)
 * {
 *    SCIP_SEPADATA* sepadata;
 *  
 *    sepadata = SCIPsepaGetData(sepa);
 *    assert(sepadata != NULL);
 *
 *    SCIPfreeMemory(scip, &sepadata);
 *
 *    SCIPsepaSetData(sepa, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection SEPAINIT
 *
 * The SEPAINIT callback is executed after the problem was transformed.
 * The separator may, e.g., use this call to initialize his separator data.
 *
 * @subsection SEPAEXIT
 *
 * The SEPAEXIT callback is executed before the transformed problem is freed.
 * In this method, the separator should free all resources that have been allocated for the solving process in SEPAINIT.
 *
 * @subsection SEPAINITSOL
 *
 * The SEPAINITSOL callback is executed when the presolving was finished and the branch and bound process is about to 
 * begin. The separator may use this call to initialize its branch and bound specific data.
 *
 * @subsection SEPAEXITSOL
 *
 * The SEPAEXITSOL callback is executed before the branch and bound process is freed. The separator should use this call 
 * to clean up its branch and bound data, in particular to release all LP rows that he has created or captured.
 *
 * @subsection SEPAEXECLP
 *
 * The SEPAEXECLP callback is executed during the price-and-cut loop of the subproblem processing.
 * It should try to generate general purpose cutting planes in order to separate the current LP solution.
 * The method is called in the LP solution loop, which means that a valid LP solution exists.
 *
 * Usually, the callback searches and produces cuts, that are added with a call to SCIPaddCut().
 * If the cut should be remembered in the global cut pool, it may also call SCIPaddPoolCut().
 * However, it may also produce domain reductions or add other constraints.
 *
 * The SEPAEXECLP callback has the following options:
 *  - detecting that the node is infeasible in the variable's bounds and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint (result SCIP_CONSADDED)
 *  - reducing a variable's domain (result SCIP_REDUCEDDOM)
 *  - adding a cutting plane to the LP (result SCIP_SEPARATED)
 *  - stating that the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the separator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the separator was skipped, but should be called again (result SCIP_DELAYED)
 *
 * @subsection SEPAEXECSOL
 *
 * The SEPAEXECSOL callback is executed during separation loop on arbitrary primal solutions.
 * It should try to generate general purpose cutting planes in order to separate the given primal solution.
 * The method is not called in the LP solution loop, which means that there is no valid LP solution.
 *
 * Usually, the callback searches and produces cuts, that are added with a call to SCIPaddCut().
 * If the cut should be remembered in the global cut pool, it may also call SCIPaddPoolCut().
 * However, it may also produce domain reductions or add other constraints.
 *
 * The SEPAEXECSOL callback has the following options:
 *  - detecting that the node is infeasible in the variable's bounds and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint (result SCIP_CONSADDED)
 *  - reducing a variable's domain (result SCIP_REDUCEDDOM)
 *  - adding a cutting plane to the LP (result SCIP_SEPARATED)
 *  - stating that the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the separator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the separator was skipped, but should be called again (result SCIP_DELAYED)
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page PROP How to add propagators
 *
 * This page is not yet written. Here we will explain how to add propagation routines to SCIP.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page BRANCH How to add branching rules
 *
 * This page is not yet written. Here we will explain how to add branching rules to SCIP.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page NODESEL How to add node selectors
 *
 * This page is not yet written. Here we will explain how to add node selection rules to SCIP.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page HEUR How to add heuristics
 *
 * This page is not yet written. Here we will explain how to add primal heuristics to SCIP.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page RELAX How to add additional relaxations
 *
 * This page is not yet written. Here we will explain how to add additional relaxations (apart from the LP relaxation)
 * to SCIP.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page READER How to add file readers
 *
 * This page is not yet written. Here we will explain how to add input file readers to SCIP.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DIALOG How to add dialog options
 *
 * This page is not yet written. Here we will explain how to extend the interactive shell by adding new dialog options
 * to SCIP.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DISP How to add display columns
 *
 * This page is not yet written. Here we will explain how to add additional display columns to SCIP.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page OBJ Creating, capturing, releasing, and adding data objects
 *
 *  Data objects (variables, constraints, rows) are subject to reference counting
 *  to avoid expensive copying operations. Creating such an object will set the
 *  reference count to one. Capturing an object increases the reference counter,
 *  releasing it decreases the counter. If the reference counter gets zero, the
 *  object is destroyed.
 *
 *  Remember that a created data object is automatically captured. If the user
 *  doesn't need the object anymore, he has to call the object's release() method.
 *
 *  When a data object is added to SCIP, it is captured again, such that a
 *  release() call does not destroy the object. If SCIP doesn't need the object
 *  anymore, it is automatically released.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page PARAM Adding additional user parameters
 *
 *  The user may add own parameters to SCIP with a call to SCIPaddXxxParam(). Using
 *  these methods, he has two possibilities where to store the actual parameter value:
 *   - If the given valueptr is NULL, SCIP stores the parameter value internally, and
 *     the user can only access the value with the SCIPgetXxxParam() and
 *     SCIPsetXxxParam() calls.
 *   - If the given valueptr is not NULL, SCIP stores the parameter value at the given
 *     address, and the user can directly manipulate the value at this address.
 *     He has to be careful with memory management in string parameters: when the
 *     SCIPaddStringParam() method is called, the given address must hold a char*
 *     pointer with value NULL. The default value is then copied into this pointer,
 *     allocating memory with BMSallocMemoryArray(). If the parameter is changed, the
 *     old string is freed with BMSfreeMemoryArray() and the new one is copied to a new
 *     memory area allocated with BMSallocMemoryArray(). When the parameter is freed,
 *     the memory is freed with BMSfreeMemoryArray().
 *     The user should not interfere with this internal memory management. Accessing
 *     the string parameter through the given valueptr is okay as long as it does not
 *     involve reallocating memory for the string.
 */

