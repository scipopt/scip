/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2008 Konrad-Zuse-Zentrum                            */
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
 * @author Timo Berthold
 * @author Mathias Kinder
 * @author Marc Pfetsch
 * @author Kati Wolter
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage SCIP (Solving Constraint Integer Programs)
 * @version  1.1.0
 * @author   Tobias Achterberg
 * @author   Timo Berthold
 * @author   Gerald Gamrath
 * @author   Stefan Heinz
 * @author   Thorsten Koch
 * @author   Alexander Martin
 * @author   Marc Pfetsch
 * @author   Christian Raack
 * @author   Michael Winkler
 * @author   Kati Wolter
 *
 * SCIP is a program and library to solve constraint integer programs (CIPs).
 *
 * SCIP is based on SIP (Solving Integer Programs) by Alexander Martin. The main developer of SCIP
 * was Tobias Achterberg (2002-2007) with contributions by Timo Berthold and Kati Wolter. The
 * persons listed above have contributed or are currently contributing to SCIP.
 *
 * <b>General Information</b>
 *
 * - \ref FAQ     "Frequently asked questions (FAQ)"
 * - \ref START   "How to start a new project"
 * - \ref DOC     "How to search the documentation for interface methods"
 * - \ref MAKE    "Makefiles"
 * - \ref DEBUG   "Debugging"
 * - \ref TEST    "How to run automated tests with SCIP"
 *
 * <b>Programming with SCIP</b>
 *
 * - \ref CODE    "Coding style guidelines"
 * - \ref CONS    "How to add constraint handlers"
 * - \ref PRICER  "How to add variable pricers"
 * - \ref PRESOL  "How to add presolvers"
 * - \ref SEPA    "How to add separators"
 * - \ref PROP    "How to add propagators"
 * - \ref BRANCH  "How to add branching rules"
 * - \ref NODESEL "How to add node selectors"
 * - \ref HEUR    "How to add primal heuristics"
 * - \ref RELAX   "How to add relaxation handlers"
 * - \ref READER  "How to add file readers"
 * - \ref DIALOG  "How to add dialogs"
 * - \ref DISP    "How to add display columns"
 * - \ref OBJ     "Creating, capturing, releasing, and adding data objects"
 * - \ref PARAM   "Adding additional user parameters"
 *
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
 * - No spaces between control structure keywords like "if", "for", "while", "switch" and the corresponding brackets.
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
/**@page MAKE Makefiles
 *
 * SCIP contains a makefile system, which allows the individual setting of several parameters. For
 * instance, the following settings are supported:
 *
 * - <code>OPT=\<dbg|opt\></code> Here <code>dbg</code> turns on the debug mode of SCIP. This enables
 *   asserts and avoids macros for several function in order to ease debugging. The default is
 *   <code>opt</code>, which enables the optimized mode.  
 *
 * - <code>LPS=\<clp|cpx|msk|spx|xprs|none\></code> This determines the LP-solver, which should have been
 *   installed separately from SCIP. The options are the following:
 *      - <code>clp</code>: COIN Clp LP-solver
 *      - <code>cpx</code>: CPLEX LP-solver
 *      - <code>msk</code>: Mosek LP-solver
 *      - <code>spx</code>: SoPlex LP-solver (default)
 *      - <code>xprs</code>: XPress LP-solver
 *      - <code>none</code>: no LP-solver (you should set the parameter \<lp/solvefreq\> to \<-1\> to avoid solving LPs)
 * - <code>LPSOPT=\<opt|dbg\></code> Chooses the optimized or debug version of the LP-solver. (currently only available for SoPlex and CLP)
 *
 * - <code>ZIMPL=\<true|false\></code> Turns direct support of ZIMPL in SCIP on or off, respectively.
 * - <code>ZIMPLOPT=\<opt|dbg\></code> Chooses the optimized or debug version of ZIMPL, if ZIMPL support is enabled.
 * 
 * - <code>READLINE=\<true|false\></code> Turns support via the readline library on or off, respectively.
 *
 * There are additional parameters for Linux/Gnu compilers:
 * - <code>OPT=noblkmem</code> turns off the internal SCIP memory.
 *   This way the code can be check via valgrind or similar tools.
 * - <code>OPT=opt-shared</code> generates a shared object of the SCIP libraries.
 *   (The binary uses these shared libraries as well.)
 * - <code>OPT=prf</code> generates a profiling version of SCIP 
 *   providing a detailed statistic of the time usage of every method of SCIP.
 *
 * The SCIP makefiles are structured as follows.
 *
 * - <code>Makefile</code> This is the basic makefile in the SCIP root directory. It loads
 *   additional makefile information depending on the parameters set.
 * - <code>make/make.project</code> This file contains definitions that are useful for all codes
 *   that use SCIP, for instance, the examples.
 * - <code>make.\<sys\>.\<machine\>.\<compiler\>.\<dbg|opt|prf\></code> These file contain system/compiler specific
 *   definitions. If you have a yet unsupported compiler, you could copy one of these and modify it
 *   accordingly.
 *
 * If your platform or compiler is not supported by SCIP you might try and copy one of the existing
 * makefile in the <code>make</code> directory and modify it. If you succeed, we are always
 * interested in including more Makefiles into the system.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page START How to start a new project
 *
 * If you want to use SCIP to write your own branch-and-cut or branch-and-cut-and-price code, below
 * you find some hints of how to get started.
 *
 * - Copy one of the examples in the <code>examples</code> directory (in the SCIP root
 *   directory). For branch-and-cut, the linear ordering example <code>LOP</code> should be a good
 *   starting point. If you use column generation have a look at <code>SamplePricer</code>.
 * - Edit the makefile according to your needs - in particular, include a correct path to the SCIP
 *   root at the top. Moreover, you should rename the targets and source file names.
 * - Once you have edited the makefile, you can use all the flags that can be used in SCIP to
 *   compile your code, see \ref MAKE.
 *
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
 * \n
 * A complete list of all constraint handles contained in this release can be found \ref CONSHDLRS "here".
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
 * Here is what you have to do (assuming your constraint handler should be named "subtour"):
 * -# Copy the template files "src/scip/cons_xxx.c" and "src/scip/cons_xxx.h" into files "cons_subtour.c"
 *    and "cons_subtour.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "subtour".
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
 * Additionally, if you are searching for a constraint handler with SCIPfindConshdlr(), this name is looked up.
 * Names have to be unique: no two constraint handlers may have the same name.
 *
 * \par CONSHDLR_DESC: the description of the constraint handler.
 * This string is printed as description of the constraint handler in the interactive shell of SCIP.
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
 * A separation frequency of -1 disables the separation method of the constraint handler.
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
 * callback of the constraint handler.
 * A propagation frequency of 0 means that propagation is only applied in preprocessing and at the root node.
 * A propagation frequency of -1 disables the propagation method of the constraint handler.
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
 * If the separation method of the constraint handler is very expensive, you may want to mark it to be delayed until all
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
 * Below the header "Data structures" you can find two structs called "struct SCIP_ConsData" and
 * "struct SCIP_ConshdlrData".
 * If you are using C++, you only need to define the "struct SCIP_ConsData".
 * The constraint handler data must be implemented as member variables of your constraint handler class.
 * \n
 * The constraint data are the information that is needed to define a single constraint of the constraint handler's
 * constraint class.
 * For example, the data of a knapsack constraint would consist of a list of variables, a list of weights, and
 * the capacity of the knapsack.
 * The data of a subtour constraint consists of the graph on which the problem is defined.
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
 * The method SCIPincludeConshdlrSubtour() has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the constraint handler by calling the method
 * SCIPincludeConshdlr().
 * It is called by the user, if he wants to include the constraint handler, i.e., if he wants to make
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
 * The method SCIPcreateConsSubtour() is called to create a single constraint of the constraint handler's constraint
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
 *    SCIP_Bool             removable,
 *    SCIP_Bool             stickingatnode
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
 *          local, modifiable, dynamic, removable, stickingatnode) );
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
 *          SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
 *          SCIPconsIsStickingAtNode(sourcecons)) );
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
 * changes, i.e., the bounds of variables at the time of the propagation, that forced the constraint to set the
 * conflict variable's bound to its current value. It can use the "inferinfo" tag to identify its own propagation
 * rule and thus identify the "reason" bounds. The bounds that form the reason of the assignment must then be provided
 * by calls to SCIPaddConflictLb() and SCIPaddConflictUb() in the propagation conflict resolving method.
 *
 * For example, the logicor constraint c = "x or y or z" fixes variable z to TRUE (i.e., changes the lower bound of z
 * to 1.0), if both, x and y, are assigned to FALSE (i.e., if the upper bounds of these variables are 0.0). It uses
 * SCIPinferVarLbCons(scip, z, 1.0, c, 0) to apply this assignment (an inference information tag is not needed by the
 * constraint handler and is set to 0).
 * In the conflict analysis, the constraint handler may be asked to resolve the lower bound change on z with
 * constraint c, that was applied at a time given by a bound change index "bdchgidx".
 * With a call to SCIPvarGetLbAtIndex(z, bdchgidx), the handler can find out, that the lower bound of variable z was
 * set to 1.0 at the given point of time, and should call SCIPaddConflictUb(scip, x, bdchgidx) and
 * SCIPaddConflictUb(scip, y, bdchgidx) to tell SCIP, that the upper bounds of x and y at this point of time were
 * the reason for the deduction of the lower bound of z.
 *
 * If conflict analysis should not be supported, the method has to set the result code to SCIP_DIDNOTFIND.
 * Although this is viable approach to circumvent the implementation of the usually rather complex conflict resolving mehod,
 * it will make the conflict analysis less effective. We suggest to first omit the conflict resolving method and check
 * how effective the propagation method is. If it produces a lot of propagations for your application, you definitely should
 * consider to implement the conflict resolving method.
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
 * The CONSPRINT callback method is called, when the user asks SCIP to display the problem to the screen or save
 * the problem into a file.
 * The constraint handler should display the data of the constraint in an appropriate form.
 * The output format that is defined by the CONSPRINT callbacks is called CIP format.
 * In later versions of SCIP, the constraint handlers should also be able to parse (i.e., read) constraints
 * which are given in CIP format.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page PRICER How to add variable pricers
 *
 * A pricer performs the dynamic generation of new variables in a column generation algorithm.
 * It is an algorithmic representation of a (usually exponential) number of variables.
 * The \ref PRICERREDCOST and \ref PRICERFARKAS methods are called after each LP solve to generate additional
 * variables which may improve the objective value or decrease the LP infeasibility, respectively. 
 * \n
 * A complete list of all pricers contained in this release can be found \ref PRICERS "here".
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
 * Additionally, if you are searching for a pricer with SCIPfindPricer(), this name is looked up.
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
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_PricerData".
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
 * It is called by the user, if he wants to include the pricer, i.e., if he wants to solve a model for which variables should
 * be generated by this pricer.
 *
 * If you are using pricer data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &pricerdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_PricerData afterwards.
 *
 * You may also add user parameters for your pricer, see the method SCIPincludeConshdlrKnapsack() in the knapsack constraint handler
 * src/scip/cons_knapsack.c for an example of how to add user parameters.
 *
 * 
 * @section PRICER_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Pricer
 *
 * The fundamental callback methods have to be implemented in order to obtain an operational algorithm.
 * In the case of a pricer, the fundamental callbacks are the two variable pricing callbacks which search and add
 * new variables to the problem.
 * In the C++ wrapper class ObjPricer, the scip_redcost() method (which corresponds to the PRICERREDCOST callback) is a virtual
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
 * Pricing in this case means to add variables \f$i\f$ with positive farkas value, i.e., \f$y^T A_i x'_i > 0\f$.
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
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page PRESOL How to add presolvers
 *
 * Presolvers are used to reduce the size of the model by removing irrelevant information like redundant constraints,
 * to strengthen the LP relaxation by exploiting integrality information, and to extract useful information in the 
 * presolving step.
 * Constraint based presolving is done in the CONSPRESOL callback methods of the constraint handlers, see \ref CONSPRESOL.
 * The presolver plugins complement the constraint based presolving by additional, usually optimality based, presolving
 * reductions. 
 * \n 
 * A complete list of all presolvers contained in this release can be found \ref PRESOLVERS "here".
 *
 * In the following, we explain how the user can add its own presolver.
 * Take the dual fixing presolver (src/scip/presol_dualfix.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the ObjPresol wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_PRESOL... callback methods.
 *
 * Additional documentation for the callback methods of a presolver, in particular for their input parameters, 
 * can be found in the file "type_presol.h".
 *
 * Here is what you have to do to implement a presolver:
 * -# Copy the template files "src/scip/presol_xxx.c" and "src/scip/presol_xxx.h" into files named "presol_mypresolver.c"
 *    and "presol_mypresolver.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "mypresolver".
 * -# Adjust the properties of the presolver (see \ref PRESOL_PROPERTIES).
 * -# Define the presolver data (see \ref PRESOL_DATA). This is optional.
 * -# Implement the interface methods (see \ref PRESOL_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref PRESOL_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref PRESOL_ADDITIONALCALLBACKS). This is optional.
 *
 * 
 * @section PRESOL_PROPERTIES Properties of a Presolver
 *
 * At the top of the new file "presol_mypresolver.c", you can find the presolver properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the presolver properties by calling the constructor
 * of the abstract base class ObjPresol from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par PRESOL_NAME: the name of the presolver.
 * This name is used in the interactive shell to address the presolver.
 * Additionally, if you are searching for a presolver with SCIPfindPresol(), this name is looked up.
 * Names have to be unique: no two presolvers may have the same name.
 *
 * \par PRESOL_DESC: the description of the presolver.
 * This string is printed as description of the presolver in the interactive shell.
 *
 * \par PRESOL_PRIORITY: the priority of the presolver.
 * In each presolving round, the presolvers and presolving methods of the constraint handlers are called in
 * a predefined order, which is given by the priorities of the presolvers and the check priorities of the
 * constraint handlers, see \ref CONS_PROPERTIES.
 * First, the presolvers with non-negative priority are called in the order of decreasing priority.
 * Next, the presolving methods of the different constraint handlers are called in the order of decreasing check
 * priority.
 * Finally, the presolvers with negative priority are called in the order of decreasing priority.
 * \n
 * The priority of the presolver should be set according to the complexity of the presolving algorithm and the impact of the reduction:
 * presolvers that provide fast algorithms that usually have a high impact (i.e., remove lots of variables or tighten 
 * bounds of many variables) should have a high priority. An easy way to list the 
 * priorities of all presolvers and constraint handlers is to type "display presolvers" and "display conshdlrs" in 
 * the interactive shell of SCIP.
 *
 * \par PRESOL_MAXROUNDS: the default maximal number of rounds the presolver participates in.
 * The presolving is conducted in rounds: the presolvers and presolving methods of the constraint handlers
 * are called iteratively until no more reductions have been found or some other abort criterion applies.
 * The "maxrounds" parameter of a presolver imposes a limit on the number of presolving rounds in which the
 * presolver is called. The PRESOL_MAXROUNDS property specifies the default value for this parameter.
 * A value of -1 represents an unlimited number of rounds.
 *
 * \par PRESOL_DELAY: the default for whether the presolver should be delayed, if other presolvers found reductions.
 * If the presolver is marked to be delayed, it is only executed if no other presolvers found a reduction during the current
 * presolving round.
 * If the presolver is very expensive, you may want to mark it to be delayed until all cheap presolving methods have been executed.
 * 
 *
 * @section PRESOL_DATA Presolver Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_PresolData".
 * In this data structure, you can store the data of your presolver. For example, you should store the adjustable parameters
 * of the presolver in this data structure.  
 * If you are using C++, you can add presolver data as usual as object variables to your class.
 * \n
 * Defining presolver data is optional. You can leave this struct empty.
 *
 *
 * @section PRESOL_INTERFACE Interface Methods
 *
 * At the bottom of "presol_mypresolver.c", you can find the interface method SCIPincludePresolMypresolver(), 
 * which also appears in "presol_mypresolver.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the presolver by calling the method
 * SCIPincludePresol(). 
 * SCIPincludePresolMypresolver() is called by the user, if he wants to include the presolver, 
 * i.e., if he wants to use the presolver in his application.
 *
 * If you are using presolver data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &presoldata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_PresolData afterwards. For freeing the 
 * presolver data, see \ref PRESOLFREE.
 *
 * You may also add user parameters for your presolver, see \ref PARAM for how to add user parameters and 
 * the method SCIPincludePresolProbing() in src/scip/presol_probing.c for an example.
 *
 * 
 * @section PRESOL_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Presolver
 *
 * The fundamental callback methods of the plugins are the ones that have to be implemented in order to obtain 
 * an operational algorithm. Presolver plugins have only one fundamental callback method, namely the PRESOLEXEC method.
 * This method has to be implemented for every presolver; the other callback methods are optional.
 * In the C++ wrapper class ObjPresol, the scip_exec() method (which corresponds to the PRESOLEXEC callback) is a virtual
 * abstract member function.
 * You have to implement it in order to be able to construct an object of your presolver class.
 *
 * Additional documentation to the callback methods, in particular to their input parameters, 
 * can be found in "type_presol.h".
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
 * If you are using presolver data (see \ref PRESOL_DATA and \ref PRESOL_INTERFACE), you have to implement this method in order to free the presolver data.
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
 * The difference between the original and the transformed problem is explained in 
 * "What is this thing with the original and the transformed problem about?" on \ref FAQ. 
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
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page SEPA How to add separators
 *
 * Separators are used to generate general purpose cutting planes. 
 * Constraint based cutting planes, the second type of cutting planes in SCIP, are separated in the CONSSEPALP and 
 * CONSSEPASOL callback methods of the constraint handlers, see \ref CONSSEPALP and \ref CONSSEPASOL. These cuts are 
 * valid inequalities or even facets of the polyhedron described by a single constraint or a subset of the constraints of
 * a single constraint class. In contrast, general purpose cuts do not require or exploit any knowledge about the 
 * underlying problem structure but use only the current LP relaxation and the integrality conditions. See also
 * "When should I implement a constraint handler, when should I implement a separator?" on \ref FAQ.
 * \n
 * A complete list of all separators contained in this release can be found \ref SEPARATORS "here".
 *
 * In the following, we explain how the user can add its own separator.
 * Take the separator for the class of Gomory mixed integer inequalities (src/scip/sepa_gomory.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the ObjSepa wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_SEPA... callback methods.
 *
 * Additional documentation for the callback methods of a separator, in particular for the input parameters, 
 * can be found in the file "type_sepa.h".
 *
 * Here is what you have to do to implement a separator:
 * -# Copy the template files "src/scip/sepa_xxx.c" and "src/scip/sepa_xxx.h" into files "sepa_myseparator.c"
 *    and "sepa_myseparator.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "myseparator".
 * -# Adjust the properties of the separator (see \ref SEPA_PROPERTIES).
 * -# Define the separator data (see \ref SEPA_DATA). This is optional.
 * -# Implement the interface methods (see \ref SEPA_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref SEPA_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref SEPA_ADDITIONALCALLBACKS).  This is optional.
 *
 *
 * @section SEPA_PROPERTIES Properties of a Separator
 *
 * At the top of the new file "sepa_myseparator.c", you can find the separator properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the separator properties by calling the constructor
 * of the abstract base class ObjSepa from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par SEPA_NAME: the name of the separator.
 * This name is used in the interactive shell to address the separator.
 * Additionally, if you are searching for a separator with SCIPfindSepa(), this name is looked up.
 * Names have to be unique: no two separators may have the same name.
 *
 * \par SEPA_DESC: the description of the separator.
 * This string is printed as description of the separator in the interactive shell.
 *
 * \par SEPA_PRIORITY: the priority of the separator.
 * In each separation round during the price-and-cut loop of the subproblem processing or the separation loop
 * of the primal solution separation, the separators and separation methods of the constraint handlers are called in
 * a predefined order, which is given by the priorities of the separators and the separation priorities 
 * of the constraint handlers (see \ref CONS_PROPERTIES).
 * First, the separators with non-negative priority are called in the order of decreasing priority.
 * Next, the separation methods of the constraint handlers are called in the order of decreasing separation
 * priority.
 * Finally, the separators with negative priority are called in the order of decreasing priority. An easy way to list the 
 * priorities of all separators and constraint handlers is to type "display separators" and "display conshdlrs" in 
 * the interactive shell.
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
 * A frequency of -1 disables the separator.
 * \n
 * The frequency can be adjusted by the user. The property of the separator only defines the default value of the frequency.
 * If you want to have a more flexible control of when to execute the separation algorithm, you have to assign
 * a frequency of 1 and implement a check at the beginning of your separation methods whether you really want to execute 
 * the separation or not. If you do not want to execute it, set the result code of 
 * \ref SEPAEXECLP and \ref SEPAEXECSOL to SCIP_DIDNOTRUN.
 *
 * \par SEPA_MAXBOUNDDIST: the default maximal relative distance from the current node's dual bound to primal bound compared to best node's dual bound for applying separation.
 * At the current branch-and-bound node, the relative distance from its dual bound (local dual bound) 
 * to the primal bound compared to the best node's dual bound (global dual bound) is considered. The separation method 
 * of the separator will only be applied at the current node if this relative distance does not exceed SEPA_MAXBOUNDDIST. 
 * \n
 * For example, if the global dual bound is 50 and the primal bound is 60, SEPA_MAXBOUNDDIST = 0.25 means that separation 
 * is only applied if the current node's dual bound is in the first quarter of the interval [50,60], i.e., if it is less 
 * than or equal to 52.5.
 * \n
 * In particular, the values 0.0 and 1.0 mean that separation is applied at the current best node only or at all 
 * nodes, respectively. Since separation seems to be most important to apply at nodes that define to the global 
 * dual bound, 0.0 is probably a good choice for SEPA_MAXBOUNDDIST.
 * Note that separators with a frequency of SEPA_FREQ = 0 are only applied at the root node.
 * Obviously, at the root node the local dual bound is equal to the global dual bound and thus, the separator is called
 * for any value of SEPA_MAXBOUNDDIST.
 *
 * \par SEPA_DELAY: the default for whether the separation method should be delayed, if other separators or constraint handlers found cuts.
 * If the separator's separation method is marked to be delayed, it is only executed after no other separator
 * or constraint handler found a cut during the price-and-cut loop. 
 * If the separation method of the separator is very expensive, you may want to mark it to be delayed until all cheap 
 * separation methods have been executed.
 *
 * @section SEPA_DATA Separator Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_SepaData".
 * In this data structure, you can store the data of your separator. For example, you should store the adjustable 
 * parameters of the separator in this data structure. In a separator, user parameters for the maximal number of 
 * separation rounds per node and for the maximal number of cuts separated per separation round might be useful. 
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
 * SCIPincludeSepaMyseparator() is called by the user, if he wants to include the separator, i.e., if he wants to use 
 * the separator in his application.
 *
 * If you are using separator data, you have to allocate the memory 
 * for the data at this point. You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &sepadata) );
 * \endcode
 * You also have to initialize the fields in "struct SCIP_SepaData" afterwards. For freeing the 
 * separator data, see \ref SEPAFREE.
 *
 * You may also add user parameters for your separator, see \ref PARAM for how to add user parameters and 
 * the method SCIPincludeSepaGomory() in src/scip/sepa_gomory.c for an example.
 *
 * 
 * @section SEPA_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Separator
 *
 * The fundamental callback methods of the plugins are the ones that have to be implemented in order to obtain 
 * an operational algorithm. Separator plugins do not have any fundamental callback methods, i.e., all 
 * callback methods are optional. But, most probably, you want to use at least the optional callback \ref SEPAEXECLP
 * to separate LP solutions.
 *
 * Additional documentation to the callback methods, in particular to their input parameters, 
 * can be found in "type_sepa.h".
 *
 * @section SEPA_ADDITIONALCALLBACKS Additional Callback Methods of a Separator
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 * The main functions of separators are the \ref SEPAEXECLP and \ref SEPAEXECSOL callbacks.
 *
 * @subsection SEPAEXECLP
 *
 * The SEPAEXECLP callback is executed during the price-and-cut loop of the subproblem processing.
 * It should try to generate general purpose cutting planes in order to separate the current LP solution.
 * The method is called in the LP solution loop, which means that a valid LP solution exists.
 *
 * Usually, the callback searches and produces cuts, that are added with a call to SCIPaddCut().
 * If the cut should be added to the global cut pool, it calls SCIPaddPoolCut().
 * In addition to LP rows, the callback may also produce domain reductions or add additional constraints.
 *
 * Overall, the SEPAEXECLP callback has the following options, which is indicated by the possible return values of
 * the 'result' variable (see "type_sepa.h"):
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
 * The SEPAEXECSOL callback is executed during the separation loop on arbitrary primal solutions.
 * It should try to generate general purpose cutting planes in order to separate the given primal solution.
 * The method is not called in the LP solution loop, which means that there is no valid LP solution.
 *
 * In the standard SCIP environment, the SEPAEXECSOL callback is not used because only LP solutions are
 * separated. The SEPAEXECSOL callback provides means to support external relaxation handlers like semidefinite
 * relaxations that want to separate an intermediate primal solution vector. Thus, if you do not want to support
 * such external plugins, you do not need to implement this callback method.
 *
 * Usually, the callback searches and produces cuts, that are added with a call to SCIPaddCut().
 * If the cut should be added to the global cut pool, it calls SCIPaddPoolCut().
 * In addition to LP rows, the callback may also produce domain reductions or add other constraints.
 *
 * Overall, the SEPAEXECSOL callback has the following options, which is indicated by the possible return values of
 * the 'result' variable (see "type_sepa.h"):
 *  - detecting that the node is infeasible in the variable's bounds and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint (result SCIP_CONSADDED)
 *  - reducing a variable's domain (result SCIP_REDUCEDDOM)
 *  - adding a cutting plane to the LP (result SCIP_SEPARATED)
 *  - stating that the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the separator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the separator was skipped, but should be called again (result SCIP_DELAYED)
 *
 * @subsection SEPAFREE
 *
 * If you are using separator data (see \ref SEPA_DATA and \ref SEPA_INTERFACE), you have to implement this method 
 * in order to free the separator data. This can be done by the following procedure:
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
 * The difference between the original and the transformed problem is explained in 
 * "What is this thing with the original and the transformed problem about?" on \ref FAQ. 
 *
 * @subsection SEPAEXIT
 *
 * The SEPAEXIT callback is executed before the transformed problem is freed.
 * In this method, the separator should free all resources that have been allocated for the solving process in SEPAINIT.
 *
 * @subsection SEPAINITSOL
 *
 * The SEPAINITSOL callback is executed when the presolving was finished and the branch-and-bound process is about to 
 * begin. The separator may use this call to initialize its branch-and-bound specific data.
 *
 * @subsection SEPAEXITSOL
 *
 * The SEPAEXITSOL callback is executed before the branch-and-bound process is freed. The separator should use this call 
 * to clean up its branch-and-bound data, in particular to release all LP rows that he has created or captured.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page PROP How to add propagators
 *
 * Propagators are used to tighten the domains of the variables. Like for cutting planes, there are two different types of 
 * domain propagations. Constraint based (primal) domain propagation algorithms are part of the corresponding constraint
 * handlers, see \ref CONSPROP. In contrast, domain propagators usually provide dual propagations, i.e., propagations that can be
 * applied due to the objective function and the current best known primal solution. 
 * \n
 * A complete list of all propagators contained in this release can be found \ref PROPAGATORS "here".
 *
 * In the following, we explain how the user can add an own propagator.
 * Take the pseudo objective function propagator (src/scip/prop_pseudoobj.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the ObjProp wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_PROP... callback methods.
 *
 * Additional documentation for the callback methods of a propagator can be found in the file "type_prop.h".
 *
 * Here is what you have to do to implement a propagator:
 * -# Copy the template files "src/scip/prop_xxx.c" and "src/scip/prop_xxx.h" into files named "prop_mypropagator.c"
 *    and "prop_mypropagator.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "mypropagator".
 * -# Adjust the properties of the propagator (see \ref PROP_PROPERTIES).
 * -# Define the propagator data (see \ref PROP_DATA).
 * -# Implement the interface methods (see \ref PROP_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref PROP_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref PROP_ADDITIONALCALLBACKS).
 *
 * @section PROP_PROPERTIES Properties of a Propagator
 *
 * At the top of the new file "prop_mypropagator.c" you can find the propagator properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the propagator properties by calling the constructor
 * of the abstract base class ObjProp from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par PROP_NAME: the name of the propagator.
 * This name is used in the interactive shell to address the propagator.
 * Additionally, if you are searching for a propagator with SCIPfindProp(), this name is looked up.
 * Names have to be unique: no two propagators may have the same name.
 *
 * \par PROP_DESC: the description of the propagator.
 * This string is printed as description of the propagator in the interactive shell.
 *
 * \par PROP_PRIORITY: the priority of the propagator.
 * In each propagation round, the propagators and propagation methods of the constraint handlers are called in
 * a predefined order, which is given by the priorities of the propagators and the check priorities of the
 * constraint handlers.
 * First, the propagators with non-negative priority are called in the order of decreasing priority.
 * Next, the propagation methods of the different constraint handlers are called in the order of decreasing check
 * priority.
 * Finally, the propagators with negative priority are called in the order of decreasing priority.
 * \n
 * The priority of the propagators should be set according to the complexity of the propagation algorithm and the impact 
 * of the domain propagations: propagators that provide fast algorithms that usually have a high impact (i.e., tighten 
 * many bounds) should have a high priority.
 *
 * \par PROP_FREQ: the default frequency for propagating domains.
 * The frequency defines the depth levels at which the propagation method \ref PROPEXEC is called.
 * For example, a frequency of 7 means, that the propagation callback is executed for subproblems that are in depth 
 * 0, 7, 14, ... of the branching tree. A frequency of 0 means that propagation is only applied in preprocessing and 
 * at the root node. A frequency of -1 disables the propagator.
 * \n
 * The frequency can be adjusted by the user. The property of the propagator only defines the default value of the 
 * frequency. If you want to have a more flexible control of when to execute the propagation algorithm, you have to assign
 * a frequency of 1 and implement a check at the beginning of your propagation algorithm whether you really 
 * want to execute the domain propagation or not. If you do not want to execute it, set the result code to SCIP_DIDNOTRUN.
 *
 * \par PROP_DELAY: the default for whether the propagation method should be delayed, if other propagators or constraint handlers found domain reductions?
 * If the propagator's propagation method is marked to be delayed, it is only executed after no other propagator or 
 * constraint handler found a domain reduction in the current iteration of the domain propagation loop.
 * If the propagation method of the propagator is very expensive, you may want to mark it to be delayed until all cheap 
 * propagation methods have been executed.
 *
 *
 * @section PROP_DATA Propagator Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_PropData".
 * In this data structure, you can store the data of your propagator. For example, you should store the adjustable 
 * parameters of the propagator in this data structure.
 * If you are using C++, you can add propagator data as usual as object variables to your class.
 * \n
 * Defining propagator data is optional. You can leave the struct empty.
 * 
 *
 * @section PROP_INTERFACE Interface Methods
 *
 * At the bottom of "prop_mypropagator.c" you can find the interface method SCIPincludePropMypropagator(), which also 
 * appears in "prop_mypropagator.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the propagator by calling the method
 * SCIPincludeProp().
 * It is called by the user, if he wants to include the propagator, i.e., if he wants to use the propagator in his 
 * application.
 *
 * If you are using propagator data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &propdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_PropData afterwards.
 *
 * You may also add user parameters for your propagator, see the method SCIPincludePropPseudoobj() in 
 * src/scip/prop_pseudoobj.c for an example.
 *
 *
 * @section PROP_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Propagator
 *
 * Propagator plugins have two fundamental callback methods, namely the PROPEXEC method and the PROPRESPROP method.
 * These methods have to be implemented for every propagator; the other callback methods are optional.
 * In the C++ wrapper class ObjProp, the scip_exec() method and the scip_resprop() method (which correspond to the 
 * PROPEXEC callback and PROPRESPROP callback, respectively) are virtual abstract member functions.
 * You have to implement them in order to be able to construct an object of your propagator class.
 *
 * Additional documentation to the callback methods can be found in "type_prop.h".
 *
 * @subsection PROPEXEC
 *
 * The PROPEXEC callback is called during presolving and during the subproblem processing. 
 * It should perform the actual domain propagation, which means that it should tighten the variables' bounds.
 * The technique of domain propagation, which is the main workhorse of constraint programming, is called "node preprocessing" in the
 * Integer Programming community.
 *
 * The PROPEXEC callback has the following options:
 *  - detecting that the node is infeasible in the variable's bounds and can be cut off (result SCIP_CUTOFF)
 *  - reducing (i.e, tightening) the domains of some variables (result SCIP_REDUCEDDOM)
 *  - stating that the propagator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the propagator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the propagator was skipped, but should be called again (result SCIP_DELAYED)
 * 
 * @subsection PROPRESPROP
 *
 * If the propagator wants to support conflict analysis, it has to supply the PROPRESPROP method.
 * It also should call SCIPinferVarLbProp() or SCIPinferVarUbProp() in the domain propagation instead of SCIPchgVarLb() 
 * or SCIPchgVarUb() in order to deduce bound changes on variables.
 * In the SCIPinferVarLbProp() and SCIPinferVarUbProp() calls, the propagator provides a pointer to itself and an integer 
 * value "inferinfo" that can be arbitrarily chosen.
 *  
 * The propagation conflict resolving method PROPRESPROP must then be implemented, to provide the "reasons" for the bound
 * changes, i.e., the bounds of variables at the time of the propagation, that forced the propagator to set the
 * conflict variable's bound to its current value. It can use the "inferinfo" tag to identify its own propagation
 * rule and thus identify the "reason" bounds. The bounds that form the reason of the assignment must then be provided
 * by calls to SCIPaddConflictLb() and SCIPaddConflictUb() in the propagation conflict resolving method.
 *
 * See the description of the propagation conflict resolving method \ref CONSRESPROP of constraint handlers for 
 * further details.
 *
 * If conflict analysis should not be supported, the method has to set the result code to SCIP_DIDNOTFIND.
 * Although this is viable approach to circumvent the implementation of the usually rather complex conflict resolving mehod,
 * it will make the conflict analysis less effective. We suggest to first omit the conflict resolving method and check
 * how effective the propagation method is. If it produces a lot of propagations for your application, you definitely should
 * consider to implement the conflict resolving method.
 *
 *
 * @section PROP_ADDITIONALCALLBACKS Additional Callback Methods of a Propagator
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 *
 * @subsection PROPFREE
 *
 * If you are using propagator data, you have to implement this method in order to free the propagator data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_PROPFREE(propFreeMypropagator)
 * { 
 *    SCIP_PROPDATA* propdata;
 * 
 *    propdata = SCIPpropGetData(prop);
 *    assert(propdata != NULL);
 *
 *   SCIPfreeMemory(scip, &propdata);
 *
 *   SCIPpropSetData(prop, NULL);
 *
 *   return SCIP_OKAY;
 * }
 * \endcode
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection PROPINIT
 *
 * The PROPINIT callback is executed after the problem was transformed.
 * The propagator may, e.g., use this call to initialize his propagator data.
 *
 * @subsection PROPEXIT
 *
 * The PROPEXIT callback is executed before the transformed problem is freed.
 * In this method, the propagator should free all resources that have been allocated for the solving process in PROPINIT.
 *
 * @subsection PROPINITSOL
 *
 * The PROPINITSOL callback is executed when the presolving was finished and the branch and bound process is about to
 * begin.
 * The propagator may use this call to initialize its branch and bound specific data.
 * 
 * @subsection PROPEXITSOL
 *
 * The PROPEXITSOL callback is executed before the branch and bound process is freed.
 * The propagator should use this call to clean up its branch and bound data.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page BRANCH How to add branching rules
 *
 * Branching rules are used to split the current problem into subproblems. If the LP solution of the current problem is 
 * fractional, the integrality constraint handler calls the branching rules. Additionally, branching rules are called 
 * as a last resort on integral solutions that violate one or more constraints for which the associated constraint handlers
 * were not able to resolve the infeasibility in a more sophisticated way, see the constraint handlers' callback methods 
 * \ref CONSENFOLP and \ref CONSENFOPS. Note that in these resolving methods the constraint handlers may also apply 
 * constraint specific branching rules, like the so-called special ordered set branching of the set 
 * partitioning/packing/covering constraint handler.
 * \n
 * Usually, a branching rule creates two subproblems by splitting a single variable's domain. It is also possible to 
 * implement much more general branching schemes, for example by creating more than two subproblems, or by adding 
 * additional constraints to the subproblems instead of tightening the domains of the variables. 
 * \n 
 * A complete list of all branching rules contained in this release can be found \ref BRANCHINGRULES "here".
 *
 * In the following, we explain how the user can add an own branching rule.
 * Take the most infeasible LP branching rule (src/scip/branch_mostinf.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the ObjBranchrule wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_BRANCH... callback methods.
 *
 * Additional documentation for the callback methods of a branching rule can be found in the file "type_branch.h".
 *
 * Here is what you have to do to implement a branching rule:
 * -# Copy the template files "src/scip/branch_xxx.c" and "src/scip/branch_xxx.h" into files named 
 *    "branch_mybranchingrule.c" and "branch_mybranchingrule.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "mybranchingrule".
 * -# Adjust the properties of the branching rule (see \ref BRANCHRULE_PROPERTIES).
 * -# Define the branching rule data (see \ref BRANCHRULE_DATA).
 * -# Implement the interface methods (see \ref BRANCHRULE_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref BRANCHRULE_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref BRANCHRULE_ADDITIONALCALLBACKS).
 *
 *
 * @section BRANCHRULE_PROPERTIES Properties of a Branching Rule
 *
 * At the top of the new file "branch_mybranchingrule.c" you can find the branching rule properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the branching rule properties by calling the constructor
 * of the abstract base class ObjBranchrule from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par BRANCHRULE_NAME: the name of the branching rule.
 * This name is used in the interactive shell to address the branching rule.
 * Additionally, if you are searching for a branching rule with SCIPfindBranchrule(), this name is looked up.
 * Names have to be unique: no two branching rules may have the same name.
 *
 * \par BRANCHRULE_DESC: the description of the branching rule.
 * This string is printed as description of the branching rule in the interactive shell.
 *
 * \par BRANCHRULE_PRIORITY: the default value for the priority of the branching rule.
 * In the subproblem processing, the branching rules are called in decreasing order of their priority until 
 * one succeeded to branch. Since most branching rules are able to generate a branching in all situations,
 * only the rule of highest priority is used. In combination with the BRANCHRULE_MAXDEPTH and
 * BRANCHRULE_MAXBOUNDDIST settings, however, interesting strategies can be easily employed. For example,
 * the user can set the priority of the "full strong branching" strategy to the highest value and assign the
 * second highest value to the "reliable pseudo cost" rule. If he also sets the maximal depth for the
 * "full strong branching" to 5, in the top 5 depth levels of the search tree the "full strong branching" is
 * applied, while in the deeper levels "reliable pseudo cost branching" is used.
 * \n
 * Note that the BRANCHRULE_PRIORITY property only specifies the default value of the priority. The user can
 * change this value arbitrarily.
 *
 * \par BRANCHRULE_MAXDEPTH: the default value for the maximal depth level of the branching rule.
 * This parameter denotes the maximal depth level in the branch-and-bound tree up to which the branching method of the 
 * branching rule will be applied. Use -1 for no limit.
 * \n
 * Note that this property only specifies the default value. The user can change this value arbitrarily.
 *
 * \par BRANCHRULE_MAXBOUNDDIST: the default value for the maximal relative distance from current node's dual bound to 
 * primal bound compared to best node's dual bound for applying branching.
 * At the current branch-and-bound node, the relative distance from its dual bound (local dual bound) 
 * to the primal bound compared to the best node's dual bound (global dual bound) is considered. The branching method of 
 * the branching rule will only be applied at the node if this relative distance does not exceed BRANCHRULE_MAXBOUNDDIST. 
 * \n
 * For example, if the global dual bound is 50 and the primal bound is 60, BRANCHRULE_MAXBOUNDDIST = 0.25 means that 
 * branching is only applied if the current node's dual bound is in the first quarter of the interval [50,60], i.e., if it 
 * is less than or equal to 52.5. In particular, the values 0.0 and 1.0 mean that the branching rule is applied at the
 * current best node only or at all nodes, respectively.
 * \n
 * Note that the BRANCHRULE_MAXBOUNDDIST property only specifies the default value of the maximal bound distance.
 * The user can change this value arbitrarily.
 * 
 *
 * @section BRANCHRULE_DATA Branching Rule Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_BranchruleData".
 * In this data structure, you can store the data of your branching rule. For example, you should store the adjustable 
 * parameters of the branching rule in this data structure.
 * If you are using C++, you can add branching rule data as usual as object variables to your class.
 * \n
 * Defining branching rule data is optional. You can leave the struct empty.
 *
 *
 * @section BRANCHRULE_INTERFACE Interface Methods
 *
 * At the bottom of "branch_mybranchingrule.c" you can find the interface method SCIPincludeBranchruleMybranchingrule(), 
 * which also appears in "branch_mybranchingrule.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the branching rule by calling the method
 * SCIPincludeBranchrule().
 * It is called by the user, if he wants to include the branching rule, i.e., if he wants to use the branching rule in his 
 * application.
 *
 * If you are using branching rule data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_BranchruleData afterwards.
 *
 * You may also add user parameters for your branching rule, see the method SCIPincludeBranchruleRelpscost() in 
 * src/scip/branch_relpscost.c for an example.
 *
 * 
 * @section BRANCHRULE_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Branching Rule
 *
 * Branching rules do not have any fundamental callback methods, i.e., all callback methods are optional.
 * In most cases, however, you want to implement the \ref BRANCHEXECLP method and sometimes the \ref BRANCHEXECPS method.
 *
 *
 * @section BRANCHRULE_ADDITIONALCALLBACKS Additional Callback Methods of a Branching Rule
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 * The most important callback methods are the \ref BRANCHEXECLP and \ref BRANCHEXECPS methods, which perform the
 * actual task of generating a branching.
 *
 * Additional documentation to the callback methods can be found in "type_branch.h".
 *
 * @subsection BRANCHEXECLP
 *
 * The BRANCHEXECLP callback is executed during node processing if a fractional LP solution is available. It should 
 * split the current problem into subproblems. Usually, the branching is done in a way such that the current fractional
 * LP solution is no longer feasible in the relaxation of the subproblems.
 * It is, however, possible to create a child node for which the fractional LP solution is still feasible in the relaxation,
 * for example, by branching on a variable with integral LP value.
 * In every case, you have to make sure that each subproblem is a proper restriction of the current problem.
 * Otherwise, you risk to produce an infinite path in the search tree.
 *
 * The user gains access to the branching candidates, i.e., to the fractional variables, and their LP solution values
 * by calling the method SCIPgetLPBranchCands(). Furthermore, SCIP provides two methods for performing the actual 
 * branching, namely SCIPbranchVar() and SCIPcreateChild(). 
 *
 * Given an integral variable \f$x\f$ with fractional LP solution value 
 * \f$x^*\f$, the method SCIPbranchVar() creates two child nodes; one contains the bound \f$x \le \lfloor x^* \rfloor\f$ 
 * and the other one contains the bound \f$x \ge \lceil x^* \rceil\f$, see the BRANCHEXECLP callback in 
 * src/scip/branch_mostinf.c for an example. In addition, if a proven lower objective bound of a created child node is
 * known, like after strong branching has been applied, the user may call the method SCIPupdateNodeLowerbound() in order
 * to update the child node's lower bound.   
 *
 * In order to apply more general branching schemes, one should use the method SCIPcreateChild(). For each child node
 * that the branching rule wants to generate, it has to call SCIPcreateChild() once. The branching rule has
 * to assign two values to the new nodes: a node selection priority for each node and an estimate for the objective value 
 * of the best feasible solution contained in the subtree after applying the branching. If the method SCIPbranchVar() 
 * is used, these values are automatically assigned. Here, the user can calculate these values by calling the method 
 * SCIPcalcNodeselPriority(), which takes into account the preferred branch direction of the branching variable, and
 * the method SCIPcalcChildEstimate(), which is based on pseudo costs, and pass them to the SCIPcreateChild() call.     
 * \n
 * After having created a child node, the additional restrictions of the child node have to be added with calls
 * to SCIPaddConsNode(), SCIPchgVarLbNode(), and SCIPchgVarUbNode().
 *
 * In some cases, the branching rule can tighten the current subproblem instead of producing a branching.
 * Therefore, the BRANCHEXECLP callback may also produce domain reductions or add additional constraints to the current subproblem.
 *
 * The BRANCHEXECLP callback has the following options:
 *  - detecting that the node is infeasible and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint (e.g. a conflict constraint) (result SCIP_CONSADDED; note that this action
 *    must not be performed if the input "allowaddcons" is FALSE)
 *  - reducing the domain of a variable (result SCIP_REDUCEDDOM)
 *  - adding a cutting plane to the LP (result SCIP_SEPARATED)
 *  - applying a branching (result SCIP_BRANCHED)
 *  - stating that the branching rule was skipped (result SCIP_DIDNOTRUN), which means that the branching rule with the
 *    next largest priority is called.
 *
 * @subsection BRANCHEXECPS
 *
 * The BRANCHEXECPS callback is executed during node processing if no LP solution is available and at least one of the
 * integer variables is not yet fixed. It should split the current problem into subproblems.
 *
 * The user gains access to the branching candidates, i.e., to the non-fixed integer variables, by calling the method 
 * SCIPgetPseudoBranchCands(). Furthermore, SCIP provides two methods for performing the actual branching, namely 
 * SCIPbranchVar() and SCIPcreateChild(). 
 *
 * Given an integral variable \f$x\f$ with pseudo solution value \f$x^*\f$, the method SCIPbranchVar() creates 
 * three child nodes. Two of them contain the bounds \f$x \le x^* - 1 \f$ and \f$x \ge x^* + 1 \f$, respectively, such that the current 
 * pseudo solution is cut off. In the third node, the variable \f$x\f$ is fixed to \f$ x^*\f$ and this hopefully reduces 
 * other variables' domains and allows to cut off the current pseudo solution. See the BRANCHEXECPS callback in 
 * src/scip/branch_random.c for an example. In addition, if a proven lower bound of a created child node is known the user 
 * may call the method SCIPupdateNodeLowerbound() in order to update the child node's lower bound.   
 *
 * In order to apply more general branching schemes, one should use the method SCIPcreateChild(). The branching rule has
 * to assign two values to the new nodes: a node selection priority for each node and an estimate for the objective value 
 * of the best feasible solution contained in the subtree after applying the branching. If the method SCIPbranchVar() 
 * is used, these values are automatically assigned. Here, the user can calculate these values by calling the method 
 * SCIPcalcNodeselPriority(), which takes into account the preferred branch direction of the branching variable, and
 * the method SCIPcalcChildEstimate(), which is based on pseudo costs, and pass them to the SCIPcreateChild() call.     
 * \n
 * After having created a child node, the additional restrictions of the child node have to be added with calls
 * to SCIPaddConsNode(), SCIPchgVarLbNode(), and SCIPchgVarUbNode().
 *
 * In some cases, the branching rule can tighten the current subproblem instead of producing a branching. For example,
 * strong branching might have proven that rounding up a variable would lead to an infeasible LP relaxation and thus,
 * the variable must be rounded down. Therefore, the BRANCHEXECPS callback may also produce domain reductions or add
 * additional constraints to the current subproblem.
 *
 * The BRANCHEXECPS callback has the following options:
 *  - detecting that the node is infeasible and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint (e.g. a conflict constraint) (result SCIP_CONSADDED; note that this action
 *    must not be performed if the input "allowaddcons" is FALSE)
 *  - reducing the domain of a variable that rendered the current LP solution to be infeasible (result SCIP_REDUCEDDOM)
 *  - applying a branching (result SCIP_BRANCHED)
 *  - stating that the branching rule was skipped (result SCIP_DIDNOTRUN).
 *
 * Note that in contrast to the BRANCHEXECLP callback, the BRANCHEXECPS callback cannot add cutting planes to the current
 * LP relaxation.
 *
 * @subsection BRANCHFREE
 *
 * If you are using branching rule data, you have to implement this method in order to free the branching rule data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_BRANCHFREE(branchFreeMybranchingrule)
 * {
 *    SCIP_BRANCHRULEDATA* branchruledata;
 * 
 *    branchruledata = SCIPbranchruleGetData(branchrule);
 *    assert(branchruledata != NULL);
 *    
 *    SCIPfreeMemory(scip, &branchruledata);
 *  
 *    SCIPbranchruleSetData(branchrule, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection BRANCHINIT
 *
 * The BRANCHINIT callback is executed after the problem was transformed.
 * The branching rule may, e.g., use this call to initialize his branching rule data.
 *
 * @subsection BRANCHEXIT
 *
 * The BRANCHEXIT callback is executed before the transformed problem is freed.
 * In this method, the branching rule should free all resources that have been allocated for the solving process in 
 * BRANCHINIT.
 *
 * @subsection BRANCHINITSOL
 *
 * The BRANCHINITSOL callback is executed when the presolving was finished and the branch and bound process is about to
 * begin.
 * The branching rule may use this call to initialize its branch and bound specific data.
 *
 * @subsection BRANCHEXITSOL
 *
 * The BRANCHEXITSOL callback is executed before the branch and bound process is freed.
 * The branching rule should use this call to clean up its branch and bound data.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page NODESEL How to add node selectors
 *
 * Node selectors are used to decide which of the leaves in the current branching tree is selected as next subproblem 
 * to be processed. The ordering relation of the tree's leaves for storing them in the leave priority queue is also
 * defined by the node selectors.  
 * \n
 * A complete list of all propagators contained in this release can be found \ref NODESELECTORS "here".
 *
 * In the following, we explain how the user can add an own node selector.
 * Take the node selector for depth first search (src/scip/nodesel_dfs.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the ObjNodesel wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_NODESEL... callback methods.
 *
 * Additional documentation for the callback methods of a node selector can be found in the file "type_nodesel.h".
 *
 * Here is what you have to do to implement a node selector:
 * -# Copy the template files "src/scip/nodesel_xxx.c" and "src/scip/nodesel_xxx.h" into files named "nodesel_mynodeselector.c"
 *    and "nodesel_mynodeselector.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "mynodeselector".
 * -# Adjust the properties of the node selector (see \ref NODESEL_PROPERTIES).
 * -# Define the node selector data (see \ref NODESEL_DATA).
 * -# Implement the interface methods (see \ref NODESEL_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref NODESEL_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref NODESEL_ADDITIONALCALLBACKS).
 *
 *
 * @section NODESEL_PROPERTIES Properties of a Node Selector
 *
 * At the top of the new file "nodesel_mynodeselector.c" you can find the node selector properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the node selector properties by calling the constructor
 * of the abstract base class ObjNodesel from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par NODESEL_NAME: the name of the node selector.
 * This name is used in the interactive shell to address the node selector.
 * Additionally, if you are searching for a node selector with SCIPfindNodesel(), this name is looked up.
 * Names have to be unique: no two node selectors may have the same name.
 * 
 * \par NODESEL_DESC: the description of the node selector.
 * This string is printed as description of the node selector in the interactive shell.
 *
 * \par NODESEL_STDPRIORITY: the default priority of the node selector in the standard mode.
 *
 * The first step of each iteration of the main solving loop is the selection of the next subproblem to be processed. 
 * The node selector of highest priority (the active node selector) is called to do this selection. 
 * Note that SCIP has two different operation modes: the standard mode and the memory saving mode. If the memory 
 * limit - given as a parameter by the user - is almost reached, SCIP switches from the standard mode to the memory saving 
 * mode in which different priorities for the node selectors are applied. NODESEL_STDPRIORITY is the priority of the 
 * node selector used in the standard mode.
 * \n
 * Note that this property only defines the default value of the priority. The user may change this value arbitrarily by
 * adjusting the corresponding parameter setting.
 *
 * \par NODESEL_MEMSAVEPRIORITY: the default priority of the node selector in the memory saving mode.
 *
 * The priority NODESEL_MEMSAVEPRIORITY of the node selector has the same meaning as the priority NODESEL_STDPRIORITY, but 
 * is used in the memory saving mode.
 * Usually, you want the best performing node selector, for example best estimate search, to have maximal
 * standard priority, while you want a node selector which tends to keep the growth of the search tree limited, for example
 * depth first search, to have maximal memory saving priority.
 * \n
 * Note that this property only defines the default value of the priority. The user may change this value arbitrarily by
 * adjusting the corresponding parameter setting.
 *
 *
 * @section NODESEL_DATA Node Selector Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_NodeselData".
 * In this data structure, you can store the data of your node selector. For example, you should store the adjustable 
 * parameters of the node selector in this data structure.
 * If you are using C++, you can add node selector data as usual as object variables to your class.
 * \n
 * Defining node selector data is optional. You can leave the struct empty.
 *
 *
 * @section NODESEL_INTERFACE Interface Methods
 *
 * At the bottom of "nodesel_mynodeselector.c" you can find the interface method SCIPincludeNodeselMynodeselector(), 
 * which also appears in "nodesel_mynodeselector.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the node selector by calling the method
 * SCIPincludeNodesel().
 * It is called by the user, if he wants to include the node selector, i.e., if he wants to use the node selector in 
 * his application.
 *
 * If you are using node selector data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &nodeseldata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_NodeselData afterwards.
 *
 * You may also add user parameters for your node selector, see the method SCIPincludeNodeselRestartdfs() in 
 * src/scip/nodesel_restartdfs.c for an example.
 *
 * 
 * @section NODESEL_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Node Selector
 *
 * Node selector plugins have two fundamental callback methods, namely the NODESELSELECT method and the NODESELCOMP method.
 * These methods have to be implemented for every node selector; the other callback methods are optional.
 * In the C++ wrapper class ObjNodesel, the scip_select() method and the scip_comp() method (which correspond to the 
 * NODESELSELECT callback and the NODESELCOMP callback, respectively) are virtual abstract member functions.
 * You have to implement them in order to be able to construct an object of your node selector class.
 *
 * Additional documentation to the callback methods can be found in "type_nodesel.h".
 *
 * @subsection NODESELSELECT
 *
 * The NODESELSELECT callback is the first method called in each iteration in the main solving loop. It should decide 
 * which of the leaves in the current branching tree is selected as next subproblem to be processed. It can decide between 
 * the current node's children and siblings, and the "best" of the remaining leaves stored in the tree. This choice can 
 * have a large impact on the solver's performance, because it influences the finding of feasible solutions and the 
 * development of the global dual bound. 
 *
 * The following methods provide access to the various leaf nodes:
 * - SCIPgetPrioChild() returns the child of the current node with the largest node selection priority, as assigned by the
 *   branching rule, see the \ref BRANCHEXECLP and \ref BRANCHEXECPS callbacks of the branching rules. If no child is
 *   available (for example, because the current node was pruned), a NULL pointer is returned.
 * - SCIPgetBestChild() returns the best child of the current node with respect to the node selector's ordering relation as
 *   defined by the \ref NODESELCOMP callback. If no child is available, a NULL pointer is returned.
 * - SCIPgetPrioSibling() returns the sibling of the current node with the largest node selection priority.
 *   If no sibling is available (for example, because all siblings of the current node have already been processed), a NULL
 *   pointer is returned.
 * - SCIPgetBestSibling() returns the best sibling of the current node with respect to the node selector's ordering relation
 *   as defined by the \ref NODESELCOMP callback. If no sibling is available, a NULL pointer is returned.
 * - SCIPgetBestNode() returns the best leaf from the tree (children, siblings, or other leaves) with respect to the node
 *   selector's ordering relation as defined by the \ref NODESELCOMP callback. If no open leaf exists, a NULL pointer is
 *   returned. In this case, the optimization is finished, and the node selector should return a NULL pointer as 'selnode'.
 * - SCIPgetBestboundNode() returns a leaf from the tree (children, siblings, or other leaves) with the smallest lower (dual)
 *   objective bound. If the queue is empty, a NULL pointer is returned. In this case, the optimization is finished, and the
 *   node selector should return a NULL pointer as 'selnode'.
 *   
 *
 * @subsection NODESELCOMP
 *
 * The NODESELCOMP callback is called to compare two leaves of the current branching tree (say node 1 and node 2) 
 * regarding their ordering relation.
 *
 * The NODESELCOMP should return the following values:
 *  - value < 0, if node 1 comes before (is better than) node 2
 *  - value = 0, if both nodes are equally good
 *  - value > 0, if node 2 comes after (is worse than) node 2.
 *
 * @section NODESEL_ADDITIONALCALLBACKS Additional Callback Methods of a Node Selector
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 *
 * @subsection NODESELFREE
 *
 * If you are using node selector data, you have to implement this method in order to free the node selector data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_NODESELFREE(nodeselFreeMynodeselector)
 * { 
 *    SCIP_NODESELDATA* nodeseldata;
 *
 *    nodeseldata = SCIPnodeselGetData(nodesel);
 *    assert(nodeseldata != NULL);
 *
 *    SCIPfreeMemory(scip, &nodeseldata);
 *
 *    SCIPnodeselSetData(nodesel, NULL);
 * 
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection NODESELINIT
 *
 * The NODESELINIT callback is executed after the problem was transformed.
 * The node selector may, e.g., use this call to initialize his node selector data.
 *
 * @subsection NODESELEXIT
 *
 * The NODESELEXIT callback is executed before the transformed problem is freed.
 * In this method, the node selector should free all resources that have been allocated for the solving process 
 * in NODESELINIT.
 *
 * @subsection NODESELINITSOL
 *
 * The NODESELINITSOL callback is executed when the presolving was finished and the branch and bound process is about to
 * begin.
 * The node selector may use this call to initialize its branch and bound specific data.
 *
 * @subsection NODESELEXITSOL
 *
 * The NODESELEXITSOL callback is executed before the branch and bound process is freed.
 * The node selector should use this call to clean up its branch and bound data.
 */


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page HEUR How to add primal heuristics
 *
 * Feasible solutions can be found in two different ways during the traversal of the branch-and-bound tree. On the one 
 * hand, the solution of a node's relaxation may be feasible with respect to the constraints. On the other hand, feasible
 * solutions can be discovered by primal heuristics.  
 * \n
 * A complete list of all primal heuristics contained in this release can be found \ref PRIMALHEURISTICS "here".
 *
 * In the following, we explain how the user can add an own primal heuristic.
 * Take the simple and fast LP rounding heuristic (src/scip/heur_simplerounding.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the ObjHeur wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_HEUR... callback methods.
 *
 * Additional documentation for the callback methods of a primal heuristic can be found in the file "type_heur.h".
 *
 * Here is what you have to do to implement a primal heuristic:
 * -# Copy the template files "src/scip/heur_xxx.c" and "src/scip/heur_xxx.h" into files named "heur_myheuristic.c"
 *    and "heur_myheuristic.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "myheuristic".
 * -# Adjust the properties of the primal heuristic (see \ref HEUR_PROPERTIES).
 * -# Define the primal heuristic data (see \ref HEUR_DATA).
 * -# Implement the interface methods (see \ref HEUR_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref HEUR_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref HEUR_ADDITIONALCALLBACKS).
 *
 *
 * @section HEUR_PROPERTIES Properties of a Primal Heuristic
 *
 * At the top of the new file "heur_myheuristic.c" you can find the primal heuristic properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the primal heuristic properties by calling the constructor
 * of the abstract base class ObjHeur from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par HEUR_NAME: the name of the primal heuristic.
 * This name is used in the interactive shell to address the primal heuristic.
 * Additionally, if you are searching for a primal heuristic with SCIPfindHeur(), this name is looked up.
 * Names have to be unique: no two primal heuristics may have the same name.
 *
 * \par HEUR_DESC: the description of the primal heuristic.
 * This string is printed as description of the primal heuristic in the interactive shell.
 *
 * \par HEUR_DISPCHAR: the display character of the primal heuristic.
 * In the interactive shell, this character is printed in the first column of a status information row, if the primal 
 * heuristic found the feasible solution belonging to the primal bound. Note that a star stands for an integral 
 * LP-relaxation.
 * In order to avoid confusion, display characters should be unique: no two primal heuristics should have the same display character.
 * You can get a list of all primal heuristics along with their display characters by entering "display heuristics" in the
 * SCIP interactive shell.
 *
 * \par HEUR_PRIORITY: the priority of the primal heuristic.
 * At each of the different entry points of the primal heuristics during the solving process (see HEUR_TIMING), they are 
 * called in decreasing order of their priority. 
 * \n
 * The priority of a primal heuristic should be set according to the complexity of the heuristic and the likelihood to find
 * feasible solutions: primal heuristics that provide fast algorithms that often succeed in finding a feasible solution should have
 * a high priority. In addition, the interaction between different types of primal heuristics should be taken into account.
 * For example, improvement heuristics, which try to generate improved solutions  by inspecting one or more of the feasible
 * solutions that have already been found, should have a small priority.
 *
 * \par HEUR_FREQ: the default frequency for executing the primal heuristic.
 * The frequency together with the frequency offset (see HEUR_FREQOFS) defines the depth levels at which the execution
 * method of the primal heuristic \ref HEUREXEC is called. For example, a frequency of 7 together with a frequence offset 
 * of 5 means, that the \ref HEUREXEC callback is executed for subproblems that are in depth 5, 12, 19, ... of the branching tree. A 
 * frequency of 0 together with a frequence offset of 3 means, that the execution method is only called at those nodes that are in
 * depth level 3 (i.e., at most for \f$2^3 = 8\f$ nodes if binary branching is applied). A frequency of 0 and an offset of 0 means that
 * the heuristic is only called at the root node. 
 * A frequency of -1 disables the heuristic.
 * \n
 * The frequency can be adjusted by the user. The property of the primal heuristic only defines the default value of the 
 * frequency. If you want to have a more flexible control of when to execute the primal heuristic, you have to assign
 * a frequency of 1 and implement a check at the beginning of your execution method whether you really want to search for feasible
 * solutions or not. If you do not want to execute the method, set the result code to SCIP_DIDNOTRUN.
 *
 * \par HEUR_FREQOFS: the frequency offset for executing the primal heuristic.
 * The frequency offset defines the depth of the branching tree at which the primal heuristic is executed for the first 
 * time. For example, a frequency of 7 (see HEUR_FREQ) together with a frequency offset of 10 means, that the 
 * callback is executed for subproblems that are in depth 10, 17, 24, ... of the branching tree. In particular, assigning 
 * different offset values to heuristics of the same type, like diving heuristics, can be useful for evenly spreading the 
 * application of these heuristics across the branch-and-bound tree.
 * Note that if the frequency is equal to 1, the heuristic is applied for all nodes with depth level larger or equal to
 * the frequency offset.
 *
 * \par HEUR_MAXDEPTH: the maximal depth level for executing the primal heuristic.
 * This parameter denotes the maximal depth level in the branching tree up to which the execution method of the primal 
 * heuristic is called. Use -1 for no limit. 
 *
 * \par HEUR_TIMING: the execution timing of the primal heuristic.
 * Primal heuristics have different entry points during the solving process and the execution timing parameter defines the
 * entry point at which the primal heuristic is executed first. 
 * \n
 * The primal heuristic can be called first:
 * - before the processing of the node starts (SCIP_HEURTIMING_BEFORENODE)
 * - after each LP solving during the cut-and-price loop (SCIP_HEURTIMING_DURINGLPLOOP) 
 * - after the cut-and-price loop was finished (SCIP_HEURTIMING_AFTERLPLOOP) 
 * - after the processing of a node <em>with solved LP</em>  was finished (SCIP_HEURTIMING_AFTERLPNODE)
 * - after the processing of a node <em>without solved LP</em> was finished (SCIP_HEURTIMING_AFTERPSEUDONODE)
 * - after the processing of the last node in the current plunge was finished, <em>and only if the LP was solved for 
 *   this node</em> (SCIP_HEURTIMING_AFTERLPPLUNGE) 
 * - after the processing of the last node in the current plunge was finished, <em>and only if the LP was not solved 
 *   for this node</em> (SCIP_HEURTIMING_AFTERPSEUDOPLUNGE).
 * \par
 * These flags can be combined as or concatenations to call the heuristic at multiple times. Two useful combinations
 * are already predefined:
 * - after the processing of a node was finished (SCIP_HEURTIMING_AFTERNODE; combines SCIP_HEURTIMING_AFTERLPNODE and
 *   SCIP_HEURTIMING_AFTERPSEUDONODE)
 * - after the processing of the last node in the current plunge was finished (SCIP_HEURTIMING_AFTERPLUNGE; combines
 *   SCIP_HEURTIMING_AFTERLPPLUNGE and SCIP_HEURTIMING_AFTERPSEUDOPLUNGE)
 * \par
 * Calling a primal heuristic "before the processing of the node starts" is particularly useful for heuristics 
 * that do not need to access the LP solution of the current node. If such a heuristic finds a feasible solution, the 
 * leaves of the branching tree exceeding the new primal bound are pruned. It may happen that even the current node can 
 * be cut off without solving the LP relaxation. Combinatorial heuristics, like the farthest insert heuristic for the TSP 
 * (see examples/TSP/src/HeurFarthestInsert.cpp), are often applicable at this point.
 * \n
 * Very fast primal heuristics that require an LP solution can also be called "after each LP solving during the 
 * cut-and-price loop". Rounding heuristics, like the simple and fast LP rounding heuristic 
 * (src/scip/heur_simplerounding.c), belong to this group of primal heuristics. 
 * \n
 * Most heuristics, however, are called either after the cut-and-price loop was finished, after a node was completely
 * processed, or even only after a full plunge was finished. These points correspond to the situation where a node was completely
 * processed, but differ by the type of node and by whether an LP solution is needed. For example, diving heuristics,
 * like the LP diving heuristic (src/scip/heur_fracdiving.c), are often executed after the last node in the current
 * plunge has been processed. A plunge is the successive solving of child and sibling nodes in the search tree. A plunge is
 * finished, if the node selection selects an unrelated node from the queue of open subproblems, see \ref NODESEL.
 *
 * Computational experiments seem to indicate that for the overall performance of a MIP solver, it is important to evenly 
 * spread the application of the heuristics across the branch-and-bound tree. Thus, the assignment of the parameters 
 * HEUR_FREQ, HEUR_FREQOFS, and HEUR_TIMING should contribute to this aim.
 *
 * Note that all diving heuristics in the SCIP distribution (see, e.g., src/scip/heur_guideddiving.c) check whether other diving
 * heuristics have already been called at the current node. This can be done by comparing SCIPgetLastDivenode(scip) and
 * SCIPgetNNodes(scip). If the two are equal, and if the current node is not the root node (SCIPgetDepth(scip) > 0), diving
 * heuristics should be delayed by returning the result code 'SCIP_DELAYED'. This is an additional contribution to the goal of
 * not calling multiple heuristics at the same node.
 *
 *
 * @section HEUR_DATA Primal Heuristic Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_HeurData".
 * In this data structure, you can store the data of your primal heuristic. For example, you should store the adjustable 
 * parameters of the primal heuristic in this data structure. 
 * If you are using C++, you can add primal heuristic data as usual as object variables to your class.
 * \n
 * Defining primal heuristic data is optional. You can leave the struct empty.
 *
 *
 * @section HEUR_INTERFACE Interface Methods
 *
 * At the bottom of "heur_myheuristic.c" you can find the interface method SCIPincludeHeurMyheuristic(), which also 
 * appears in "heur_myheuristic.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the primal heuristic by calling the method SCIPincludeHeur().
 * It is called by the user, if he wants to include the primal heuristic, i.e., if he wants to use the primal heuristic
 * in his application.
 *
 * If you are using primal heuristic data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_HeurData afterwards.
 *
 * You may also add user parameters for your primal heuristic, see the method SCIPincludeHeurFeaspump() in 
 * src/scip/heur_feaspump.c for an example.
 *
 * 
 * @section HEUR_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Primal Heuristic
 *
 * Primal heuristic plugins have only one fundamental callback method, namely the HEUREXEC method.
 * This method has to be implemented for every primal heuristic; the other callback methods are optional.
 * In the C++ wrapper class ObjHeur, the scip_exec() method (which corresponds to the HEUREXEC callback) is a virtual
 * abstract member function. You have to implement it in order to be able to construct an object of your primal heuristic 
 * class.
 *
 * Additional documentation to the callback methods can be found in "type_heur.h".
 *
 * @subsection HEUREXEC
 *
 * The HEUREXEC callback is called at different positions during the node processing loop, see HEUR_TIMING. It should
 * search for feasible solutions and add them to the solution pool. For creating a new feasible solution, the 
 * methods SCIPcreateSol() and SCIPsetSolVal() can be used. Afterwards, the solution can be added to the storage by 
 * calling the method SCIPtrySolFree() (or SCIPtrySol() and SCIPfreeSol()).
 *
 * The HEUREXEC callback has the following options:
 *  - finding at least one feasible solution (result SCIP_FOUNDSOL)
 *  - stating that the primal heuristic searched, but did not find a feasible solution (result SCIP_DIDNOTFIND)
 *  - stating that the primal heuristic was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the primal heuristic was skipped, but should be called again (result SCIP_DELAYED).
 *
 *
 * @section HEUR_ADDITIONALCALLBACKS Additional Callback Methods of a Primal Heuristic
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 *
 * @subsection HEURFREE
 *
 * If you are using primal heuristic data, you have to implement this method in order to free the primal heuristic data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_HEURFREE(heurFreeMyheuristic)
 * {
 *    SCIP_HEURDATA* heurdata;
 *  
 *    heurdata = SCIPheurGetData(heur);
 *    assert(heurdata != NULL);
 *
 *    SCIPfreeMemory(scip, &heurdata);
 *
 *    SCIPheurSetData(heur, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection HEURINIT
 *
 * The HEURINIT callback is executed after the problem was transformed.
 * The primal heuristic may, e.g., use this call to initialize his primal heuristic data.
 *
 * @subsection HEUREXIT
 *
 * The HEUREXIT callback is executed before the transformed problem is freed.
 * In this method, the primal heuristic should free all resources that have been allocated for the solving process in 
 * HEURINIT.
 *
 * @subsection HEURINITSOL
 *
 * The HEURINITSOL callback is executed when the presolving was finished and the branch and bound process is about to 
 * begin. The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 * @subsection HEUREXITSOL
 *
 * The HEUREXITSOL callback is executed before the branch and bound process is freed. The primal heuristic should use this
 * call to clean up its branch and bound data, which was allocated in HEURINITSOL.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page RELAX How to add relaxation handlers
 *
 * SCIP provides specific support for LP relaxations of constraint integer programs. In addition, relaxation handlers, 
 * also called relaxators, can be used to include other relaxations, e.g., Lagrange relaxations or semidefinite 
 * relaxations. The relaxation handler manages the necessary data structures and calls the relaxation solver to generate dual 
 * bounds and primal solution candidates.
 * \n
 * However, the data to define a single relaxation must either be extracted by the relaxation handler itself (e.g., from
 * the user defined problem data, the LP information, or the integrality conditions), or be provided by the constraint
 * handlers. In the latter case, the constraint handlers have to be extended to support this specific relaxation. 
 * \n
 * A complete list of all primal heuristics contained in this release can be found \ref RELAXATORS "here".
 *
 * In the following, we explain how the user can add an own relaxation handler using the C interface. It is very easy to 
 * transfer the C explanation to C++: whenever a method should be implemented using the SCIP_DECL_RELAX... notion, 
 * reimplement the corresponding virtual member function of the abstract ObjRelax wrapper base class.
 * Unfortunately, SCIP does not contain a default relaxation handler plugin, which could be used as an example.
 *
 * Additional documentation for the callback methods of a relaxation handler can be found in the file "type_relax.h".
 *
 * Here is what you have to do to implement a relaxation handler:
 * -# Copy the template files "src/scip/relax_xxx.c" and "src/scip/relax_xxx.h" into files named "relax_myrelaxator.c"
 *    and "relax_myrelaxator.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "myrelaxator".
 * -# Adjust the properties of the relaxation handler (see \ref RELAX_PROPERTIES).
 * -# Define the relaxation handler data (see \ref RELAX_DATA).
 * -# Implement the interface methods (see \ref RELAX_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref RELAX_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref RELAX_ADDITIONALCALLBACKS).
 *
 * 
 * @section RELAX_PROPERTIES Properties of a Relaxation Handler
 *
 * At the top of the new file "relax_myrelaxator.c" you can find the relaxation handler properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the relaxation handler properties by calling the constructor
 * of the abstract base class ObjRelax from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par RELAX_NAME: the name of the relaxation handler.
 * This name is used in the interactive shell to address the relaxation handler.
 * Additionally, if you are searching for a relaxation handler with SCIPfindRelax(), this name is looked up.
 * Names have to be unique: no two relaxation handlers may have the same name.
 *
 * \par RELAX_DESC: the description of the relaxation handler.
 * This string is printed as description of the relaxation handler in the interactive shell.
 *
 * \par RELAX_PRIORITY: the priority of the relaxation handler.
 * In each relaxation solving round during the subproblem processing, the included relaxation handlers and the 
 * price-and-cut loop for solving the LP relaxation are called in a predefined order, which is given by the priorities 
 * of the relaxation handlers. 
 * First, the relaxation handlers with non-negative priority are called in the order of decreasing priority.
 * Next, the price-and-cut loop for solving the LP relaxation is executed. 
 * Finally, the relaxation handlers with negative priority are called in the order of decreasing priority.
 * \n
 * Usually, you will have only one relaxation handler in your application and thus only have to decide whether it should 
 * be called before or after solving the LP relaxation. For this decision you should consider the complexity of 
 * the relaxation solving algorithm and the impact of the resulting solution: if your relaxation handler provides a fast 
 * algorithm that usually has a high impact (i.e., the relaxation is a good approximation of the convex hull of the 
 * feasible region of the subproblem and the solution severely reduces the primal-dual gap), it should have a non-negative 
 * priority.
 * \n
 * Note that for certain applications, it is useful to disable the LP relaxation and only use your custom relaxation.
 * This can easily be achieved by setting the "lp/solvefreq" parameter to -1.
 *
 * \par RELAX_FREQ: the default frequency for solving the relaxation.
 * The frequency defines the depth levels at which the relaxation solving method \ref RELAXEXEC is called.
 * For example, a frequency of 7 means, that the relaxation solving callback is executed for subproblems that are in depth 
 * 0, 7, 14, ... of the branching tree. A frequency of 0 means that the callback is only executed at the root node, i.e., 
 * only the relaxation of the root problem is solved. A frequency of -1 disables the relaxation handler.
 *
 *
 * @section RELAX_DATA Relaxation Handler Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_RelaxData".
 * In this data structure, you can store the data of your relaxation handler. For example, you should store the adjustable
 * parameters of the relaxation handler in this data structure.
 * If you are using C++, you can add relaxation handler data as usual as object variables to your class.
 * \n
 * Defining relaxation handler data is optional. You can leave the struct empty.
 *
 *
 * @section RELAX_INTERFACE Interface Methods
 *
 * At the bottom of "relax_myrelaxator.c" you can find the interface method SCIPincludeRelaxMyrelaxator(), which also 
 * appears in "relax_myrelaxator.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the relaxation handler by calling the method SCIPincludeRelax().
 * It is called by the user, if he wants to include the relaxation handler, i.e., if he wants to use the relaxation 
 * handler in his application.
 *
 * If you are using relaxation handler data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_RelaxData afterwards.
 *
 * You may also add user parameters for your relaxation handler, see the method SCIPincludeConshdlrKnapsack() in 
 * the knapsack constraint handler src/scip/cons_knapsack.c for an example of how to add user parameters.
 *
 *
 * @section RELAX_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Relaxation Handler
 *
 * Relaxation handler plugins have only one fundamental callback method, namely the \ref RELAXEXEC method.
 * This method has to be implemented for every relaxation handler; the other callback methods are optional.
 * In the C++ wrapper class ObjRelax, the scip_exec() method (which corresponds to the \ref RELAXEXEC callback) is a virtual
 * abstract member function.
 * You have to implement it in order to be able to construct an object of your relaxation handler class.
 *
 * Additional documentation to the callback methods can be found in "type_relax.h".
 *
 * @subsection RELAXEXEC
 * The RELAXEXEC is called in each relaxation solving round during node processing. It should solve the current 
 * subproblem's relaxation.  
 *
 * Note that, like the LP relaxation, the relaxation handler should only operate on variables for which the corresponding 
 * column does exist in the transformed problem. Typical methods called by a relaxation handler are SCIPconstructLP() to
 * make sure that the LP of the current node is constructed and its data can be accessed via calls to SCIPgetLPRowsData()
 * and SCIPgetLPColsData(), SCIPseparateSol() to call the cutting plane separators for a given primal solution, and
 * SCIPupdateLocalLowerbound() to update the current node's dual bound after having solved the relaxation.
 * In addition, you may want to call SCIPtrySolFree() if you think that you have found a feasible primal solution.
 *
 * Note that the primal solution of the relaxation cannot be stored inside the data structures of SCIP, which means in
 * particular, that the branching rules cannot take the solution as a guide on how to split the problem into subproblems.
 * If you want to branch with respect to your relaxation solution, you have to implement your own branching rule and
 * extract the primal solution vector from the relaxation directly.
 *
 * Usually, the callback only solves the relaxation and provides a lower (dual) bound with a call to SCIPupdateLocalLowerbound().
 * However, it may also produce domain reductions, add additional constraints or generate cutting planes. It has the
 * following options:
 *  - detecting that the node is infeasible in the variable's bounds and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint and stating that the relaxation handler should not be called again on the same 
 *    relaxation (result SCIP_CONSADDED)
 *  - reducing a variable's domain and stating that the relaxation handler should not be called again on the same 
 *    relaxation (result SCIP_REDUCEDDOM)
 *  - adding a cutting plane to the LP and stating that the relaxation handler should not be called again on the same 
 *    relaxation (result SCIP_SEPARATED)
 *  - stating that the relaxation handler solved the relaxation and should not be called again on the same relaxation 
 *    (result SCIP_SUCCESS)
 *  - interrupting the solving process to wait for additional input, e.g., cutting planes (SCIP_SUSPENDED) 
 *  - stating that the separator was skipped (result SCIP_DIDNOTRUN).
 *
 * In the above criteria, "the same relaxation" means that the LP relaxation stayed unmodified. This means in particular
 * that no row has been added and no bounds have been modified. For example, changing the bounds of a variable will, as
 * long as it was a COLUMN variable, lead to a modification in the LP such that the relaxation handler is called again
 * after it returned with the result code SCIP_REDUCEDDOM.
 * 
 *
 * @section RELAX_ADDITIONALCALLBACKS Additional Callback Methods of a Relaxation Handler
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 *
 * @subsection RELAXFREE
 *
 * If you are using relaxation handler data, you have to implement this method in order to free the relaxation handler 
 * data. This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_RELAXFREE(relaxFreeMyrelaxator)
 * {
 *    SCIP_RELAXDATA* relaxdata;
 *  
 *    relaxdata = SCIPrelaxGetData(relax);
 *    assert(relaxdata != NULL);
 *
 *    SCIPfreeMemory(scip, &relaxdata);
 *
 *    SCIPrelaxSetData(relax, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection RELAXINIT
 *
 * The RELAXINIT callback is executed after the problem was transformed.
 * The relaxation handler may, e.g., use this call to initialize his relaxation handler data.
 *
 * @subsection RELAXEXIT
 *
 * The RELAXEXIT callback is executed before the transformed problem is freed.
 * In this method, the relaxation handler should free all resources that have been allocated for the solving process in 
 * RELAXINIT.
 *
 * @subsection RELAXINITSOL
 *
 * The RELAXINITSOL callback is executed when the presolving was finished and the branch and bound process is about to
 * begin. The relaxation handler may use this call to initialize its branch and bound specific data.
 *
 * @subsection REALXEXITSOL
 *
 * The RELAXEXITSOL callback is executed before the branch and bound process is freed.
 * The relaxation handler should use this call to clean up its branch and bound data.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page READER How to add file readers
 *
 * Mainly, file readers are called to parse an input file and generate a constraint integer programming model. They create 
 * constraints and variables and activate variable pricers if necessary. However, they can also be called, for example, to parse an 
 * input file containing information about a primal solution or fixing of variables. 
 * \n
 * A complete list of all file readers contained in this release can be found \ref FILEREADERS "here".
 *
 * In the following, we explain how the user can add an own file reader.
 * Take the file reader for MIPs in IBM's Mathematical Programming System format (src/scip/reader_mps.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the ObjReader wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_READER... callback methods.
 *
 * Additional documentation for the callback methods of a file reader can be found in the file "type_reader.h".
 *
 * Here is what you have to do to implement a file reader:
 * -# Copy the template files "src/scip/reader_xxx.c" and "src/scip/reader_xxx.h" into files named "reader_myreader.c"
 *    and "reader_myreader.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "myreader".
 * -# Adjust the properties of the file reader (see \ref READER_PROPERTIES).
 * -# Define the file reader data (see \ref READER_DATA).
 * -# Implement the interface methods (see \ref READER_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref READER_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref READER_ADDITIONALCALLBACKS).
 *
 * 
 * @section READER_PROPERTIES Properties of a File Reader
 *
 * At the top of the new file "reader_myreader.c" you can find the file reader properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the file reader properties by calling the constructor
 * of the abstract base class ObjReader from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par READER_NAME: the name of the file reader.
 * This name is used in the interactive shell to address the file reader.
 * Additionally, if you are searching for a file reader with SCIPfindReader(), this name is looked up.
 * Names have to be unique: no two file readers may have the same name.
 * 
 * \par READER_DESC: the description of the file reader.
 * This string is printed as description of the file reader in the interactive shell.
 *
 * \par READER_EXTENSION: the file name extension of the file reader. 
 * Each file reader is hooked to a single file name extension. It is automatically called if the user wants to read in a 
 * file of corresponding name. The extensions of the different file readers have to be unique.
 * Note that the additional extension '.gz', '.z', or '.Z' (indicating a gzip compressed file) are ignored for assigning
 * an input file to a reader.
 *
 *
 * @section READER_DATA File Reader Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_ReaderData".
 * In this data structure, you can store the data of your file reader. For example, you should store the adjustable 
 * parameters of the file reader in this data structure.
 * If you are using C++, you can add file reader data as usual as object variables to your class.
 * \n
 * Defining file reader data is optional. You can leave the struct empty.
 *
 *
 * @section READER_INTERFACE Interface Methods
 *
 * At the bottom of "reader_myreader.c" you can find the interface method SCIPincludeReaderMyreader(), which also 
 * appears in "reader_myreader.h".
 * It is responsible for notifying SCIP of the presence of the file reader by calling the method
 * SCIPincludeReader().
 * It is called by the user, if he wants to include the file reader, i.e., if he wants to use the file reader in his 
 * application.
 *
 * If you are using file reader data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &readerdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_ReaderData afterwards.
 *
 * You may also add user parameters for your file reader, see the method SCIPincludeReaderLp() in 
 * src/scip/reader_lp.c for an example.
 *
 *
 * @section READER_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a File Reader
 *
 * File reader plugins have only one fundamental callback method, namely the \ref READERREAD method.
 * This method has to be implemented for every file reader; the other callback method is optional.
 * In the C++ wrapper class ObjReader, the scip_read() method (which corresponds to the \ref READERREAD callback) is a virtual
 * abstract member function.
 * You have to implement it in order to be able to construct an object of your file reader class.
 *
 * Additional documentation to the callback methods can be found in "type_reader.h".
 *
 * @subsection READERREAD
 *
 * The READERREAD callback is called when the user invokes SCIP to read in a file with file name extension 
 * corresponding to the READER_EXTENSION property of the file reader. This is usually triggered by a call to the method
 * SCIPreadProb() or by an interactive shell command. The READERREAD callback should parse the input file and perform
 * the desired action, which usually means to generate a constraint integer programming model, to add a primal solution, or
 * to fix variables in an existing model.
 * \n
 * Typical methods called by a file reader which is used to generate constraint integer programming models are, 
 * for example, 
 *
 * - creating an empty problem: SCIPcreateProb() 
 * - creating the variables: SCIPcreateVar(), SCIPchgVarType(), SCIPchgVarLb(), SCIPchgVarUb(), SCIPaddVar(), and 
 *   SCIPreleaseVar()      
 * - modifying the objective function: SCIPchgVarObj() and SCIPsetObjsense().
 * - creating the constraints: SCIPcreateConsLinear(), SCIPaddCoefLinear(), SCIPchgLhsLinear(), SCIPchgRhsLinear(), 
 *   SCIPaddCons(), and SCIPreleaseCons()
 *
 * Primal solutions can only be created for the transformed problem. Therefore, the user has to call SCIPtransformProb()
 * before he reads in the file containing the solution and adds it to the solution pool via the method SCIPreadSol(). 
 *
 * 
 * @section READER_ADDITIONALCALLBACKS Additional Callback Methods of a File Reader
 *
 * File reader plugins have only one additional callback method, namely the READERFREE method. It needs not to be 
 * implemented in every case.
 * 
 * @subsection READERFREE
 *
 * If you are using file reader data, you have to implement this method in order to free the file reader data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_READERFREE(readerFreeMyreader)
 * {
 *    SCIP_READERDATA* readerdata;
 *  
 *    readerdata = SCIPreaderGetData(reader);
 *    assert(readerdata != NULL);
 *
 *    SCIPfreeMemory(scip, &readerdata);
 *
 *    SCIPreaderSetData(reader, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DIALOG How to add dialogs
 *
 * SCIP comes with a command line shell which allows the user to read in problem instances, modify the solver's 
 * parameters, initiate the optimization, and display certain statistics and solution information. This shell consists 
 * of dialogs, which are organized as a tree in SCIP. A node of this tree which is not a leaf represents a menu in 
 * the shell and the children of this node correspond to the entries of this menu (which can again be menus). All 
 * different dialogs are managed by a dialog handler, which, in particular, is responsible for executing the dialog 
 * corresponding to the user's command in the shell. That is, the concept of a dialog handler is different to that 
 * of a constraint handler, which is used to manage objects of the same structure, see \ref CONS. In particular, SCIP 
 * features only one dialog handler, whereas there may exist different constraint handlers. 
 * \n
 * A complete list of all dialogs contained in this release can be found \ref DIALOGS "here".
 *
 * In the following, we explain how the user can extend the interactive shell by adding an own dialog.
 * We give the explanation for creating an own source file for each additional dialog. Of course, you can collect 
 * different dialogs in one source file. Take "src/scip/dialog_default.c", where all default dialog plugins are collected, as an 
 * example.
 * As all other default plugins, the default dialog plungins and the template dialog are written in C. C++ users can easily
 * adapt the code by using the ObjDialog wrapper base class and implement the scip_...() virtual methods instead of the
 * SCIP_DECL_DIALOG... callback methods.
 *
 * Additional documentation for the callback methods of a dialog can be found in the file "type_dialog.h".
 *
 * Here is what you have to do to add a dialog:
 * -# Copy the template files "src/scip/dialog_xxx.c" and "src/scip/dialog_xxx.h" into files named "dialog_mydialog.c"
 *    and "dialog_mydialog.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "mydialog".
 * -# Adjust the properties of the dialog (see \ref DIALOG_PROPERTIES).
 * -# Define the dialog data (see \ref DIALOG_DATA).
 * -# Implement the interface methods (see \ref DIALOG_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref DIALOG_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref DIALOG_ADDITIONALCALLBACKS).
 *
 *
 * @section DIALOG_PROPERTIES Properties of a Dialog
 *
 * At the top of the new file "dialog_mydialog.c" you can find the dialog properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the dialog properties by calling the constructor
 * of the abstract base class ObjDialog from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par DIALOG_NAME: the name of the dialog.
 * In the interactive shell, this name appears as command name of the dialog in the parent dialog. 
 * Additionally, if you are searching an entry in a menu with SCIPdialogFindEntry(), this name is looked up.
 * Names within one menu have to be unique: no two dialogs in the same menu may have the same name.
 *
 * \par DIALOG_DESC: the description of the dialog.
 * This string is printed as description of the dialog in the interactive shell if the DIALOGDESC callback
 * is not implemented.
 *
 * \par DIALOG_ISSUBMENU: whether the dialog is a (sub)menu.
 * This parameter states whether the dialog is a menu in the interactive shell, i.e., is the parent of further 
 * dialogs.
 * 
 *
 * @section DIALOG_DATA Dialog Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_DialogData".
 * In this data structure, you can store the data of your dialog. 
 * If you are using C++, you can add dialog data as usual as object variables to your class.
 * \n
 * Defining dialog data is optional. You can leave the struct empty.
 *
 *
 * @section DIALOG_INTERFACE Interface Methods
 *
 * At the bottom of "dialog_mydialog.c" you can find the interface method SCIPincludeDialogMydialog(), which also appears
 * in "dialog_mydialog.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the dialog, which can be done by the following lines of code:
 * \code
 * if( !SCIPdialogHasEntry(parentdialog, DIALOG_NAME) )
 * {
 *    SCIP_CALL( SCIPcreateDialog(scip, &dialog, dialogExecXxx, dialogDescXxx, dialogFreeXxx,
 *          DIALOG_NAME, DIALOG_DESC, DIALOG_ISSUBMENU, dialogdata) );
 *
 *    SCIP_CALL( SCIPaddDialogEntry(scip, parentdialog, dialog) );
 * 
 *    SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
 * }
 * \endcode
 * Here "parentdialog" has to be an existing dialog which is defined to be a menu (see DIALOG_ISSUBMENU), e.g., 
 * the default root dialog.   
 *
 * The interface method is called by the user, if he wants to include the dialog, i.e., if he wants to use the dialog in 
 * his application. 
 * Note that in order to be able to link the new dialog to an existing default dialog it has to be included <b>after the 
 * default dialogs plugin</b>, i.e., the SCIPincludeDialogMydialog() call has to occure after the 
 * SCIPincludeDialogDefault() call. The SCIPincludeDialogDefault() method is called from within the SCIPincludeDefaultPlugins()
 * method. Therefore, it suffices to include your dialog plugins after you have called SCIPincludeDefaultPlugins().
 *
 * If you are using dialog data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &dialogdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_DialogData afterwards.
 *
 * Consider the following example. The user wants to add a "drawgraph" command to the root menu of SCIP.
 * He copies the "dialog_xxx.c" and "dialog_xxx.h" files into files "dialog_drawgraph.c" and "dialog_drawgraph.h", respectively.
 * Then, he puts the following code into the SCIPincludeDialogDrawgraph() method, compare SCIPincludeDialogDefault() in
 * src/scip/dialog_default.c:
 * \code
 * SCIP_RETCODE SCIPincludeDialogDrawgraph(
 *    SCIP*                 scip
 *    )
 * {
 *    SCIP_DIALOG* root;
 *    SCIP_DIALOG* dialog;
 * 
 *    root = SCIPgetRootDialog(scip);
 *    if( root == NULL )
 *    {
 *       SCIP_CALL( SCIPcreateDialog(scip, &root, SCIPdialogExecMenuLazy, NULL, NULL,
 *             "SCIP", "SCIP's main menu", TRUE, NULL) );
 *       SCIP_CALL( SCIPsetRootDialog(scip, root) );
 *       SCIP_CALL( SCIPreleaseDialog(scip, &root) );
 *       root = SCIPgetRootDialog(scip);
 *    }
 * 
 *    if( !SCIPdialogHasEntry(root, "drawgraph") )
 *    {
 *       SCIP_CALL( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDrawgraph, NULL, NULL,
 *             "drawgraph", "draws the graph for the current problem instance", FALSE, NULL) );
 *       SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
 *       SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
 *    }
 * 
 *    return SCIP_OKAY;
 * }
 * \endcode
 *
 * Using this code, it is even possible to call SCIPincludeDialogDrawgraph() before including the default dialog plugins,
 * and you can also call it multiple times without causing inconsistencies in the dialog structure.
 *
 * 
 * @section DIALOG_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Dialog
 *
 * Dialogs have only one fundamental callback method, namely the \ref DIALOGEXEC method.
 * This method has to be implemented for every dialog; the other callback methods are optional.
 * In the C++ wrapper class ObjDialog, the scip_exec() method (which corresponds to the \ref DIALOGEXEC callback) is a virtual
 * abstract member function.
 * You have to implement it in order to be able to construct an object of your dialog class.
 * 
 * Additional documentation to the callback methods can be found in "type_dialog.h".
 *
 * @subsection DIALOGEXEC
 *
 * The DIALOGEXEC method is invoked, if the user selected the dialog's command name in the parent's menu. It should 
 * execute what is stated in DIALOG_DESC, e.g., the display constraint handlers dialog should display information about 
 * the constraint handlers included in SCIP, see "src/scip/dialog_default.c". 
 *
 * For typical methods called by the execution method, have a look at "src/scip/dialog_default.c".
 * 
 * The callback has to return which dialog should be processed next. This can be, for example, the root dialog 
 * (SCIPdialoghdlrGetRoot()), the parent dialog (SCIPdialogGetParent()) or NULL, which stands for closing the interactive 
 * shell.
 *
 *
 * @section DIALOG_ADDITIONALCALLBACKS Additional Callback Methods of a Dialog
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to free private data.
 *
 * @subsection DIALOGPFREE
 *
 * If you are using dialog data, you have to implement this method in order to free the dialog data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_DIALOGFREE(dialogFreeMydialog)
 * {
 *    SCIP_DIALOGDATA* dialogdata;
 *  
 *    dialogdata = SCIPdialogGetData(dialog);
 *    assert(dialogdata != NULL);
 *
 *    SCIPfreeMemory(scip, &dialogdata);
 *
 *    SCIPdialogSetData(dialog, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection DIALOGDESC
 *
 * The method is called, when the help menu of the parent is displayed. It should output (usually a single line of) 
 * information describing the meaning of the dialog. 
 * \n
 * If this callback is not implemented, the description string of the dialog (DIALOG_DESC) is displayed instead.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DISP How to add display columns
 *
 * While solving a constraint integer program, SCIP displays status information in a column-like fashion. The current 
 * number of processed branching tree nodes, the solving time, and the relative gap between primal and dual bound are 
 * examples of such display columns. There already exists a wide variety of display columns which can be activated or 
 * deactivated on demand, see "src/scip/disp_default.c". Additionally, the user can implement his own display columns
 * in order to track problem or algorithm specific values.  
 * \n
 * A complete list of all displays contained in this release can be found \ref DISPLAYS "here".
 *
 * In the following, we explain how the user can add an own display column. 
 * We give the explanation for creating an own source file for each additional display column. Of course, you can collect 
 * different additional display columns in one source file.
 * Take "src/scip/disp_default.c", where all default display columns are collected, as an example.
 * As all other default plugins, the default display column plugins and the display column template are written in C.
 * C++ users can easily adapt the code by using the ObjDisp wrapper base class and implement the scip_...() virtual methods
 * instead of the SCIP_DECL_DISP... callback methods.
 * 
 *
 * Additional documentation for the callback methods of a display column can be found in the file "type_disp.h".
 *
 * Here is what you have to do to implement a display column:
 * -# Copy the template files "src/scip/disp_xxx.c" and "src/scip/disp_xxx.h" into files named "disp_mydisplaycolumn.c"
 *    and "disp_mydisplaycolumn.h".
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xxx" by "mydisplaycolumn".
 * -# Adjust the properties of the display column (see \ref DISP_PROPERTIES).
 * -# Define the display column data (see \ref DISP_DATA).
 * -# Implement the interface methods (see \ref DISP_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref DISP_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref DISP_ADDITIONALCALLBACKS).
 *
 *
 * @section DISP_PROPERTIES Properties of a Display Column
 *
 * At the top of the new file "disp_mydisplaycolumn.c" you can find the display column properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the display column properties by calling the constructor
 * of the abstract base class ObjDisp from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par DISP_NAME: the name of the display column.
 * This name is used in the interactive shell to address the display column.
 * Additionally, if you are searching for a display column with SCIPfindDisp(), this name is looked up.
 * Names have to be unique: no two display columns may have the same name.
 *
 * \par DISP_DESC: the description of the display column.
 * This string is printed as description of the display column in the interactive shell.
 *
 * \par DISP_HEADER: the header of the display column.
 * This string is printed as header of the display column in the status information display.
 *
 * \par DISP_WIDTH: the width of the display column.
 * This parameter defines the width (number of characters) of the display column. The value of the parameter has to be
 * greater than or equal to the number of characters in the header string. 
 *
 * \par DISP_PRIORITY: the priority of the display column.
 * The total width of status information lines is bounded by the parameter "display width". The display columns actually contained
 * in the status information display are selected in decreasing order of their priority. Furthermore, the user can force 
 * columns to be displayed or not to be displayed in the status information display. For that, he has to switch the value 
 * of the display column's parameter "active" from "auto" (its default value) to "on" or "off", respectively. 
 *
 * \par DISP_POSITION: the relative position of the display column.
 * In the status information display, the display columns are arranged from left to right in increasing order of their 
 * relative position.  
 *
 * \par DISP_STRIPLINE: the default for whether the display column should be separated with a line from its right neighbor.
 * This parameter states whether the display column should be separated with the string "|" from its right neighbor. In so 
 * doing, the clearness of the status information display may improve.  
 * 
 * @section DISP_DATA Display Column Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_DispData".
 * In this data structure, you can store the data of your display column. For example, you should store the adjustable 
 * parameters of the display column in this data structure.
 * If you are using C++, you can add display column data as usual as object variables to your class.
 * \n
 * Defining display column data is optional. You can leave the struct empty.
 *
 *
 * @section DISP_INTERFACE Interface Methods
 *
 * At the bottom of "disp_mydisplaycolumn.c" you can find the interface method SCIPincludeDispMydisplaycolumn(), which also
 * appears in "disp_mydisplaycolumn.h".
 * \n
 * This method has only to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the display column by calling the method
 * SCIPincludeDisp().
 *
 * The interface method is called by the user, if he wants to include the display column, i.e., if he wants to use the display column in his 
 * application. 
 * Note that additional display column plugins have to be included <b>before the default display columns plugin</b>, i.e.,
 * the SCIPincludeDispMydisplaycolumn() call has to occure before the SCIPincludeDispDefault() call.  
 *
 * If you are using display column data, you have to allocate the memory for the data at this point.
 * You can do this by calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &dispdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_DispData afterwards.
 *
 * Although this is very uncommon, you may also add user parameters for your display column, see the method
 * SCIPincludeConshdlrKnapsack() in the knapsack constraint handler src/scip/cons_knapsack.c for an example.
 *
 * 
 * @section DISP_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Display Column
 *
 * Display column plugins have only one fundamental callback method, namely the DISPOUTPUT method.
 * This method has to be implemented for every display column; the other callback methods are optional.
 * In the C++ wrapper class ObjDisp, the scip_output() method (which corresponds to the DISPOUTPUT callback) is a virtual
 * abstract member function.
 * You have to implement it in order to be able to construct an object of your display column class.
 *
 * Additional documentation to the callback methods can be found in "type_disp.h".
 *
 * @subsection DISPOUTPUT
 *
 * The DISPOUTPUT callback is called after each pricing loop during node processing and after a node has been processed. 
 * In addition, at the root node, the callback is executed after each iteration of the price-and-cut loop. 
 * It should write the display column information for the current node to a given output file stream. 
 *
 * Typical methods called by a display column are, for example, SCIPdispLongint(), SCIPdispInt(), SCIPdispTime(), and 
 * SCIPinfoMessage().
 *
 *
 * @section DISP_ADDITIONALCALLBACKS Additional Callback Methods of a Display Column
 *
 * The additional callback methods need not to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 *
 * @subsection DISPFREE
 *
 * If you are using display column data, you have to implement this method in order to free the display column data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_DISPFREE(dispFreeMydisplaycolumn)
 * {
 *    SCIP_DISPDATA* dispdata;
 *  
 *    dispdata = SCIPdispGetData(disp);
 *    assert(dispdata != NULL);
 *
 *    SCIPfreeMemory(scip, &dispdata);
 *
 *    SCIPdispSetData(disp, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection DISPINIT
 *
 * The DISPINIT callback is executed after the problem was transformed.
 * The display column may, e.g., use this call to initialize its display column data.
 *
 * @subsection DISPEXIT
 *
 * The DISPEXIT callback is executed before the transformed problem is freed.
 * In this method, the display column should free all resources that have been allocated for the solving process in 
 * DISPINIT.
 *
 * @subsection DISPINITSOL
 *
 * The DISPINITSOL callback is executed when the presolving was finished and the branch and bound process is about to 
 * begin. The display column may use this call to initialize its branch and bound specific data.
 *
 * @subsection DISPEXITSOL
 *
 * The DISPEXITSOL callback is executed before the branch and bound process is freed. The display column should use this 
 * call to clean up its branch and bound data specific data.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page OBJ Creating, capturing, releasing, and adding data objects
 *
 *  Data objects (variables, constraints, rows, ... ) are subject to reference counting
 *  to avoid expensive copying operations. This concept is similar to smart pointers.
 *  Creating such an object (e.g., by calling SCIPcreateVar()) will set the
 *  reference counter to one. Capturing an object (e.g., by calling SCIPcaptureVar()) increases the reference counter,
 *  releasing it (e.g., by calling SCIPreleaseVar()) decreases the counter. If the reference counter gets zero, the
 *  object will be destroyed automatically.
 *
 *  Remember that a created data object is automatically captured. If the user
 *  doesn't need the object anymore, he has to call the object's release method.
 *
 *  When a data object is added to SCIP (e.g., by calling SCIPaddVar()) , it is captured again, such that a
 *  release call does not destroy the object. If SCIP doesn't need the object
 *  anymore, it is automatically released.
 *  
 *  E.g., if the user calls
 * \code
 *  SCIPcreateVar(); // reference counter 1
 *  SCIPaddVar(); // reference counter 2
 *  SCIPreleaseVar(); // reference counter 1
 * \endcode
 *  the reference counter will be 1 afterwards, and the variable will be destroyed, if SCIP frees the problem.
 *  If the user wants to use this variable, e.g. for extracting statistics after SCIP was finished, the user must not call
 *  SCIPreleaseVar() right after adding the variable, but before terminating the program.
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

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DEBUG Debugging
 *
 *  If you want to debug your own code that uses SCIP, there are some tricks that we collect as
 *  follows - check whether one of these helps in your case.
 * 
 *  - Use the debug mode (<code>make OPT=dbg</code>, see \ref MAKE) and run the code.
 *  - Use asserts in your code (see \ref CODE).
 *  - Turn on additional debug output by placing <code>\#define SCIP_DEBUG</code> at the top of SCIP files you
 *    want to analyze. This will output messages included in the code using <code>SCIPdebugMessage()</code> (see \ref EXAMPLE_1).
 *    We recommend to also use <code>SCIPdebugMessage()</code> in your own code for being able to activate 
 *    debug output in the same way. 
 *  - If available on your system, we recommend to use a debugger like gdb to trace all function calls on the stack,
 *    display values of certain expressions, manually break the running code, and so forth.
 *  - If available on your system, you can use software like valgrind to check for uninitialized
 *    values or segmentation faults.
 *  - For checking the usage of SCIP memory, you can use
 *    <code>SCIPprintMemoryDiagnostic()</code>. This outputs memory that is currently in use. This is
 *    almost always only useful after a <code>SCIPfree()</code> call.
 *  - If your code cuts off a feasible solution, but you do not know which component is responsible,
 *    you define <code>SCIP_DEBUG_SOLUTION</code> in the file <code>debug.h</code> to be a filename
 *    containing a solution in SCIP format (see \ref EXAMPLE_2). 
 *    This solution is then read and it is checked for every cut, whether the solution violates the cut.
 * 
 *  @section EXAMPLE_1 How to activate debug messages
 *    For example, if we include a <code>\#define SCIP_DEBUG</code> at the top of heur_oneopt.c and recompile in DBG mode, 
 *    and run the scip interactive shell to solve p0033.mps from the miplib, we get some output like:
 * \code
 * SCIP version 1.1.0 [precision: 8 byte] [memory: block] [mode: debug] [LP solver: SoPlex 1.4.0]
 * Copyright (c) 2002-2008 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
 * 
 * user parameter file <scip.set> not found - using default parameters
 * 
 * SCIP> read check/IP/miplib/p0033.mps
 * original problem has 33 variables (33 bin, 0 int, 0 impl, 0 cont) and 16 constraints
 * SCIP> optimize
 * ...
 *  0.1s|     1 |     0 |   132 | 257k|   0 |  14 |  30 |  13 |  13 |  30 |  51 |  39 |   0 |   0 | 3.026472e+03 | 3.347000e+03 |  10.59%
 * [src/scip/heur_oneopt.c:332] debug: Row <R122> has activity 110
 * [src/scip/heur_oneopt.c:332] debug: Row <R123> has activity 216
 * ...
 * [src/scip/heur_oneopt.c:101] debug: Try to shift down variable <t_C157> with
 * [src/scip/heur_oneopt.c:102] debug:     lb:<-0> <= val:<1> <= ub:<1> and obj:<171> by at most: <1>
 * [src/scip/heur_oneopt.c:135] debug:  -> The shift value had to be reduced to <0>, because of row <R122>.
 * [src/scip/heur_oneopt.c:137] debug:     lhs:<-1e+20> <= act:<110> <= rhs:<148>, colval:<-60>
 * ...
 * [src/scip/heur_oneopt.c:383] debug:  Only one shiftcand found, var <t_C167>, which is now shifted by<-1.0>
 * k 0.1s|     1 |     0 |   132 | 258k|   0 |  14 |  30 |  13 |  13 |  30 |  51 |  39 |   0 |   0 | 3.026472e+03 | 3.164000e+03 |   4.54%
 * [src/scip/heur_oneopt.c:436] debug: found feasible shifted solution:
 * objective value:                     3164.00000000012
 * C157                                                1   (obj:171)
 * C163                                                1   (obj:163)
 * C164                                                1   (obj:69)
 * C170                                                1   (obj:49)
 * C172                                                1   (obj:258)
 * C174                                                1   (obj:250)
 * C175                                                1   (obj:500)
 * C179                                                1   (obj:318)
 * C181                                                1   (obj:318)
 * C182                                                1   (obj:159)
 * C183                                 1.00000000000038   (obj:318)
 * C184                                                1   (obj:159)
 * C185                                                1   (obj:318)
 * C186                                                1   (obj:114)
 * [src/scip/heur_oneopt.c:498] debug: Finished 1-opt heuristic
 * ...
 * \endcode 
 *
 *  @section EXAMPLE_2 How to add a debug solution
 * Continuing the example above, we finish the solving process.
 * The optimal solution can now be written to a file:
 * \code
 * SCIP> display solution
 * 
 * objective value:                                 3089
 * C157                                                1   (obj:171)
 * C163                                                1   (obj:163)
 * C164                                                1   (obj:69)
 * C166                                                1   (obj:183)
 * C170                                                1   (obj:49)
 * C174                                                1   (obj:250)
 * C177                                                1   (obj:500)
 * C179                                                1   (obj:318)
 * C181                                                1   (obj:318)
 * C182                                                1   (obj:159)
 * C183                                                1   (obj:318)
 * C184                                                1   (obj:159)
 * C185                                                1   (obj:318)
 * C186                                                1   (obj:114)
 * 
 * SCIP> write solution check/p0033.sol
 * 
 * written solution information to file <check/p0033.sol>
 * \endcode
 * 
 * If we afterwards comment in 
 * <code>\#define SCIP_DEBUG_SOLUTION "check/p0033.sol"</code> in debug.h, recompile and run SCIP,
 * it will output 
 * \code
 * SCIP> read check/IP/miplib/p0033.mps
 * original problem has 33 variables (33 bin, 0 int, 0 impl, 0 cont) and 16 constraints
 * SCIP> optimize
 * 
 * presolving:
 * ***** debug: reading solution file <check/p0033.sol>
 * ***** debug: read 15 non-zero entries
 * \endcode
 * Further debug output would only appear, if the solution was cut off in the solving process.
 * Of course, this is not the case! Hopefully...otherwise, please send a bug report ;-)
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page TEST How to run automated tests with SCIP
 *
 *  SCIP comes along with a set of useful tools that allow us to perform automated tests. The
 *  following is a step-by-step guide from setting up the test environment to evaluation and
 *  customization of test runs.
 *
 *
 *  @section SETUP Setting up the test environment
 *
 *  Suppose we have a test set of one or more "lp" files. These test instances have to accessible in
 *  the @c check directory; you can also set a link to the instances. Once this is done, we have
 *  to create a @em test file within the @c check directory that specifies the location of our test
 *  instances. In our example, we have three "lp" files, <tt>one.lp</tt>, <tt>two.lp</tt>, and
 *  <tt>three.lp</tt>, stored in the directory <code>check/testrun</code>, a sample @em test file
 *  <tt>testrun.test</tt> would look like this
 *
 *  - <tt>testrun/one.lp</tt>
 *  - <tt>testrun/two.lp</tt>
 *  - <tt>testrun/three.lp</tt>
 *
 *  Optionally, we can provide a "solu" file in the "check" directory containing the best known
 *  objective values for our test instances. SCIP can use these values to verify the results. The
 *  file has to have the same basename as the @em test file, i.e. in our case
 *  <tt>testrun.solu</tt>. The file content may look as follows
 * 
 *  - <tt>=opt=    one.lp     73389.09</tt>
 *  - <tt>=opt=    two.lp     153183.89</tt>
 *  - <tt>=best=  three.lp   234234.21</tt>
 *
 *  The string <tt>=best=</tt> indicates that the objective value may not be optimal.
 *
 *
 *  @section STARTING Starting a test run
 *
 *  We are now in a position to start the test run. The easiest way to do so, is by calling
 *
 *  <code>make TEST=testrun test</code>
 *
 *  in the SCIP root directory. Note that <tt>testrun</tt> is exactly the basename of our @em test
 *  file. This will cause SCIP to solve our test instances one after another and to create various
 *  output files (see @ref EVAL "Evaluating a test run for details").
 *
 * 
 *  @section EVAL Evaluating a test run
 *
 *  During computation, SCIP automatically creates the following output files in the <code>check/result</code>
 *  directory.
 *
 *  \arg <tt>*.out</tt> - output of <tt>stdout</tt>
 *  \arg <tt>*.err</tt> - output of <tt>stderr</tt>
 *  \arg <tt>*.set</tt> - copy of the settings file used
 *
 *  In order to obtain a summary of the computational results, we call the <tt>evalcheck.sh</tt>
 *  script in the @c check directory with the "out" file as argument. This produces the following
 *  files in the <code>check/result</code> directory"
 *
 *  \arg <tt>*.res</tt> - ASCII table containing a summary of the computational results
 *  \arg <tt>*.tex</tt> - TeX table containing a summary of the computational results
 *  \arg <tt>*.pav</tt> - <a href="http://www.gamsworld.org/performance/paver">PAVER</a> output
 *
 *  @b Note The basename of each output file allows us to reconstruct the test setting. Suppose, for
 *  example, the output file is
 *
 *      <tt>check.testrun.scip_binary.machine.fast.out</tt>
 *
 *  Then <tt>testrun</tt> is the basename of the @em test file and <tt>scip_binary</tt> is the name of
 *  the SCIP binary file. Furthermore, <tt>machine</tt> denotes the machine on which the test was
 *  run and <tt>fast</tt> is the basename of the settings file used.
 *
 *
 *  @section USING Using customized setting files
 *
 *  It is possible to use customized settings files for the test run. These have to be placed in the
 *  @c settings directory. Note that the access to settings files in subfolders of the @c settings
 *  directory is currently not possible.
 *
 *  To run SCIP with a custom settings file, say <tt>fast.set</tt>, we call
 *
 *      <tt>make TEST=testrun SETTING=fast test</tt>
 *
 *  in the SCIP root directory.
 *
 * 
 *  @section ADVANCED Advanced options
 *
 *  We can further customize the test run by specifying the following options in the <tt>make</tt>
 *  call.
 *
 *  \arg <tt>TIME</tt>  - time limit for each test instance in seconds [default: 3600]
 *  \arg <tt>NODES</tt> - node limit [default: 2100000000]
 *  \arg <tt>MEM</tt>   -  memory limit in MB [default: 1536]
 *  \arg <tt>CONTINUE</tt> - continue test run [default: "false"]
 *
 * 
 *  @section COMPARE Comparing test runs for different settings
 *
 *  Often test runs are performed on the basis of different settings. In this case, it is useful to
 *  have a performance comparison. For this purpose, we can use the <tt>allcmpres.sh</tt> script in
 *  the @c check directory.
 *
 *  Suppose, we performed our test run with two different settings, say <tt>fast.set</tt> and
 *  <tt>slow.set</tt>. Assuming that all other parameters (including the SCIP binary), were the same,
 *  we may have the following "res" files in the <code>check/result</code> directory
 *
 *  <tt>check.testrun.scip_binary.machine.fast.res</tt> and @n
 *  <tt>check.testrun.scip_binary.machine.slow.res</tt>
 *
 *  For a comparison of both computations we simply call
 *
 *  <code>allcmpres.sh results/check.testrun.scip_binary.machine.fast.res results/check.testrun.scip_binary.machine.slow.res</code>
 *
 *  in the @c check directory. This produces two ASCII tables on the console that provide a detailed
 *  performance comparison of both test runs. Note that the first "res" file serves as the reference
 *  computation. In the following list explains the output, where we use the term "solver" for the
 *  combination of SCIP with a specific settings file.
 *
 *  \arg <tt>Nodes</tt> - Number of nodes processed.
 *  \arg <tt>Time</tt>  - Computation time in seconds.
 *  \arg <tt>F</tt>     - If no feasible solutions were found, then '#', empty otherwise.
 *  \arg <tt>NodQ</tt>  - Equals Nodes(i) / Nodes(0), where 'i' denotes the current solver and '0' stands for the reference solver.
 *  \arg <tt>TimQ</tt>  - Equals Time(i) / Time(0).
 *  \arg <tt>bounds check</tt> - Status of the primal and dual bound check.
 *
 *  \arg <tt>proc</tt> - Number of instances processed.
 *  \arg <tt>eval</tt> - Number of instances evaluated (bounds check = "ok"). Only these instances are used in the calculation of the mean values.
 *  \arg <tt>fail</tt> - Number of instances with bounds check = "fail".
 *  \arg <tt>time</tt> - Number of instances with timeout.
 *  \arg <tt>solv</tt> - Number of instances correctly solved within the time limit.
 *  \arg <tt>wins</tt> - Number of instances on which the solver won (i.e., the
 *      solver was at most 10% slower than the fastest solver OR had the best
 * 	primal bound in case the instance was not solved by any solver within
 *	the time limit).
 *  \arg <tt>bett</tt>    - Number of instances on which the solver was better than the
 *	reference solver (i.e. more than 10% faster).
 *  \arg <tt>wors</tt>    - Number of instances on which the solver was worse than the
 *	reference solver (i.e. more than 10% slower).
 *  \arg <tt>bobj</tt>    - Number of instances on which the solver had a better primal
 *	bound than the reference solver (i.e. a difference larger than 10%).
 *  \arg <tt>wobj</tt>    - Number of instances on which the solver had a worse primal
 *	bound than the reference solver (i.e. a difference larger than 10%).
 *  \arg <tt>feas</tt>    - Number of instances for which a feasible solution was found.
 *  \arg <tt>nodes</tt>   - Geometric mean of the processed nodes over all evaluated instances.
 *  \arg <tt>shnodes</tt> - Shifted geometric mean of the processed nodes over all evaluated instances.
 *  \arg <tt>nodesQ</tt>  - Equals nodes(i) / nodes(0), where 'i' denotes the current
 *	solver and '0' stands for the reference solver.
 *  \arg <tt>shnodesQ</tt> - Equals shnodes(i) / shnodes(0).
 *  \arg <tt>time</tt>    - Geometric mean of the computation time over all evaluated instances.
 *  \arg <tt>shtime</tt>  - Shifted geometric mean of the computation time over all evaluated instances.
 *  \arg <tt>timeQ</tt>   - Equals time(i) / time(0).
 *  \arg <tt>shtimeQ</tt> - Equals shtime(i) / shtime(0).
 *  \arg <tt>score</tt>   - N/A
 *
 *  \arg <tt>all</tt>   - All solvers.
 *  \arg <tt>diff</tt>  - Solvers with instances that differ from the reference solver in the number of processed nodes or number of iterations.
 *  \arg <tt>equal</tt> - Solvers with instances whose number of processed nodes and number of
 *       iteration is equal to the reference solver (including a 10% tolerance) and where no timeout
 *       occured.
 *  \arg <tt>optimal auto settings</tt> - Theoretical result for a solver that performed 'best of all' for every instance.
 */


/**@page FAQ Frequently Asked Questions (FAQ)
 * \htmlinclude faqcss.inc  
 * \htmlinclude faq.inc  
 */


/**@defgroup BRANCHINGRULES Branching Rules 
 * @brief In the following you find a list of all branching rule which are currently available.
 */

/**@defgroup CONSHDLRS  Constraint Handler 
 * @brief In the following you find a list of all constraint handlers which are currently available.
 */

/**@defgroup DIALOGS Dialogs 
 * @brief In the following you find a list of all dialogs which are currently available.
 */

/**@defgroup DISPLAYS Displays 
 * @brief In the following you find a list of all displays (output columns)  which are currently available.
 */

/**@defgroup FILEREADERS File Readers 
 * @brief In the following you find a list of all file readers which are currently available.
 */
 
/**@defgroup LPIS LP Interfaces
 * @brief In the following you find a list of all LP instances which are currently available.
 */

/**@defgroup NODESELECTORS Node Selectors
 * @brief In the following you find a list of all node selectors which are currently available.
 */

/**@defgroup PRESOLVERS Presolvers
 * @brief In the following you find a list of all presolvers which are currently available.
 */

/**@defgroup PRICERS Pricers
 * @brief In the following you find a list of all pricers which are currently available.
 */

/**@defgroup PRIMALHEURISTICS Primal Heuristics
 * @brief In the following you find a list of all primal heuristics which are currently available.
 */

/**@defgroup PROPAGATORS Propagators
 * @brief In the following you find a list of all propagators which are currently available.
 */

/**@defgroup RELAXATORS Relaxators
 * @brief In the following you find a list of all relaxators which are currently available.
 */

/**@defgroup SEPARATORS Separators
 * @brief In the following you find a list of all separators  which are currently available.
 */

