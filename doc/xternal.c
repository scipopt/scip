/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  this file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*                  2002-2013 Konrad-Zuse-Zentrum                            */
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
 * @author Gerald Gamrath
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Mathias Kinder
 * @author Marc Pfetsch
 * @author Stefan Vigerske
 * @author Robert Waniek
 * @author Kati Wolter
 * @author Michael Winkler
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Overview (\OTHERDOCU)
 *
 * \OTHERDOCUTEXT
 *
 *
 *
 * @section WHATISSCIP What is SCIP?
 *
 * SCIP is a framework to solve constraint integer programs (CIPs). In particular,
 *
 * - SCIP is a branch-and-cut-and-price framework,
 * - incorporates a full-scale mixed integer programming (MIP) solver, and
 * - incorporates a full-scale mixed integer quadratically constrained programming (MIQCP) solver.
 *
 * See the web site of <a href="http://scip.zib.de">SCIP</a> for more information about licensing and to download SCIP.
 *
 * SCIP is developed together with <a href="http://www3.mathematik.tu-darmstadt.de/ags/optimierung/research/discrete-optimization.html">TU Darmstadt</a> and
 * <a href="http://www.am.uni-erlangen.de/wima/">University of Erlangen-N&uuml;rnberg (Chair of EDOM)</a>
 * and has more than 500,000 lines of C code.
 *
 * @section GETTINGSTARTED Getting started
 *
 * - \ref MAKE    "Installation information / Makefiles"
 * - \ref LICENSE "License"
 *
 * - \ref SHELL       "Tutorial: the interactive shell"
 * - \ref FILEREADERS "Readable file formats"
 * - \ref START       "How to start a new project"
 * - \ref EXAMPLES    "Examples"
 *
 * @section FURTHERINFORMATION References
 *
 * - \ref PUBLICMETHODS "List of callable functions"
 * - \ref PARAMETERS "List of all SCIP parameters"
 *
 * - \ref DOC     "How to search the documentation for interface methods"
 * - \ref FAQ     "Frequently asked questions (FAQ)"
 *
 *
 * @section PROGRAMMING Programming with SCIP
 *
 * @subsection CODINGBASICS Coding basics for SCIP
 *
 *   - \ref CODE    "Coding style guidelines"
 *   - \ref OBJ     "Creating, capturing, releasing, and adding data objects"
 *   - \ref DEBUG   "Debugging"
 *
 * @subsection HOWTOADD How to add ...
 *   - \ref CONS    "Constraint handlers"
 *   - \ref PRICER  "Variable pricers"
 *   - \ref PRESOL  "Presolvers"
 *   - \ref SEPA    "Separators"
 *   - \ref PROP    "Propagators"
 *   - \ref BRANCH  "Branching rules"
 *   - \ref NODESEL "Node selectors"
 *   - \ref HEUR    "Primal heuristics"
 *   - \ref RELAX   "Relaxation handlers"
 *   - \ref READER  "File readers"
 *   - \ref DIALOG  "Dialogs"
 *   - \ref DISP    "Display columns"
 *   - \ref EVENT   "Event handler"
 *   - \ref NLPI    "Interfaces to NLP solvers"
 *   - \ref EXPRINT "Interfaces to expression interpreters"
 *   - \ref CONF    "Conflict analysis"
 *   - \ref PARAM   "additional user parameters"
 *
 * - \ref TEST     "How to run automated tests with SCIP"
 * - \ref COUNTER  "How to use SCIP to count feasible solutions"
 *
 *
 * @section FURTHERINFO Further information
 *
 * @subsection CHG Changes between different versions of SCIP
 * - \ref CHANGELOG    "Change log"
 * - \ref RELEASENOTES "Release notes"
 * - \ref CHG6         "Interface changes between version 2.1 and 3.0"
 * - \ref CHG5         "Interface changes between version 2.0 and 2.1"
 * - \ref CHG4         "Interface changes between version 1.2 and 2.0"
 * - \ref CHG3         "Interface changes between version 1.1 and 1.2"
 * - \ref CHG2         "Interface changes between version 1.0 and 1.1"
 * - \ref CHG1         "Interface changes between version 0.9 and 1.0"
 *
 * @subsection AUTHORS SCIP Authors
 * - <a class="el" href="AUTHORS.shtml#main">Current main developers</a>
 * - <a class="el" href="AUTHORS.shtml#further">Further developers</a>
 * - <a class="el" href="AUTHORS.shtml#contributors">Contributors</a>
 *
 * @version  3.0.2
 *
 * \image html scippy.png
 *
 */

/** @page EXAMPLES Examples projects
 *
 *  SCIP contains several examples that demonstrate its usage. They are contained in the &quot;examples&quot; directory
 *  in the source code distribution.
 *
 *  @section BRANCHANDPRICE Branch-and-price
 *
 *  <table>
 *  <tr>
 *  <td>
 *  <a href="http://scip.zib.de/doc/examples/Binpacking/index.shtml"><b>Binpacking</b></a>
 *  </td>
 *  <td>
 *  An implementation of the column generation approach for the binpacking problem. It includes a customized reader,
 *  Ryan/Foster branching rule, (global) problem data, variable data, and constraint handler.
 *  </td>
 *  </tr>
 *  <tr>
 *  <td>
 *  <a href="http://scip.zib.de/doc/examples/Coloring/index.shtml"><b>Coloring</b></a>
 *  </td>
 *  <td>
 *  An implemenation of the column generation approach for graph coloring of Mehrotra and Trick.
 *  </td>
 *  </tr>
 *  <tr>
 *  <td>
 *  <a href="http://scip.zib.de/doc/examples/VRP/index.shtml"><b>VRP</b></a>
 *  </td>
 *  <td>
 *  A solver for a simple capacity-constrained vehicle routing problem, which is based on pricing tours via a dynamic
 *  programming algorithm.
 *  </td>
 *  </tr>
 *  </table>
 *
 *  @section BRANCHANDCUT Branch-and-cut
 *
 *  <table>
 *  <tr>
 *  <td>
 *  <a href="http://scip.zib.de/doc/examples/LOP/index.shtml"><b>LOP</b></a>
 *  </td>
 *  <td>
 *  An example for implementing a constraint handler.
 *  </td>
 *  </tr>
 *  <tr>
 *  <td>
 *  <a href="http://scip.zib.de/doc/examples/TSP/index.shtml"><b>TSP</b></a>
 *  </td>
 *  <td>
 *  A short implementations of a constraint handler, two easy combinatorial heuristics, a file reader, etc. which
 *  demonstrate the usage of SCIP as a branch-and-cut-framework for solving traveling salesman problem instances.
 *  </td>
 *  </tr>
 *  </table>
 *
 *  @section CALLABLELIBRARY Callable library
 *
 *  <table>
 *  <tr>
 *  <td>
 *  <a href="http://scip.zib.de/doc/examples/CallableLibrary/index.shtml"><b>CallableLibrary</b></a>
 *  </td>
 *  <td>
 *  An example showing how to setup constraints (esp. nonlinear ones) when using SCIP as callable library.
 *  </td>
 *  </tr>
 *  <tr>
 *  <td>
 *  <a href="http://scip.zib.de/doc/examples/MIPSolver/index.shtml"><b>MIPSolver</b></a>
 *  </td>
 *  <td>
 *  A minimal implementation for using SCIP included into another source code
 *  </td>
 *  </tr>
 *  <tr>
 *  <td>
 *  <b>Queen</b>
 *  </td>
 *  <td>
 *  An example showing the use of SCIP as callable library.
 *  </td>
 *  </tr>
 *  </table>
 *
 *
 *  @section OTHERPLUGINS Other plugins
 *
 *  <table>
 *  <tr>
 *  <td>
 *  <a href="http://scip.zib.de/doc/examples/Eventhdlr/index.shtml"><b>Eventhdlr</b></a>
 *  </td>
 *  <td>
 *  A small example illustrating the use of an event handler.
 *  </td>
 *  </tr>
 *  <tr>
 *  <td>
 *  <a href="http://scip.zib.de/doc/examples/Scheduler/index.shtml"><b>Scheduler</b></a>
 *  </td>
 *  <td>
 *  An example containing three readers and one primal heuristic for scheduling problems.
 *  </td>
 *  </tr>
 *  </table>
 *
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CODE Coding style guidelines
 *
 * We follow the following coding style guidelines and recommend them for all developers.
 *
 * - Indentation is 3 spaces. No tabs anywhere in the code.
 * - Always only one declaration in a line.
 * - Braces are on a new line and not indented.
 * - Spaces around all operators.
 * - No spaces between control structure keywords like "if", "for", "while", "switch" and the corresponding brackets.
 * - No spaces between a function name and the parenthesis in both the definition and function calls.
 * - Use assert() to show preconditions for the parameters, invariants and postconditions.
 * - All global functions start with "SCIP". In the usual naming scheme this is followed by the object and a method name
 *   like in SCIPlpAddRow(). Functions return TRUE or FALSE should be named like SCIPisFeasEQ().
 * - Make all functions that are not used outside the module 'static'. Naming should start with a lower case letter.
 * - Variable names should be all lower case.
 * - For each structure there is a typedef with the name in all upper case.
 * - Defines should be named all upper case.
 * - Document functions, parameters, and variables in a doxygen conformed way.
 *
 * As an example, have a look at tree.c and see the examples below. We also provide settings for
 * \ref XEMACS "(x)emacs" and \ref ECLIPSE "eclipse".
 *
 * @section CODEEXAMPLES Examples
 *
 * In this section we state a few examples illustrating the \SCIP code style.
 *
 * \code
 * #ifdef __cplusplus
 * extern "C" {
 * #endif
 *
 * /** SCIP operation stage */
 * enum SCIP_Stage
 * {
 *    SCIP_STAGE_INIT         =  0,        /**< SCIP datastructures are initialized, no problem exists */
 *    SCIP_STAGE_PROBLEM      =  1,        /**< the problem is being created and modified */
 *    SCIP_STAGE_TRANSFORMING =  2,        /**< the problem is being transformed into solving data space */
 *    SCIP_STAGE_TRANSFORMED  =  3,        /**< the problem was transformed into solving data space */
 *    SCIP_STAGE_PRESOLVING   =  4,        /**< the problem is being presolved */
 *    SCIP_STAGE_PRESOLVED    =  5,        /**< the problem was presolved */
 *    SCIP_STAGE_INITSOLVE    =  6,        /**< the solving process data is being initialized */
 *    SCIP_STAGE_SOLVING      =  7,        /**< the problem is being solved */
 *    SCIP_STAGE_SOLVED       =  8,        /**< the problem was solved */
 *    SCIP_STAGE_FREESOLVE    =  9,        /**< the solving process data is being freed */
 *    SCIP_STAGE_FREETRANS    = 10         /**< the transformed problem is being freed */
 * };
 * typedef enum SCIP_Stage SCIP_STAGE;
 *
 * /** possible settings for enabling/disabling algorithms and other features */
 * enum SCIP_Setting
 * {
 *    SCIP_UNDEFINED = 0,                  /**< undefined setting */
 *    SCIP_DISABLED  = 1,                  /**< feature is disabled */
 *    SCIP_AUTO      = 2,                  /**< feature is set to automatic mode */
 *    SCIP_ENABLED   = 3                   /**< feature is enabled */
 * };
 * typedef enum SCIP_Setting SCIP_SETTING;
 *
 * #ifdef __cplusplus
 * }
 * #endif
 * \endcode
 *
 * @section XEMACS Customize (x)emacs
 *
 * If you are using (x)emacs, you can use the following customization for the c++-mode. These settings satisfy the
 * coding guidelines of \SCIP.
 *
 * \verbatim
  (add-hook 'c++-mode-hook
    (function
      (lambda ()
    ;; SCIP customizations for c-mode and c++-mode
    (setq-default c-basic-offset 3)
    (c-set-offset 'substatement-open 0)
    (c-set-offset 'statement-case-open 0)
    (c-set-offset 'brace-list-open '-)
    (c-set-offset 'inextern-lang '0)
    (c-set-offset 'arglist-intro '+)
    (c-set-offset 'arglist-cont 0)
    (c-set-offset 'arglist-cont-nonempty '+)
    (c-set-offset 'arglist-close '+)
    (set-variable 'fill-column 120)
   ;; this will make sure spaces are used instead of tabs
    (setq tab-width 8 indent-tabs-mode nil)
    )))\endverbatim
 *
 * @section ECLIPSE Customize eclipse
 *
 *
 * Eclipse user can use the profile below. This profile does not match the \SCIP coding guideline completely.
 *
 * \code
 *
 * <?xml version="1.0" encoding="UTF-8" standalone="no"?>
 * <profiles version="1">
 * <profile kind="CodeFormatterProfile" name="scip" version="1">
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_paren_in_method_declaration" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_paren_in_for" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_new_line_in_empty_block" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.lineSplit" value="124"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_comma_in_base_types" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.keep_else_statement_on_same_line" value="false"/>
 * <setting id="org.eclipse.cdt.core.formatter.indent_switchstatements_compare_to_switch" value="false"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_brace_in_array_initializer" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_comma_in_method_declaration_parameters" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_paren_in_if" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_paren_in_exception_specification" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_paren_in_parenthesized_expression" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_comma_in_base_types" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.indent_body_declarations_compare_to_access_specifier" value="true"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_paren_in_exception_specification" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_comma_in_template_arguments" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_brace_in_block" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_paren_in_method_declaration" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.use_tabs_only_for_leading_indentations" value="false"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_colon_in_labeled_statement" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_colon_in_case" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_comma_in_array_initializer" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_comma_in_enum_declarations" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.alignment_for_expressions_in_array_initializer" value="16"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_comma_in_declarator_list" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_bracket" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_paren_in_for" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_prefix_operator" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.tabulation.size" value="3"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_new_line_before_else_in_if_statement" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.alignment_for_enumerator_list" value="48"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_paren_in_parenthesized_expression" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_between_empty_parens_in_method_declaration" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.alignment_for_declarator_list" value="16"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_paren_in_switch" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_paren_in_parenthesized_expression" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.indent_empty_lines" value="false"/>
 * <setting id="org.eclipse.cdt.core.formatter.indent_switchstatements_compare_to_cases" value="true"/>
 * <setting id="org.eclipse.cdt.core.formatter.keep_empty_array_initializer_on_one_line" value="false"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_brace_in_method_declaration" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.put_empty_statement_on_new_line" value="true"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_brace_in_switch" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_paren_in_cast" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_between_empty_braces_in_array_initializer" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.brace_position_for_method_declaration" value="next_line"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_paren_in_while" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_question_in_conditional" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_semicolon" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_closing_angle_bracket_in_template_arguments" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_colon_in_base_clause" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.indent_breaks_compare_to_cases" value="true"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_unary_operator" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_comma_in_declarator_list" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.alignment_for_arguments_in_method_invocation" value="16"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_paren_in_while" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_between_empty_brackets" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_bracket" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.alignment_for_parameters_in_method_declaration" value="48"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_new_line_before_closing_brace_in_array_initializer" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.number_of_empty_lines_to_preserve" value="1"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_paren_in_method_invocation" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_brace_in_array_initializer" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_semicolon_in_for" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.brace_position_for_block" value="next_line"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_colon_in_conditional" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.brace_position_for_type_declaration" value="next_line"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_assignment_operator" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_angle_bracket_in_template_arguments" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_comma_in_expression_list" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_angle_bracket_in_template_parameters" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.continuation_indentation" value="1"/>
 * <setting id="org.eclipse.cdt.core.formatter.alignment_for_expression_list" value="0"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_paren_in_method_declaration" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_comma_in_template_parameters" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_colon_in_default" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_binary_operator" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.alignment_for_conditional_expression" value="16"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_between_empty_parens_in_method_invocation" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_comma_in_array_initializer" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_paren_in_if" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.format_guardian_clause_on_one_line" value="false"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_paren_in_cast" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.indent_access_specifier_compare_to_type_header" value="false"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_brace_in_type_declaration" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_colon_in_labeled_statement" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.continuation_indentation_for_array_initializer" value="1"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_comma_in_method_declaration_parameters" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_semicolon_in_for" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_paren_in_method_invocation" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.indent_body_declarations_compare_to_namespace_header" value="false"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_closing_brace_in_block" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_assignment_operator" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.alignment_for_compact_if" value="0"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_brace_in_array_initializer" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_new_line_at_end_of_file_if_missing" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_comma_in_template_parameters" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_comma_in_expression_list" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_question_in_conditional" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_paren_in_exception_specification" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_binary_operator" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_new_line_before_identifier_in_function_declaration" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.alignment_for_base_clause_in_type_declaration" value="80"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_comma_in_method_declaration_throws" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_between_empty_parens_in_exception_specification" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_comma_in_method_invocation_arguments" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.indent_declaration_compare_to_template_header" value="false"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_unary_operator" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_paren_in_switch" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.indent_statements_compare_to_body" value="true"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_comma_in_method_declaration_throws" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.indent_statements_compare_to_block" value="true"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_comma_in_template_arguments" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_new_line_before_catch_in_try_statement" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.alignment_for_throws_clause_in_method_declaration" value="48"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_paren_in_method_invocation" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_closing_paren_in_cast" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_paren_in_catch" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_angle_bracket_in_template_parameters" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.tabulation.char" value="space"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_angle_bracket_in_template_parameters" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_paren_in_while" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_comma_in_method_invocation_arguments" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.brace_position_for_block_in_case" value="next_line"/>
 * <setting id="org.eclipse.cdt.core.formatter.compact_else_if" value="true"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_postfix_operator" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_colon_in_base_clause" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_new_line_after_template_declaration" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_paren_in_catch" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.keep_then_statement_on_same_line" value="false"/>
 * <setting id="org.eclipse.cdt.core.formatter.brace_position_for_switch" value="next_line"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_paren_in_if" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_paren_in_switch" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.keep_imple_if_on_one_line" value="false"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_new_line_after_opening_brace_in_array_initializer" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.indentation.size" value="3"/>
 * <setting id="org.eclipse.cdt.core.formatter.brace_position_for_namespace_declaration" value="end_of_line"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_colon_in_conditional" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_comma_in_enum_declarations" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_prefix_operator" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_angle_bracket_in_template_arguments" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.brace_position_for_array_initializer" value="end_of_line"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_colon_in_case" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_paren_in_catch" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_brace_in_namespace_declaration" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_postfix_operator" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_closing_bracket" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_new_line_before_while_in_do_statement" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_before_opening_paren_in_for" value="do not insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_closing_angle_bracket_in_template_parameters" value="insert"/>
 * <setting id="org.eclipse.cdt.core.formatter.insert_space_after_opening_angle_bracket_in_template_arguments" value="do not insert"/>
 * </profile>
 * </profiles>
 * \endcode
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page MAKE Makefiles / Installation information
 *
 *
 * In most cases (LINUX and MAC) it is quite easy to compile and install SCIP. Therefore, reading the section
 * \ref BRIEFINSTALL "Brief installation description" should usually be enough. If this is not the case you find
 * \ref DETAILEDINSTALL "Detailed installation description" below as well as \ref EXAMPLE1 "Examples".
*
 * @section BRIEFINSTALL Brief installation description
 *
 * The easiest way to install SCIP is to use the SCIP Optimization Suite which contains SCIP, SoPlex, and ZIMPL. For
 * that we refer to the INSTALL file of the SCIP Optimization Suite (main advantage: there is no need
 * to specify any directories, the compiling process is fully automated).
 *
 * Compiling SCIP directly can be done as follows:
 *
 * -# unpack the tarball <code>tar xvf scip-x.y.z.tgz</code>
 * -# change to the directory <code>cd scip-x.y.z</code>
 * -# start compiling SCIP by typing <code>make</code>
 * -# (optional) install the header, libraries, and binary <code>make install INSTALLDIR="/usr/local/</code>
 *
 * During your first compilation you will be asked for some soft-link targets,
 * depending on the LP solver you want to use. Usually, SCIP needs the
 * following information
 * -# the directory where the include files of the LP solver lie
 * -# the library file(s) "lib*.a" or/and "lib*.so"
 *
 * Besides that, SCIP needs some soft-link targets, for ZIMPL
 * -# the directory where the include files of ZIMPL lie
 * -# the library file(s) "lib*.a" or/and "lib*.so"
 *
 * You will need either the .a or the .so files and can skip the others by
 * just pressing return.
 *
 * The most common compiling issue is that some libraries are missing
 * on your system or that they are outdated. SCIP per default requires
 * zlib, gmp and readline.  Try compiling with: <code> make ZLIB=false
 * READLINE=false ZIMPL=false</code> or, better, install them. Note
 * that under Linux-based systems, you need to install the
 * developer-versions of gmp/zlib/readline, in order to also have the
 * header-files available.
 *
 @section DETAILEDINSTALL Detailed
 * installation description
 *
 * In this section we describe the use, and a few features, of the SCIP Makefile. We also give two examples for how to install
 * SCIP. The \ref EXAMPLE1 "first example" illustrates the default installation. This means, with SoPleX and ZIMPL. The
 * \ref EXAMPLE2 "second example" shows how to get CPLEX linked to SCIP without ZIMPL. This is followed by a section which
 * gives some hints on what to do if the \ref COMPILERPROBLEMS "compilation throws an error". We give some comments on
 * how to install SCIP under \ref WINDOWS "WINDOWS" and show \ref RUN "how to start SCIP".
 *
 * If you experience any problems during the installation, you will find help in the \ref INSTALL "INSTALL" file.
 *
 * SCIP contains a makefile system, which allows the individual setting of several parameters. For
 * instance, the following settings are supported:
 *
 * - <code>OPT=\<dbg|opt|opt-gccold\></code> Here <code>dbg</code> turns on the debug mode of SCIP. This enables asserts
 *   and avoids macros for several function in order to ease debugging. The default is <code>opt</code>, which enables
 *   the optimized mode. The third option <code>opt-gccold</code> will work with older GCC compilers before version
 *   4.2. We recommend using newer GCC versions.
 *
 * - <code>LPS=\<clp|cpx|grb|msk|qso|spx|xprs|none\></code> This determines the LP-solver, which should have been
 *   installed separately from SCIP. The options are the following:
 *      - <code>clp</code>: COIN-OR Clp LP-solver
 *      - <code>cpx</code>: CPLEX LP-solver
 *      - <code>grb</code>: Gurobi LP-solver (interface is in beta stage)
 *      - <code>msk</code>: Mosek LP-solver
 *      - <code>qso</code>: QSopt LP-solver
 *      - <code>spx</code>: SoPlex LP-solver (default)
 *      - <code>xprs</code>: XPress LP-solver
 *      - <code>none</code>: no LP-solver (you should set the parameter \<lp/solvefreq\> to \<-1\> to avoid solving LPs)
 *
 * - <code>LPSOPT=\<dbg|opt|opt-gccold\></code> Chooses the debug or optimized version (or old GCC optimized) version of
 *   the LP-solver. (currently only available for SoPlex and CLP)
 *
 * - <code>ZIMPL=\<true|false\></code> Turns direct support of ZIMPL in SCIP on (default) or off, respectively.
 * - <code>ZIMPLOPT=\<dbg|opt|opt-gccold\></code> Chooses the debug or optimized (default) (or old GCC optimized)
 *   version of ZIMPL, if ZIMPL support is enabled. \n
 *   If the ZIMPL-support is disabled, the GMP-library is no longer needed for SCIP and therefore not linked to SCIP.
 *
 * - <code>READLINE=\<true|false\></code> Turns support via the readline library on (default) or off, respectively.
 *
 * - <code>IPOPT=\<true|false\></code> to enable/disable(default) IPOPT interface (needs IPOPT)
 *
 * - <code>EXPRINT=\<cppad|none\></code>   to use CppAD as expressions interpreter or no expressions interpreter (default)
 *
 * There are additional parameters for Linux/Gnu compilers:
 *
 * - <code>OPT=noblkmem</code> turns off the internal SCIP memory.  This way the code can be checked by valgrind or
 *   similar tools.
 * - <code>OPT=opt-shared</code> generates a shared object of the SCIP libraries.  (The binary uses these shared
 *   libraries as well.)
 * - <code>OPT=prf</code> generates a profiling version of SCIP providing a detailed statistic of the time usage of
 *   every method of SCIP.
 *
 * You can use other compilers - depending on the system:
 *
 * - <code>COMP=intel</code> Uses of the Intel compiler which is only available with the main optimization flags
 *   <code>OPT=\<dbg|opt\></code>. (Default is gcc/g++ represented through <code>COMP=gnu</code>.)
 *
 * There is the possibility to watch the compilation more precisely:
 *
 * - <code>VERBOSE=\<true|false\></code> Turns the extensive output on or off (default).
 *
 * The SCIP makefile supports several targets (used via <code>make ... "target"</code>):
 *
 * - <code>links</code> Reconfigures the links in the "lib" directory.
 * - <code>doc</code> Creates documentation in the "doc" directory.
 * - <code>clean</code> Removes all object files.
 * - <code>depend</code> Creates dependencies files. This is only needed if you add files to SCIP.
 * - <code>check</code> Runs the check script, see \ref TEST.
 *
 * The SCIP makefiles are structured as follows.
 *
 * - <code>Makefile</code> This is the basic makefile in the SCIP root directory. It loads
 *   additional makefile information depending on the parameters set.
 * - <code>make/make.project</code> This file contains definitions that are useful for all codes
 *   that use SCIP, for instance, the examples.
 * - <code>make.\<sys\>.\<machine\>.\<compiler\>.\<dbg|opt|prf|opt-gccold\></code> These file contain system/compiler specific
 *   definitions. If you have an unsupported compiler, you can copy one of these and modify it
 *   accordingly.
 *
 * If your platform or compiler is not supported by SCIP you might try and copy one of the existing
 * makefiles in the <code>make</code> directory and modify it. If you succeed, we are always
 * interested in including more Makefiles into the system.
 *
 *
 * @section EXAMPLE1 Example 1 (defaults: SoPlex, with ZIMPL support):
 *
 * Typing <code>make</code> uses SoPlex as LP solver and includes support for the modeling language ZIMPL. You will be asked the
 * following questions on the first call to "make" (example answers are already given):
 *
 * \verbatim
  > make
  make[1]: Entering directory `scip-1.2'

  - Current settings: LPS=spx OSTYPE=linux ARCH=x86_64 COMP=gnu SUFFIX= ZIMPL=true ZIMPLOPT=opt IPOPT=false IPOPTOPT=opt

  * SCIP needs some softlinks to external programs, in particular, LP-solvers.
  * Please insert the paths to the corresponding directories/libraries below.
  * The links will be installed in the 'lib' directory.
  * For more information and if you experience problems see the INSTALL file.

    -> "spxinc" is the path to the SoPlex "src" directory, e.g., "../../soplex/src".
    -> "libsoplex.*" is the path to the SoPlex library, e.g., "../../soplex/lib/libsoplex.linux.x86.gnu.opt.a"
    -> "zimplinc" is a directory containing the path to the ZIMPL "src" directory, e.g., "../../zimpl/src".
    -> "libzimpl.*" is the path to the ZIMPL library, e.g., "../../zimpl/lib/libzimpl.linux.x86.gnu.opt.a"

  - preparing missing soft-link "lib/spxinc":
  > Enter soft-link target file or directory for "lib/spxinc" (return if not needed):
  > ../../soplex/src/
  -> creating softlink "lib/spxinc" -> "../../soplex/src"


  - preparing missing soft-link "lib/libsoplex.linux.x86_64.gnu.opt.a":
  > Enter soft-link target file or directory for "lib/libsoplex.linux.x86_64.gnu.opt.a" (return if not needed):
  > ../../soplex/lib/libsoplex.linux.x86_64.gnu.opt.a
  -> creating softlink "lib/libsoplex.linux.x86_64.gnu.opt.a" -> "../../soplex/lib/libsoplex.linux.x86_64.gnu.opt.a"


  - preparing missing soft-link "lib/libsoplex.linux.x86_64.gnu.opt.so":
  * this soft-link is not necessarily needed since "lib/libsoplex.linux.x86_64.gnu.opt.a" already exists - press return to skip
  > Enter soft-link target file or directory for "lib/libsoplex.linux.x86_64.gnu.opt.so" (return if not needed):
  >
  * skipped creation of softlink "lib/libsoplex.linux.x86_64.gnu.opt.so". Call "make links" if needed later.


  - preparing missing soft-link "lib/zimplinc/zimpl":
  > Enter soft-link target file or directory for "lib/zimplinc/zimpl" (return if not needed):
  ../../zimpl/src/
   creating softlink "lib/zimplinc/zimpl" -> "../../zimpl/src"


  - preparing missing soft-link "lib/libzimpl.linux.x86_64.gnu.opt.a":
  > Enter soft-link target file or directory for "lib/libzimpl.linux.x86_64.gnu.opt.a" (return if not needed):
  > ../../zimpl/lib/libzimpl.linux.x86_64.gnu.opt.a
  -> creating softlink "lib/libzimpl.linux.x86_64.gnu.opt.a" -> "../../zimpl/lib/libzimpl.linux.x86_64.gnu.opt.a"


  - preparing missing soft-link "lib/libzimpl.linux.x86_64.gnu.opt.so":
  * this soft-link is not necessarily needed since "lib/libzimpl.linux.x86_64.gnu.opt.a" already exists - press return to skip
  > Enter soft-link target file or directory for "lib/libzimpl.linux.x86_64.gnu.opt.so" (return if not needed):
  >
  * skipped creation of softlink "lib/libzimpl.linux.x86_64.gnu.opt.so". Call "make links" if needed later.

  ...

  -> generating library lib/libobjscip-1.2.0.linux.x86_64.gnu.opt.a
  -> generating library lib/liblpispx-1.2.0.linux.x86_64.gnu.opt.a
  -> generating library lib/libscip-1.2.0.linux.x86_64.gnu.opt.a
  -> linking bin/scip-1.2.0.linux.x86_64.gnu.opt.spx

 * \endverbatim
 *
 * @section EXAMPLE2 Example 2 (CPLEX, with no ZIMPL support):
 *
 * Typing <code>make LPS=cpx ZIMPL=false</code>  uses CPLEX as LP solver. You will be asked the following questions on
 * the first call to "make" (example answers are already given):
 *
 * \verbatim
  > make LPS=cpx ZIMPL=false
  make[1]: Entering directory `scip-1.2'

  - Current settings: LPS=cpx OSTYPE=linux ARCH=x86_64 COMP=gnu SUFFIX= ZIMPL=false ZIMPLOPT=opt IPOPT=false IPOPTOPT=opt

  * SCIP needs some softlinks to external programs, in particular, LP-solvers.
  * Please insert the paths to the corresponding directories/libraries below.
  * The links will be installed in the 'lib' directory.
  * For more information and if you experience problems see the INSTALL file.

    -> "cpxinc" is the path to the CPLEX "include" directory, e.g., "<CPLEX-path>/include/ilcplex".
    -> "libcplex.*" is the path to the CPLEX library, e.g., "<CPLEX-path>/lib/x86_rhel4.0_3.4/static_pic/libcplex.a"

  - preparing missing soft-link "lib/cpxinc":
  > Enter soft-link target file or directory for "lib/cpxinc" (return if not needed):
  > ../../cplex121/include
  -> creating softlink "lib/cpxinc" -> "../../cplex121/include"


  - preparing missing soft-link "lib/libcplex.linux.x86_64.gnu.a":
  > Enter soft-link target file or directory for "lib/libcplex.linux.x86_64.gnu.a" (return if not needed):
  > ../../cplex121/lib/x86-64_sles9.0_3.3/static_pic/libcplex.a
  -> creating softlink "lib/libcplex.linux.x86_64.gnu.a" -> "../../../../adm_cple/cplex121/lib/x86-64_sles9.0_3.3/static_pic/libcplex.a"


  - preparing missing soft-link "lib/libcplex.linux.x86_64.gnu.so":
  > Enter soft-link target file or directory for "lib/libcplex.linux.x86_64.gnu.so" (return if not needed):
  >
  * skipped creation of softlink "lib/libcplex.linux.x86_64.gnu.so". Call "make links" if needed later.

  ...

  -> generating library lib/libobjscip-1.2.0.linux.x86_64.gnu.opt.a
  -> generating library lib/liblpicpx-1.2.0.linux.x86_64.gnu.opt.a
  -> generating library lib/libscip-1.2.0.linux.x86_64.gnu.opt.a
  -> linking bin/scip-1.2.0.linux.x86_64.gnu.opt.cpx

 * \endverbatim
 *
 * @section COMPILERPROBLEMS Compilation problems:
 *
 * - If the soft-link query script does not work on your machine, read step 2 in the \ref INSTALL "INSTALL" file for
 * instructions on manually creating the soft-links.
 *
 * - If you get an error message of the type\n
 * <code>make: *** No rule to make target `lib/???', needed by `obj/O.linux.x86.gnu.opt/lib/scip/???.o'.  Stop.</code>\n
 * the corresponding soft-link was not created or points to a wrong location.  Check the soft-link targets in the "lib/"
 * subdirectory. Try to delete all soft-links from the "lib/" directory\n and call "make links" to generate them
 * again. If this still fails, read step 2 for instructions on manually\n creating the soft-links.
 *
 * - If you get an error message of the type\n
 * <code>make: *** No rule to make target `make/make.?.?.?.?.?'.  Stop.</code>,\n
 * the corresponding machine dependent makefile for your architecture and compiler is missing.\n Create one of the given
 * name in the "make/" subdirectory. You may take\n "make/make.linux.x86.gnu.opt" or any other file in the make
 * subdirectory as example.\n
 *
 * - The readline library seems to differ slightly on different OS distributions. Some versions do
 * not support the <code>remove_history()</code> call.  In this case, you have to either add
 * <code>-DNO_REMOVE_HISTORY</code> to the FLAGS in the appropriate "make/make.*" file, or to
 * compile with <code>make USRFLAGS=-DNO_REMOVE_HISTORY</code>.  Make sure, the file
 * "src/scip/dialog.c" is recompiled.  If this doesn't work either, disable the readline library
 * with <code>make READLINE=false</code>.
 *
 * - On some systems, the <code>sigaction()</code> method is not available. In this case, you have
 * to either add <code>-DNO_SIGACTION</code> to the FLAGS in the appropriate "make/make.*" file, or
 * to compile with <code>make USRFLAGS=-DNO_SIGACTION</code>.  Make sure, the file
 * "src/scip/interrupt.c" is recompiled.
 *
 * - On some systems, the <code>rand_r()</code> method is not available.  In this case, you have to either add
 * <code>-DNO_RAND_R</code> to the FLAGS in the appropriate "make/make.*" file, or to compile with
 * <code>make USRFLAGS=-DNO_RAND_R</code>.  Make sure, the file "src/scip/misc.c" is recompiled.
 *
 * - On some systems, the <code>strtok_r()</code> method is not available.  In this case, you have
 * to either add <code>-DNO_STRTOK_R</code> to the FLAGS in the appropriate make/make.* file, or to
 * compile with <code>make USRFLAGS=-DNO_STRTOK_R</code>.  Make sure, the file "src/scip/misc.c" is
 * recompiled.
 *
 * - On some systems, the <code>strerror_r()</code> method is not available.  In this case, you have
 * to either add <code>-DNO_STRERROR_R</code> to the FLAGS in the appropriate "make/make.*" file, or
 * to compile with <code>make USRFLAGS=-DNO_STRERROR_R</code>.  Make sure, the file
 * "src/scip/misc.c" is recompiled.
 *
 * - On some systems, the option [-e] is not available for the read command.  You have to compile with READ=read.
 *
 * - If you encounter other compiler or linker errors, you should recompile with <code>make
 * VERBOSE=true ...</code> in order to get the full compiler invocation. This might help to fix the
 * corresponding machine dependent makefile in the make subdirectory.
 *
 * @section WINDOWS Remarks on Installing under Windows using MinGW
 *
 * To build your own windows binaries under windows we recommend using the MinGW-Compiler with MSYS
 * from <a href="http://www.mingw.org">www.mingw.org</a> .
 *
 * First install MSYS, then MinGW to the mingw folder inside the msys folder.
 * Now you need to install the following packages to the mingw folder:
 * - zlib (or use ZLIB=false)
 * - pcre (here suffices the pcre7.0-lib.zip (or equivalent) to be extracted into the mingw-folder)
 *
 * After calling <code>make clean</code> in the ZIMPL folder you will also need flex and bison to
 * remake ZIMPL. We recommend NOT to use <code>"make clean"</code> inside the ZIMPL-folder if you do
 * not have these packages installed.
 *
 * You can download these additional packages from <a href="http://gnuwin32.sourceforge.net/packages.html">here</a>
 * or compile the source on your own from their homepages.
 *
 * Second you need to copy the file <code>sh.exe</code> to <code>bash.exe</code> otherwise various
 * scripts (including makefiles) will not work.  Normally <code>unistd.h</code> covers also the
 * getopt-options, but for mingw you need to add the entry <code>\#include <getopt.h></code> into
 * "/mingw/include/unistd.h" after the other include-entries (if not present).
 *
 * Finally, there is one package you need to compile if you want to use ZIMPL and ZIMPL-support in
 * SCIP (otherwise use <code>ZIMPL=false</code> as parameter with the make-call): the
 * <code>gmplib</code> from <a href="http://www.gmplib.org">gmplib.org</a>. The command
 * <code>./configure --prefix=/mingw ; make ; make install</code> should succeed without problems
 * and installs the gmplib to the mingw folder.
 *
 * Now <code>make READLINE=false</code> should be compiling without errors.  Please note that we
 * do NOT support creating the doxygen documentation and readline-usage under windows.
 *
 *
 * @section RUN How to run SCIP after successfully compiling SCIP
 *
 * To run the program, enter <code>bin/scip</code> for the last compiled version. If you have more than one compiled
 * binary (i. e., one in debug and one in optimized mode) and wish to specify the binary, type
 * <code>bin/scip.\$(OSTYPE).\$(ARCH).\$(COMP).\$(OPT).\$(LPS)</code>
 * (e.g. <code>bin/scip.linux.x86_64.gnu.opt.spx</code>).
 *
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page START How to start a new project
 *
 * Once you succeeded installing SCIP together with an LP-solver on your system,
 * you have a powerful tool for solving MIPs, MIQCPs,
 * MINLPs, etc... at hand. SCIP can also be customized to the type of problems you
 * are working on by additional plugins.
 * Instructions on how to write a new plugin and include it in SCIP can be found in the corresponding
 * \ref HOWTOADD "How to add ... pages".
 *
 * SCIP can also be used for writing your own branch-and-cut or branch-and-cut-and-price code. SCIP already
 * provides a number of existing code examples which we suggest as both reference and starting point
 * for these kinds of projects.
 * Below, you find some hints of how to start such a project.
 *
 * The example should be chosen
 *     depending on the programming language (<b>C</b> or <b>C++</b>) and the purpose
 *   (<b>branch-and-cut</b> or <b>branch-and-cut-and-price</b>) of your project.
 *
 *    We suggest the use one of the following examples:
 *     - The <a href="http://scip.zib.de/doc/examples/VRP/index.shtml"><b>VRP</b></a>-example is a <b>branch-and-cut-and-price</b> (column generation)-code
 *       in <b>C++</b>.
 *     - The <a href="http://scip.zib.de/doc/examples/Coloring/index.shtml"><b>Coloring</b></a>
 *        and the <a href="http://scip.zib.de/doc/examples/Binpacking/index.shtml"><b>Binpacking</b></a>-example are
 *       <b>branch-and-cut-and-price</b> (column generation)-codes in <b>C</b>.
 *     - The <a href="http://scip.zib.de/doc/examples/TSP/index.shtml"><b>TSP</b></a>-example
 *        is a <b>branch-and-cut</b>-code in <b>C++</b>.
 *     - The <a href="http://scip.zib.de/doc/examples/LOP/index.shtml"><b>LOP</b></a>-example
 *         is a <b>branch-and-cut</b>-code in <b>C</b>.

 * - Copy one of the examples in the <code>examples</code> directory (in the SCIP root
 *   directory). For instance, type
 *   \verbatim
 > cp -r examples/Coloring/ ../SCIPProject/ ; cd ../SCIPProject
 *   \endverbatim
 *
 *   from the SCIP root directory for copying the content of the <code>Coloring</code>-example into a fresh
 *   directory named SCIPProject in the parent directory of the SCIP root directory and jumping to
 *   the new SCIPProject directory rightafter.
 *
 *  - Open the <code>Makefile</code>  via
 *    \verbatim
 > kate Makefile
     \endverbatim
 *
 *    and edit the following variables at the top to have a compilable code:
 *
 *    - specify a correct path to the SCIP root (<code>SCIPDIR</code>)
 *    - rename the targets name (<code>MAINNAME</code>)
 *    - adjust the source file names (<code>MAINOBJ</code>).
 *
 * - Once you have edited the makefile, you can use all the flags that can be used in SCIP to
 *   compile your code, see \ref MAKE.
 *   Note that you need to update the dependency files before compiling your project via <code>make depend</code>.
 *
 *
 *
 *
 */


/**@page SHELL Tutorial: the interactive shell
 *
 * If are using SCIP as a black box solver, here you will find some tips and tricks what you can do.
 *
 * First of all, we need a SCIP binary and an example problem file to work with.  Therefore, you can either download the
 * SCIP standard distribution (which includes problem files) and compile it on your own or you can download a
 * precompiled binary and an example problem separately. SCIP can read files in LP, MPS, ZPL, WBO, FZN, PIP, OSiL, and other formats (see \ref FILEREADERS).
 *
 * If you want to download the source code of the SCIP standard distribution, we recommend to go to the <a
 * href="http://scip.zib.de/download.shtml">SCIP download section</a>, download the latest release (version 3.0 as
 * of this writing), inflate the tarball (e.g., with "tar xzf scipoptsuite-[version].tgz"), and follow the instructions
 * in the INSTALL file. The instance stein27, which will serve as an example in this tutorial, can be found under
 * scipoptsuite-[version]/scip-[version]/check/instances/MIP/stein27.mps.
 *
 * If you want to download a precompiled binary, go to the <a href="http://scip.zib.de/download.shtml">SCIP download
 * section</a> and download an appropriate binary for your operating system. To follow this tutorial, we recommend downloading the instance
 * <a href="http://miplib.zib.de/miplib3/miplib3/stein27.mps.gz">stein27</a> from
 * the <a href="http://miplib.zib.de/miplib3/miplib.html">MIPLIB 3.0</a> homepage.
 *
 * Now start your binary, without any arguments. This opens the interactive shell, which should look somehow like this:
 *
 * \code
 * SCIP version 2.0.1 [precision: 8 byte] [memory: block] [mode: optimized] [LP solver: SoPlex 1.5.0]
 * Copyright (c) 2002-2013 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
 *
 * External codes:
 *   SoPlex 1.5.0         Linear Programming Solver developed at Zuse Institute Berlin (soplex.zib.de)
 *   ZIMPL 3.1.0          Zuse Institute Mathematical Programming Language developed by T. Koch (zimpl.zib.de)
 *
 * user parameter file <scip.set> not found - using default parameters
 *
 * SCIP>
 * \endcode
 *
 * First of all "help" shows you a list of all available shell commands. Brackets indicate a submenu with further options.
 * \code
 * SCIP> help

 *  <display>             display information
 *  <set>                 load/save/change parameters
 * ...
 *  read                  read a problem
 * \endcode
 *
 * Okay, let's solve some MIPs... use "read <path/to/file>" to parse a problem file, "optimize" to solve it and "display
 * solution" to show the nonzero variables of the best found solution.

 * \code
 * SCIP> read check/instances/MIP/stein27.mps
 * original problem has 27 variables (27 bin, 0 int, 0 impl, 0 cont) and 118 constraints
 * SCIP> optimize
 *
 * feasible solution found by trivial heuristic, objective value  2.700000e+01
 * presolving:
 * (round 1) 0 del vars, 0 del conss, 0 chg bounds, 0 chg sides, 0 chg coeffs, 118 upgd conss, 0 impls, 0 clqs
 * presolving (2 rounds):
 *  0 deleted vars, 0 deleted constraints, 0 tightened bounds, 0 added holes, 0 changed sides, 0 changed coefficients
 *  0 implications, 0 cliques
 * presolved problem has 27 variables (27 bin, 0 int, 0 impl, 0 cont) and 118 constraints
 *       1 constraints of type <knapsack>
 *     117 constraints of type <logicor>
 * transformed objective value is always integral (scale: 1)
 * Presolving Time: 0.00
 *
 *  time | node  | left  |LP iter|LP it/n| mem |mdpt |frac |vars |cons |cols |rows |cuts |confs|strbr|  dualbound   | primalbound  |  gap
 * t 0.0s|     1 |     0 |    34 |     - | 337k|   0 |  21 |  27 | 118 |  27 | 118 |   0 |   0 |   0 | 1.300000e+01 | 2.700000e+01 | 107.69%
 * R 0.0s|     1 |     0 |    34 |     - | 338k|   0 |  21 |  27 | 118 |  27 | 118 |   0 |   0 |   0 | 1.300000e+01 | 2.600000e+01 | 100.00%
 * s 0.0s|     1 |     0 |    34 |     - | 339k|   0 |  21 |  27 | 118 |  27 | 118 |   0 |   0 |   0 | 1.300000e+01 | 2.500000e+01 |  92.31%
 *   0.0s|     1 |     0 |    44 |     - | 392k|   0 |  21 |  27 | 118 |  27 | 120 |   2 |   0 |   0 | 1.300000e+01 | 2.500000e+01 |  92.31%
 * b 0.0s|     1 |     0 |    44 |     - | 393k|   0 |  21 |  27 | 118 |  27 | 120 |   2 |   0 |   0 | 1.300000e+01 | 1.900000e+01 |  46.15%
 * ...
 *   0.1s|     1 |     2 |   107 |     - | 920k|   0 |  24 |  27 | 118 |  27 | 131 |  13 |   0 |  24 | 1.300000e+01 | 1.900000e+01 |  46.15%
 * R 0.1s|    14 |    10 |   203 |   7.4 | 935k|  13 |   - |  27 | 118 |  27 | 124 |  13 |   0 | 164 | 1.300000e+01 | 1.800000e+01 |  38.46%
 *   0.1s|   100 |    54 |   688 |   5.9 | 994k|  13 |  20 |  27 | 118 |  27 | 124 |  13 |   0 | 206 | 1.300000e+01 | 1.800000e+01 |  38.46%
 *   0.1s|   200 |    86 |  1195 |   5.5 |1012k|  13 |   - |  27 | 119 |  27 | 124 |  13 |   1 | 207 | 1.300000e+01 | 1.800000e+01 |  38.46%
 *  time | node  | left  |LP iter|LP it/n| mem |mdpt |frac |vars |cons |cols |rows |cuts |confs|strbr|  dualbound   | primalbound  |  gap
 *   0.2s|   300 |   106 |  1686 |   5.3 |1024k|  13 |   - |  27 | 119 |  27 | 124 |  13 |   1 | 207 | 1.350000e+01 | 1.800000e+01 |  33.33%
 * ...
 *   0.7s|  4100 |    50 | 18328 |   4.4 |1033k|  16 |   8 |  27 | 119 |  27 | 124 |  13 |  15 | 207 | 1.650000e+01 | 1.800000e+01 |   9.09%
 *
 * SCIP Status        : problem is solved [optimal solution found]
 * Solving Time (sec) : 0.73
 * Solving Nodes      : 4192
 * Primal Bound       : +1.80000000000000e+01 (283 solutions)
 * Dual Bound         : +1.80000000000000e+01
 * Gap                : 0.00 %
 *
 * SCIP> display solution
 *
 * objective value:                                   18
 * x0001                                               1   (obj:1)
 * x0003                                               1   (obj:1)
 * ...
 * x0027                                               1   (obj:1)
 *
 * SCIP>
 * \endcode
 *
 * What do we see here? After "optimize", SCIP first goes into presolving. Not much is happening for this instance, just
 * the linear constraints get upgraded to more specific types. Each round of presolving will be displayed in a single
 * line, with a short summary at the end. Here, there has only been one round with actual changes, the second round did
 * not bring any further reductions.  Thus, it is not displayed and presolving is stopped. Then, we see the actual
 * solving process. The first three output lines indicate that new incumbent solutions were found by the primal
 * heuristics with display characters "t", "R", and "s"; see, how the "primalbound" column goes down from 27 to 25. In
 * the fourth line, two "cuts" are added.  Up to here, we needed 44 "LP iter"ations (34 for the first LP and 10 more to
 * resolve after adding cuts). Little later, the root node processing is finished. We see that there are now two open
 * nodes in the "left" column. From now on, we will see an output line every hundredth node or whenever a new incumbent
 * is found (e.g. at node 14 in the above output). After some more nodes, the "dualbound" starts moving, too. At one
 * point, both will be the same, and the solving process terminates, showing us some wrap-up information.
 *
 * The exact performance varies amongst different architectures, operating systems, and so on. Do not be worried if
 * your installation needs more or less time or nodes to solve. Also, this instance has more than 2000 different optimal
 * solutions. The optimal objective value always has to be 18, but the solution vector may differ. If you are interested
 * in this behavior, which is called "performance variability", you may have a look at the MIPLIB2010 paper.
 *
 * We might want to have some more information now. Which were the heuristics that found the solutions? What plugins
 *  were called during the solutions process and how much time did they spend? How did the instance that we were solving
 *  look?  Information on certain plugin types (e.g., heuristics, branching rules, separators) we get by
 *  "display <plugin-type>", information on the solution process, we get by "display statistics", and "display problem"
 *  shows us the current instance.
 *
  \code
 * SCIP> display heuristics
 *  primal heuristic     c priority freq ofs  description
 *  ----------------     - -------- ---- ---  -----------
 *  trivial              t    10000    0   0  start heuristic which tries some trivial solutions
 * ...
 *  rounding             R    -1000    1   0  LP rounding heuristic with infeasibility recovering
 *  shifting             s    -5000   10   0  LP rounding heuristic with infeasibility recovering also using continuous variables
 * ...
 * SCIP> display statistics
 * ...
 *   gomory           :       0.02          6          0          0        461          0
 *   cgmip            :       0.00          0          0          0          0          0
 *   strongcg         :       0.01          6          0          0        598          0
 * ...
 *   oneopt           :       0.01          4          1
 *   coefdiving       :       0.02         57          0
 * ...
 *   primal LP        :       0.00          0          0       0.00          -
 *   dual LP          :       0.20       4187      14351       3.43   71755.00
 * ...
 * \endcode
 *
 * We see that rounding and shifting were the heuristics producing the solutions in the beginning. Rounding is called at
 * every node, shifting only at every tenth level of the tree. The statistics are quite comprehensive, thus, we just
 * explain a few lines here. We get information for all types of plugins and for the overall solving process. Besides
 * others, we see that in six calls, the gomory cut separator and the strong Chv&aacute;tal-Gomory separator each produced
 * several hundred cuts (of which only a few entered the LP). The oneopt heuristic found one solution in 4 calls,
 * whereas coefdiving failed all 57 times it was called. All the LPs have been solved with the dual simplex algorithm, which
 * took about 0.2 seconds of the 0.7 seconds overall solving time.
 *
 * Now, we can start playing around with parameters. Rounding and shifting seem to be quite successful on this instance,
 * wondering what happens if we disable them? Or what happens, if we are even more rigorous and disable all heuristics?
 * Or if we do the opposite and use aggressive heuristics?
 *
 * \code
 * SCIP> set
 *
 *   <branching>           change parameters for branching rules
 *  ...
 *   <heuristics>          change parameters for primal heuristics
 *
 * SCIP/set> heuristics
 *
 *   <actconsdiving>       LP diving heuristic that chooses fixings w.r.t. the active constraints
 *  ...
 *   <shifting>            LP rounding heuristic with infeasibility recovering also using continuous variables
 *  ...
 *
 * SCIP/set/heuristics> shifting
 *
 *   <advanced>            advanced parameters
 *   freq                  frequency for calling primal heuristic <shifting> (-1: never, 0: only at depth freqofs) [10]
 *   freqofs               frequency offset for calling primal heuristic <shifting> [0]
 *
 * SCIP/set/heuristics/shifting> freq
 * current value: 10, new value [-1,2147483647]: -1
 * heuristics/shifting/freq = -1
 *
 * SCIP> se he rou freq -1
 * heuristics/rounding/freq = -1
 *
 * SCIP> re check/instances/MIP/stein27.mps
 * original problem has 27 variables (27 bin, 0 int, 0 impl, 0 cont) and 118 constraints
 * SCIP> o
 *
 * feasible solution found by trivial heuristic, objective value  2.700000e+01
 * ...
 * z 0.1s|     3 |     4 |   140 |  10.5 |1060k|   2 |  22 |  27 | 118 |  27 | 123 |  14 |   0 |  66 | 1.300000e+01 | 1.900000e+01 |  46.15%
 * z 0.1s|     6 |     7 |   176 |  11.4 |1063k|   5 |  18 |  27 | 118 |  27 | 123 |  14 |   0 | 118 | 1.300000e+01 | 1.900000e+01 |  46.15%
 * * 0.1s|    39 |    28 |   386 |   7.0 |1092k|  14 |   - |  27 | 118 |  27 | 123 |  14 |   0 | 199 | 1.300000e+01 | 1.800000e+01 |  38.46%
 * ...
 * SCIP Status        : problem is solved [optimal solution found]
 * Solving Time (sec) : 0.75
 * Solving Nodes      : 4253
 * Primal Bound       : +1.80000000000000e+01 (287 solutions)
 * Dual Bound         : +1.80000000000000e+01
 * Gap                : 0.00 %
 *
 * SCIP>
 * \endcode
 *
 * We can navigate through the menus step-by-step and get a list of available options and submenus. Thus, we select
 * "set" to change settings, "heuristics" to change settings of primal heuristics, "shifting" for that particular
 * heuristic. Then we see a list of parameters (and yet another submenu for advanced parameters), and disable this
 * heuristic by setting its calling frequency to -1. If we already know the path to a certain setting, we can directly
 * type it (as for the rounding heuristic in the above example). Note that we do not have to use the full names, but we
 * may use short versions, as long as they are unique.
 *
 * To solve a problem a second time, we have to read it and start the optimization process again.
 *
 * \code
 * SCIP> set default
 * reset parameters to their default values
 * SCIP> set heuristics emphasis
 *
 *   aggressive            sets heuristics <aggressive>
 *   fast                  sets heuristics <fast>
 *   off                   turns <off> all heuristics
 *
 * SCIP/set/heuristics/emphasis> aggr
 * heuristics/veclendiving/freq = 5
 * ...
 * heuristics/crossover/minfixingrate = 0.5
 * SCIP> read check/instances/MIP/stein27.mps
 * original problem has 27 variables (27 bin, 0 int, 0 impl, 0 cont) and 118 constraints

 * SCIP> opt
 * ...
 * D 0.1s|     1 |     0 |   107 |     - | 971k|   0 |  24 |  27 | 122 |  27 | 131 |  13 |   4 |   0 | 1.300000e+01 | 1.800000e+01 |  38.46%
 *   0.1s|     1 |     0 |   107 |     - | 971k|   0 |  24 |  27 | 122 |  27 | 131 |  13 |   4 |   0 | 1.300000e+01 | 1.800000e+01 |  38.46%
 *   0.1s|     1 |     0 |   119 |     - |1111k|   0 |  24 |  27 | 122 |  27 | 132 |  14 |   4 |   0 | 1.300000e+01 | 1.800000e+01 |  38.46%
 *   0.1s|     1 |     2 |   119 |     - |1112k|   0 |  24 |  27 | 122 |  27 | 132 |  14 |   4 |  24 | 1.300000e+01 | 1.800000e+01 |  38.46%
 *  time | node  | left  |LP iter|LP it/n| mem |mdpt |frac |vars |cons |cols |rows |cuts |confs|strbr|  dualbound   | primalbound  |  gap
 *   0.2s|   100 |    59 |   698 |   5.8 |1138k|  14 |  11 |  27 | 122 |  27 | 123 |  14 |   4 | 204 | 1.300000e+01 | 1.800000e+01 |  38.46%
 *   0.2s|   200 |    91 |  1226 |   5.6 |1155k|  14 |   - |  27 | 122 |  27 | 123 |  14 |   4 | 207 | 1.300000e+01 | 1.800000e+01 |  38.46%
 * ^Cpressed CTRL-C 1 times (5 times for forcing termination)
 *
 * SCIP Status        : solving was interrupted [user interrupt]
 * Solving Time (sec) : 0.32
 * Solving Nodes      : 216
 * Primal Bound       : +1.80000000000000e+01 (283 solutions)
 * Dual Bound         : +1.30000000000000e+01
 * Gap                : 38.46 %
 *
 * SCIP>
 * \endcode
 *
 * Okay, what happened here? First, we reset all parameters to their default values, using "set default". Next, we
 * loaded some meta-parameter settings (also see <a href="FAQ.shtml#Section2">the FAQ</a>), to apply primal heuristics
 * more aggressively. SCIP shows us, which single parameters it changed therefor. Now, the optimal solution is already
 * found at the root node, by a heuristic which is deactivated by default.  Then, after node 200, the user pressed
 * CTRL-C which interrupts the solving process, We see that now in the short status report, primal and dual bound are
 * different, thus, the problem is not solved yet.  Nevertheless, we could access statistics, see the current incumbent
 * solution, change parameters and so on. Entering "optimize" we continue the solving process from the point on at which
 * it has been interrupted.
 *
 * SCIP can also write information to files. E.g., we could store the incumbent solution to a file, or output the
 * problem instance in another file format (the LP format is much more human readable than the MPS format, for example).
 *
 * \code
 * SCIP> write solution stein27.sol
 *
 * written solution information to file <stein27.sol>
 *
 * SCIP> write problem stein27.lp
 * written original problem to file <stein27.lp>
 *
 * SCIP> q
 * ...
 * \endcode
 *
 * We hope this tutorial gave you an overview of what is possible using the SCIP interactive shell. Please also read our
 * \ref FAQ, in particular the section <a href="FAQ.shtml#Section2">Using SCIP as a standalone MIP/MINLP-Solver</a>.
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

 * - The files \ref PUBLICMETHODS "pub_<...>.h" contain methods that perform "easy" operations that only
 *   affect the corresponding objects.
 *   Usually, with these methods you can access the data of the object.
 *   For example, in "pub_var.h" you find methods to get information about a variable.
 *
 * The file "pub_misc.h" contains methods for data structures like priority queues, hash tables, and hash maps,
 * as well as methods for sorting, numerics, random numbers, string operations, and file operations.
 *
 * If you are looking for a description of a callback method of a plugin that you want to implement, you have to
 * look at the corresponding \ref TYPEDEFINITIONS "type_<...>.h".
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CONS How to add constraint handlers
 *
 * A constraint handler defines the semantics and the algorithms to process constraints of a certain class.  A single
 * constraint handler is responsible for all constraints belonging to its constraint class.  For example, there is
 * one \ref cons_knapsack.h "knapsack constraint handler" that ensures solutions are only accepted if they satisfy all
 * knapsack constraints in the model. \n A complete list of all constraint handlers contained in this release can be
 * found \ref CONSHDLRS "here".
 *
 * We now explain how users can add their own constraint handlers.
 * For an example, look into the subtour constraint handler (examples/TSP/src/ConshdlrSubtour.cpp) of the
 * <a href="http://scip.zib.de/doc/examples/TSP/index.shtml">TSP </a> example project.
 * The example is written in C++ and uses the C++ wrapper classes.
 * However, we will explain the implementation of a constraint handler using the C interface.
 * It is very easy to transfer the C explanation to C++; whenever a method should be implemented using the
 * SCIP_DECL_CONS... notion, reimplement the corresponding virtual member function of the abstract scip::ObjConshdlr
 * base class.
 *
 * Additional documentation for the callback methods of a constraint handler can be found in the file
 * type_cons.h.
 *
 * Here is what you have to do (assuming your constraint handler should be named "subtour"):
 * -# Copy the template files src/scip/cons_xyz.c and src/scip/cons_xyz.h into files "cons_subtour.c"
 *    and "cons_subtour.h".
 *     \n
 *    Make sure to <b>adjust your Makefile</b> such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "subtour".
 * -# Adjust the \ref CONS_PROPERTIES "properties of the constraint handler".
 * -# Define the \ref CONS_DATA "constraint data and the constraint handler data". This is optional.
 * -# Implement the \ref CONS_INTERFACE "interface methods".
 * -# Implement the \ref CONS_FUNDAMENTALCALLBACKS "fundamental callback methods".
 * -# Implement the \ref CONS_ADDITIONALCALLBACKS "additional callback methods". This is optional.
 *
 *
 * @section CONS_PROPERTIES Properties of a Constraint Handler
 *
 * At the top of the new file "cons_subtour.c" you can find the constraint handler properties.
 * These are given as compiler defines. Some of them are optional, as, e.g., separation-related properties,
 * which only have to be defined if the constraint handler supports the related callbacks.
 * In the C++ wrapper class, you have to provide the constraint handler properties by calling the constructor
 * of the abstract base class scip::ObjConshdlr from within your constructor (see the TSP example).
 * The properties you have to set have the following meaning:
 *
 * @subsection CONS_FUNDAMENTALPROPERTIES Fundamental Constraint Handler properties
 *
 * \par CONSHDLR_NAME: the name of the constraint handler.
 * This name is used in the interactive shell to address the constraint handler.
 * Additionally, if you are searching for a constraint handler with SCIPfindConshdlr(), this name is looked up.
 * Names have to be unique: no two constraint handlers may have the same name.
 *
 * \par CONSHDLR_DESC: the description of the constraint handler.
 * This string is printed as a description of the constraint handler in the interactive shell of SCIP.
 *
 * \par CONSHDLR_ENFOPRIORITY: the priority of the constraint handler for constraint enforcing.
 * Like the separation priority, the enforcement priorities define the order in which the different constraint handlers
 * are called in the constraint enforcement step of the subproblem processing.
 * The constraint enforcement is called after the price-and-cut loop is executed (in the case that the LP is solved
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
 * \par CONSHDLR_NEEDSCONS: indicates whether the constraint handler should be skipped, if no constraints are available.
 * Usually, a constraint handler is only executed if there are constraints of its corresponding class in the model.
 * For those constraint handlers, the NEEDSCONS flag should be set to TRUE.
 * However, some constraint handlers must be called without having a constraint of the class in the model, because
 * the constraint is only implicitly available.
 * For example, the integrality constraint handler has the NEEDSCONS flag set to FALSE, because there is no explicit
 * integrality constraint in the model.
 * The integrality conditions are attached to the variables, and the integrality constraint handler has to check
 * all variables that are marked to be integer for integral values.
 *
 * @subsection CONS_ADDITIONALPROPERTIES Optional Constraint Handler properties
 *
 * The following properties are optional and only need to be defined if the constraint handlers support
 * separation, presolving, propagation, and/or upgrade functionality.
 *
 * \par LINCONSUPGD_PRIORITY: priority of the constraint handler for upgrading of linear constraints
 * This property is only needed if a certain linear constraint can be upgraded to a more specific one. In one of
 * the first presolving rounds SCIP tries to upgrade linear constraints to more specialized constraints, such as
 * knapsack constraints. The upgrading calls are processed in the order of decreasing priority.
 *
 * \par NONLINCONSUPGD_PRIORITY: priority of the constraint handler for upgrading of nonlinear constraints
 * This property has the same effect as the LINCONSUPGD_PRIORITY parameter, see above, and should be set whenever
 * an upgrade functionality from a general nonlinear constraint to the more specific one is defined.
 *
 * \par CONSHDLR_SEPAFREQ: the default frequency for separating cuts.
 * The separation frequency defines the depth levels at which the constraint handler's separation methods \ref CONSSEPALP
 * and \ref CONSSEPASOL are called.
 * For example, a separation frequency of 7 means, that the separation callback is executed for subproblems that are
 * in depth 0, 7, 14, ... of the branching tree.
 * A separation frequency of 0 means, that the separation method is only called at the root node.
 * A separation frequency of -1 disables the separation method of the constraint handler.
 * \n
 * The separation frequency can be adjusted by the user.
 * This property of the constraint handler only defines the default value of the frequency.
 * If you want to have a more flexible control of when to execute the separation algorithm, you have to assign
 * a separation frequency of 1 and implement a check at the beginning of your separation algorithm whether you really
 * want to execute the separator or not.
 * If you do not want to execute the method, set the result code to SCIP_DIDNOTRUN.
 *
 * \par CONSHDLR_SEPAPRIORITY: the priority of the constraint handler for separation. (optional: to be set only if the constraint handler supports separation)
 * In each separation round during the price-and-cut loop of the subproblem processing or during the separation loop
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
 * \par CONSHDLR_DELAYSEPA: the default for whether the separation method should be delayed, if other separators found cuts.
 * If the constraint handler's separation method is marked to be delayed, it is only executed after no other separator
 * or constraint handler found a cut during the price-and-cut loop.
 * If the separation method of the constraint handler is very expensive, you may want to mark it to be delayed until all
 * cheap separation methods have been executed.
 *
 * \par CONSHDLR_PROPFREQ: the default frequency for propagating domains.
 * This default frequency has the same meaning as the CONSHDLR_SEPAFREQ with respect to the domain propagation
 * callback of the constraint handler.
 * A propagation frequency of 0 means that propagation is only applied in preprocessing and at the root node.
 * A propagation frequency of -1 disables the propagation method of the constraint handler.
 *
 * \par CONSHDLR_DELAYPROP: the default for whether the propagation method should be delayed, if other propagators found reductions.
 * This property is analogous to the DELAYSEPA flag, but deals with the propagation method of the constraint handler.
 *
 * \par CONSHDLR_PROP_TIMING: the propagation timing mask of the constraint handler.
 * SCIP calls the domain propagation routines at different places in the node processing loop.
 * This property indicates at which places the propagation routine of the constraint handler is called.
 * Possible values are defined in type_timing.h and can be concatenated, e.g., as in SCIP_PROPTIMING_ALWAYS.
 *
 * \par CONSHDLR_MAXPREROUNDS: the default maximal number of presolving rounds the constraint handler participates in.
 * The preprocessing is executed in rounds.
 * If enough changes have been applied to the model, an additional preprocessing round is performed.
 * The MAXPREROUNDS parameter of a constraint handler denotes the maximal number of preprocessing rounds the constraint
 * handler participates in.
 * A value of -1 means that there is no limit on the number of rounds.
 * A value of 0 means the preprocessing callback of the constraint handler is disabled.
 *
 * \par CONSHDLR_DELAYPRESOL: the default for whether the presolving method should be delayed, if other presolvers found reductions.
 * This property is analogous to the DELAYSEPA flag, but deals with the preprocessing method of the constraint handler.
 *
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
 * At the bottom of "cons_subtour.c" you can find three interface methods, that also appear in "cons_subtour.h".
 * These are SCIPincludeConshdlrSubtour(), SCIPcreateConsSubtour(), and SCIPcreateConsSubtourBasic().
 * \n
 * The method SCIPincludeConshdlrSubtour() only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the constraint handler by calling the method
 * SCIPincludeConshdlr().
 * It is called by the user, if (s)he wants to include the constraint handler, i.e., if (s)he wants to make
 * the constraint handler available to the model, and looks like this:
 *  -# If you are using constraint handler data, you have to <b>allocate the memory for the data</b> at this point.
 *     You also have to initialize the fields in struct SCIP_ConshdlrData afterwards.
 *  \verbatim
 * SCIP_RETCODE SCIPincludeConshdlrKnapsack(
 * ...
 * )
 * {
 *    SCIP_EVENTHDLRDATA* eventhdlrdata;
 *    SCIP_CONSHDLRDATA* conshdlrdata;
 *    SCIP_CONSHDLR* conshdlr;
 *
 *  SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
 *  ...
 *  \endverbatim
 *  -# Now, <b>SCIP gets notified</b> of the presence of the constraint handler together with its \ref CONS_FUNDAMENTALCALLBACKS "basic callbacks".
 *   \code
 *  SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
 *        CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
 *        consEnfolpKnapsack, consEnfopsKnapsack, consCheckKnapsack, consLockKnapsack,
 *        conshdlrdata) );
 *  assert(conshdlr != NULL);
 *  \endcode
 *  -# All \ref CONS_ADDITIONALCALLBACKS "additional callbacks" are added via their setter functions.
 *  \code
 *  SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyKnapsack, consCopyKnapsack) );
 *  SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransKnapsack) );
 *  \endcode
 *  -# If the constraint handler is a specialization of a general linear or nonlinear constraint, we want to include an <b>automatic
 * upgrading mechanism</b> by calling the interface method
 *  \code
 *  if( SCIPfindConshdlr(scip,"linear") != NULL )
 *  {
 *       SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdKnapsack, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
 *  }
 *  \endcode
 *  or
 * \code
 * SCIP_CALL( SCIPincludeNonlinconsUpgrade(scip, nonlinconsUpgdSubtour, NULL, NONLINCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );
 * \endcode
 * in the nonlinear case.
 * See also cons_nonlinear.h for further information about the general upgrade procedure in the nonlinear case.
 *  -# You may also add <b>user parameters</b> for your constraint handler.
 * Some parameters which are important to play with are added to every constraint automatically, as, e.g.,
 * propagation or separation frequency.
 * \code
 *  SCIP_CALL( SCIPaddIntParam(scip,
 *        "constraints/knapsack/sepacardfreq",
 *        "multiplier on separation frequency, how often knapsack cuts are separated (-1: never, 0: only at root)",
 *        &conshdlrdata->sepacardfreq, TRUE, DEFAULT_SEPACARDFREQ, -1, INT_MAX, NULL, NULL) );
 *  ...
 *  return SCIP_OKAY;
 * }
 * \endcode
 *
 *
 *
 *
 * The methods SCIPcreateConsSubtour() and SCIPcreateConsSubtourBasic() are called to create a single constraint of the constraint
 * handler's constraint class.
 * It should allocate and fill the constraint data, and call SCIPcreateCons().
 * Take a look at the following example from the \ref cons_knapsack.h "knapsack constraint handler":
 *
 * \code
 * SCIP_RETCODE SCIPcreateConsKnapsack(
 *   SCIP*                 scip,
 *   SCIP_CONS**           cons,
 *   const char*           name,
 *   int                   nvars,
 *   SCIP_VAR**            vars,
 *   SCIP_Longint*         weights,
 *   SCIP_Longint          capacity,
 *   SCIP_Bool             initial,
 *   SCIP_Bool             separate,
 *   SCIP_Bool             enforce,
 *   SCIP_Bool             check,
 *   SCIP_Bool             propagate,
 *   SCIP_Bool             local,
 *   SCIP_Bool             modifiable,
 *   SCIP_Bool             dynamic,
 *   SCIP_Bool             removable,
 *   SCIP_Bool             stickingatnode
 *   )
 * {
 *    SCIP_CONSHDLRDATA* conshdlrdata;
 *    SCIP_CONSHDLR* conshdlr;
 *    SCIP_CONSDATA* consdata;
 *
 *    conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
 *    if( conshdlr == NULL )
 *    {
 *       SCIPerrorMessage("knapsack constraint handler not found\n");
 *       return SCIP_PLUGINNOTFOUND;
 *    }
 *
 *    conshdlrdata = SCIPconshdlrGetData(conshdlr);
 *    assert(conshdlrdata != NULL);
 *    assert(conshdlrdata->eventhdlr != NULL);
 *
 *    SCIP_CALL( consdataCreate(scip, &consdata, conshdlrdata->eventhdlr, nvars, vars, weights, capacity) );
 *
 *    SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
 *          local, modifiable, dynamic, removable, stickingatnode) );
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 *
 * In this example, consdataCreate() is a local method that allocates memory for the given consdata
 * and fills the data with the given <code>vars</code> array. For allocating memory for the constraint data, you
 * can use SCIP memory allocation:
 * \code
 * SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
 * \endcode
 *
 *
 * @section CONS_CALLBACKS Callback methods of Constraint handlers
 *
 * Besides the various functions which you will implement inside your constraint handler there exists a number
 * of <b> callback methods </b> associated with your constraint handler. Callback methods can be regarded as
 * tasks which your constraint handler is able to provide to the solver. They are grouped into two
 * categories:
 *
 * \ref CONS_FUNDAMENTALCALLBACKS "Fundamental Callback methods" are mandatory to implement
 * such that your code will work. For example, every constraint handler has to provide the
 * functionality to state whether all of its constraints are
 * fulfilled by a given variable assignment. Hence, the \ref CONSCHECK "CONSCHECK" callback is
 * one of the fundamental (or \a basic) callbacks of a constraint handler.
 *
 * Callbacks which are not necessarily implemented are grouped together as
 * \ref CONS_ADDITIONALCALLBACKS "additional callbacks". Such callbacks can be used to allocate and free memory
 * at different stages of the solving process. Although not mandatory, it might be useful to implement
 * some of these callbacks, e.g., to extend your constraint handler by a
 * \ref CONSSEPALP "separation" or \ref CONSPRESOL "presolving" functionality.
 *
 * All callbacks should be passed to SCIP during the SCIPinclude\<PLUGINTYPE\>\<PLUGINNAME\> method
 * (e.g., SCIPincludeConshdlrKnapsack() for the \ref cons_knapsack.h "knapsack constraint handler").
 * Since SCIP version 3.0, two ways of setting callbacks can be used, either via SCIPincludeConshdlr()
 * (all at once, as it always was), or via SCIPincludeConshdlrBasic() and setter functions for additional callbacks.
 * Since the basic inclusion methods are very unlikely to change and will thus
 * make your code more stable towards future versions of SCIP with more callbacks,
 * we recommend the latter choice, as explained in the \ref CONS_INTERFACE "interface" section.
 *
 * @section CONS_FUNDAMENTALCALLBACKS Fundamental Callback Methods
 *
 * By implementing the fundamental callbacks, you define the semantics of the constraint class the constraint handler
 * deals with.
 * If these methods are implemented, the resulting code is already correct and finds the optimal solution to the
 * given problem instance.
 * However, it might be very slow because the additional features, like cut separation and domain propagation, are
 * missing.
 * In the C++ wrapper class scip::ObjConshdlr, the fundamental callback methods are virtual abstract member functions.
 * You have to implement them in order to be able to construct an object of your constraint handler class.
 *
 * There are three fundamental callback methods that are all dealing with the feasibility of a given solution.
 * They are called at different places in the algorithm and have slightly different meaning.
 * However, it is usually reasonable to implement a single local method that is called by all of the three callback
 * methods with slightly modified parameters.
 * The fourth method provides dual information that is used for example in preprocessing.
 *
 * Additional documentation for the callback methods can be found in type_cons.h.
 *
 * @subsection CONSCHECK
 *
 * The CONSCHECK callback gets a primal solution candidate in a SCIP_SOL* data structure
 * and has to check this solution for global feasibility.
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
 * For example, the \ref cons_knapsack.h "knapsack constraint handler" loops over its constraints and
 * calculates the scalar product \f$w^T x\f$ of weights \f$w\f$ with the solution vector \f$x\f$.
 * This scalar product is compared with the capacity of the knapsack constraint.
 * If it exceeds the capacity, the CONSCHECK method is immediately aborted with the result SCIP_INFEASIBLE.
 * If all knapsack constraints are satisfied, a result SCIP_FEASIBLE is returned.
 *
 * @subsection CONSENFOLP
 *
 * The CONSENFOLP method is called after the price-and-cut loop was finished and an LP solution is available.
 * Like the CHECK call, the ENFOLP method should return a result SCIP_FEASIBLE, if the solution satisfies all the
 * constraints.
 * However, the behavior should be different, if the solution violates some of the associated constraints.
 * The constraint handler may return a result SCIP_INFEASIBLE in this situation, but this is not the best what
 * one can do.
 * The ENFOLP method has the possibility of \em resolving the infeasibility by
 * - stating that the current subproblem is infeasible (result SCIP_CUTOFF),
 * - adding an additional constraint that resolves the infeasibility (result SCIP_CONSADDED),
 * - reducing the domain of a variable (result SCIP_REDUCEDDOM),
 * - adding a cutting plane (result SCIP_SEPARATED),
 * - performing a branching (result SCIP_BRANCHED).
 *
 * However, the solution is not given as a SCIP_SOL* data structure.
 *
 * The value of a variable <code>var</code> in the LP solution can be accessed by calling
 * \code
 * SCIPgetVarSol(scip, var)
 * \endcode
 * or by
 * \code
 * SCIPgetSolVal(scip, NULL, var)
 * \endcode
 * By using the latter method, you can have a single local method to check a solution for feasibility by passing
 * the given <code>sol</code> to the CONSCHECK call and by passing a NULL pointer as <code>sol</code> to
 * the CONSENFOLP and CONSENFOPS calls.
 *
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
 *  - If the constraint may become violated by decreasing the value of a variable, it should call
 *    SCIPaddVarLocks(scip, var, nlockspos, nlocksneg), saying that rounding down is potentially rendering the
 *    (positive) constraint infeasible and rounding up is potentially rendering the negation of the constraint
 *    infeasible.
 *  - If the constraint may become violated by increasing the value of a variable, it should call
 *    SCIPaddVarLocks(scip, var, nlocksneg, nlockspos), saying that rounding up is potentially rendering the
 *    constraint's negation infeasible and rounding down is potentially rendering the constraint itself
 *    infeasible.
 *  - If the constraint may become violated by changing the variable in any direction, it should call
 *    SCIPaddVarLocks(scip, var, nlockspos + nlocksneg, nlockspos + nlocksneg).
 *
 *  <b>Note:</b> You do not have to worry about nlockspos and nlocksneg. These integer values are given as
 *  parameter of the CONSLOCK callback (see type_cons.h). Just use these variables in the above described
 *  fashion <b>without</b> adding or subtracting anything to them. In case of the knapsack constraints this
 *  method looks like this.
 *
 *  \code
 *  static
 *  SCIP_DECL_CONSLOCK(consLockKnapsack)
 *  {
 *     SCIP_CONSDATA* consdata;
 *     int i;
 *
 *     consdata = SCIPconsGetData(cons);
 *     assert(consdata != NULL);
 *
 *     for( i = 0; i < consdata->nvars; i++)
 *     {
 *        SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlocksneg, nlockspos) );
 *     }
 *
 *     return SCIP_OKAY;
 *  }
 * \endcode
 *
 *  To give same more intuition, consider the linear constraint \f$3x -5y +2z \leq 7\f$ as an example.
 *  The CONSLOCK callback method of the linear constraint handler should call
 *  SCIPaddVarLocks(scip, x, nlocksneg, nlockspos), SCIPaddVarLocks(scip, y, nlockspos, nlocksneg),
 *  and SCIPaddVarLocks(scip, z, nlocksneg, nlockspos) to tell SCIP,  that rounding up of \f$x\f$
 *  and \f$z\f$ and rounding down of \f$y\f$ can destroy the feasibility of the constraint, while rounding
 *  down of \f$x\f$ and \f$z\f$ and rounding up of \f$y\f$ can destroy the feasibility of the
 *  constraint's negation \f$3x -5y +2z > 7\f$.
 *  \n
 *  A linear constraint \f$2 \leq 3x -5y +2z \leq 7\f$ should call
 *  SCIPaddVarLocks(scip, ..., nlockspos + nlocksneg, nlockspos + nlocksneg) on all variables,
 *  since rounding in both directions of each variable can destroy both the feasibility of the
 *  constraint and it's negation \f$3x -5y +2z < 2\f$  or  \f$3x -5y +2z > 7\f$.
 *
 *
 * @section CONS_ADDITIONALCALLBACKS Additional Callback Methods
 *
 * The additional callback methods do not need to be implemented in every case, but provide useful functionality
 * for many applications. They can be added to your constraint handler via setter functions, see
 * \ref CONS_INTERFACE "here".
 *
 * @subsection CONSFREE
 *
 * If you are using constraint handler data, you have to implement this method in order to free the
 * constraint handler data. This can be done by the following procedure (which is taken from the
 * \ref cons_knapsack.h "knapsack constraint handler"):
 *
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
 *
 * If you have allocated memory for fields in your constraint handler data, remember to free this memory
 * before freeing the constraint handler data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection CONSHDLRCOPY
 *
 * The CONSHDLRCOPY callback is executed when the SCIP instance is copied, e.g. to solve a sub-SCIP. By defining this
 * callback as <code>NULL</code> the user disables the inclusion of the specified constraint handler into all copied SCIP
 * instances. This may deteriorate the performance of primal heuristics solving sub-SCIPs, since these constitute only
 * relaxations of the original problem if constraint handlers are missing.
 *
 * A usual implementation just
 * calls the interface method which includes the constraint handler to the model. For example, this callback is
 * implemented for the knapsack constraint handler as follows:
 *
 * \code
 * static
 * SCIP_DECL_CONSHDLRCOPY(conshdlrCopyKnapsack)
 * {
 *    assert(scip != NULL);
 *    assert(conshdlr != NULL);
 *    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
 *
 *    SCIP_CALL( SCIPincludeConshdlrKnapsack(scip) );
 *
 *    *valid = TRUE;
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 *
 * <b>Note:</b> If you implement this callback, take care when setting the valid pointer. The valid pointer should be
 * set to TRUE if (and only if!) you can make sure that all necessary data of the constraint handler are copied
 * correctly. If the complete problem is validly copied, i.e. if the copy methods of all problem defining plugins
 * (constraint handlers and pricers) return <code>*valid = TRUE</code>, then dual reductions found for the copied problem can be
 * transferred to the original SCIP instance. Thus, if the valid pointer is wrongly set to TRUE, it might happen that
 * optimal solutions are cut off.
 *
 * <b>Note:</b> If you implement this callback and the constraint handler needs constraints (see CONSHDLR_NEEDSCONS),
 * then you also need to implement the callback \ref CONSCOPY.
 *
 * @subsection CONSINIT
 *
 * The CONSINIT callback is executed after the problem is transformed.
 * The constraint handler may, e.g., use this call to replace the original variables in its constraints by transformed
 * variables, or to initialize its statistical constraint handler data.
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
 * before the branch-and-bound process begins.
 * Necessary constraint modifications that have to be performed even if presolving is turned off should be done here
 * or in the presolving initialization call.
 * Besides necessary modifications and clean up, no time consuming operations should be done.
 *
 * @subsection CONSINITSOL
 *
 * The CONSINITSOL callback is executed when the presolving is finished and the branch-and-bound process is about to
 * begin.
 * The constraint handler may use this call to initialize its branch-and-bound specific data.
 *
 * @subsection CONSEXITSOL
 *
 * The CONSEXITSOL callback is executed before the branch-and-bound process is freed.
 * The constraint handler should use this call to clean up its branch-and-bound data, in particular to release
 * all LP rows that it has created or captured.
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
 * Here is an example, which is taken from the \ref cons_knapsack.h "knapsack constraint handler":
 * \code
 * static
 * SCIP_DECL_CONSTRANS(consTransKnapsack)
 * {
 *    SCIP_CONSHDLRDATA* conshdlrdata;
 *    SCIP_CONSDATA* sourcedata;
 *    SCIP_CONSDATA* targetdata;
 *
 *    assert(conshdlr != NULL);
 *    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
 *    assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
 *    assert(sourcecons != NULL);
 *    assert(targetcons != NULL);
 *
 *    sourcedata = SCIPconsGetData(sourcecons);
 *    assert(sourcedata != NULL);
 *    assert(sourcedata->row == NULL);
 *
 *    conshdlrdata = SCIPconshdlrGetData(conshdlr);
 *    assert(conshdlrdata != NULL);
 *    assert(conshdlrdata->eventhdlr != NULL);
 *
 *    SCIP_CALL( consdataCreate(scip, &targetdata, conshdlrdata->eventhdlr,
 *          sourcedata->nvars, sourcedata->vars, sourcedata->weights, sourcedata->capacity) );
 *
 *    SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
 *          SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
 *          SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
 *          SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
 *          SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 *
 * @subsection CONSINITLP
 *
 * The CONSINITLP callback is executed before the first LP relaxation is solved.
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
 *  - detecting that the node is infeasible in the variables' bounds and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint (result SCIP_CONSADDED)
 *  - reducing a variable's domain (result SCIP_REDUCEDDOM)
 *  - adding a cutting plane to the LP (result SCIP_SEPARATED)
 *  - stating that the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the separator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the separator was skipped, but should be called again (result SCIP_DELAYED)
 *  - stating that a new separation round should be started without calling the remaining separator methods (result SCIP_NEWROUND)
 *
 * Please see also the @ref CONS_ADDITIONALPROPERTIES section to learn about the properties
 * CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, and CONSHDLR_DELAYSEPA, which influence the behaviour of SCIP
 * calling CONSSEPALP.
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
 *  - detecting that the node is infeasible in the variables' bounds and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint (result SCIP_CONSADDED)
 *  - reducing a variable's domain (result SCIP_REDUCEDDOM)
 *  - adding a cutting plane to the LP (result SCIP_SEPARATED)
 *  - stating that the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the separator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the separator was skipped, but should be called again (result SCIP_DELAYED)
 *  - stating that a new separation round should be started without calling the remaining separator methods (result SCIP_NEWROUND)
 *
 * Please see also the @ref CONS_ADDITIONALPROPERTIES section to learn about the properties
 * CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, and CONSHDLR_DELAYSEPA, which influence the behaviour of SCIP
 * calling CONSSEPASOL.
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
 *  - detecting that the node is infeasible in the variables' bounds and can be cut off (result SCIP_CUTOFF)
 *  - reducing a variable's domain (result SCIP_REDUCEDDOM)
 *  - stating that the propagator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the propagator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the propagator was skipped, but should be called again (result SCIP_DELAYED)
 *
 * Please see also the @ref CONS_ADDITIONALPROPERTIES section to learn about the properties
 * CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, and CONSHDLR_PROP_TIMING, which influence the behaviour of SCIP
 * calling CONSPROP.
 *
 * @subsection CONSRESPROP
 *
 * If the constraint handler should support \ref CONF "conflict analysis", it has to supply a CONSRESPROP method.
 * It also should call SCIPinferVarLbCons() or SCIPinferVarUbCons() in domain propagation instead of SCIPchgVarLb() or
 * SCIPchgVarUb() in order to deduce bound changes on variables.
 * In the SCIPinferVarLbCons() and SCIPinferVarUbCons() calls, the handler provides the constraint that deduced the
 * variable's bound change, and an integer value <code>inferinfo</code> that can be arbitrarily chosen.
 *
 * The propagation conflict resolving method CONSRESPROP must then be implemented to provide the "reasons" for the bound
 * changes, i.e., the bounds of variables at the time of the propagation, which forced the constraint to set the
 * conflict variable's bound to its current value. It can use the <code>inferinfo</code> tag to identify its own propagation rule
 * and thus identify the "reason" bounds. The bounds that form the reason of the assignment must then be provided by
 * calls to SCIPaddConflictLb() and SCIPaddConflictUb() in the propagation conflict resolving method.
 *
 * <b>Note:</b> The fact that <code>inferinfo</code> is an integer, as opposed to an arbitrary data object, is a compromise between space and speed. Sometimes a propagator would
 * need more information to efficiently infer the original propagation steps that lead to the conflict. This would,
 * however, require too much space. In the extreme, the original propagation steps have to be repeated.
 *
 * For example, the \ref cons_logicor.h "logicor constraint" \f$c = x \vee y \vee z\f$ fixes variable \f$z\f$ to TRUE (i.e., changes the lower
 * bound of \f$z\f$ to 1.0), if both, \f$x\f$ and \f$y\f$, are assigned to FALSE (i.e., if the upper bounds of these
 * variables are 0.0). It uses <code>SCIPinferVarLbCons(scip, z, 1.0, c, 0)</code> to apply this assignment (an
 * inference information tag is not needed by the constraint handler and is set to 0).  In the conflict analysis, the
 * constraint handler may be asked to resolve the lower bound change on \f$z\f$ with constraint \f$c\f$, that was
 * applied at a time given by a bound change index "bdchgidx".  With a call to <code>SCIPvarGetLbAtIndex(z,
 * bdchgidx)</code>, the handler can find out, that the lower bound of variable \f$z\f$ was set to 1.0 at the given
 * point of time, and should call <code>SCIPaddConflictUb(scip, x, bdchgidx)</code> and <code>SCIPaddConflictUb(scip, y,
 * bdchgidx)</code> to tell SCIP, that the upper bounds of \f$x\f$ and \f$y\f$ at this point of time were the reason for
 * the deduction of the lower bound of \f$z\f$.
 *
 * If conflict analysis should not be supported, the method has to set the result code to SCIP_DIDNOTFIND.  Although
 * this is a viable approach to circumvent the implementation of the usually rather complex conflict resolving method, it
 * will make the conflict analysis less effective. We suggest to first omit the conflict resolving method and check how
 * effective the \ref CONSPROP "propagation method" is. If it produces a lot of propagations for your application, you definitely should
 * consider implementing the conflict resolving method.
 *
 * @subsection CONSPRESOL
 *
 * The CONSPRESOL callback is called during preprocessing.
 * It should try to tighten the domains of the variables, tighten the coefficients of the constraints of the constraint
 * handler, delete redundant constraints, aggregate and fix variables if possible, and upgrade constraints to more
 * specific types.
 *
 * If the CONSPRESOL callback applies changes to the constraint data, you also have to implement the \ref CONSTRANS callback
 * in order to copy the constraint data to the transformed problem space and protect the original problem from the
 * preprocessing changes.
 *
 * To inform SCIP that the presolving method found a reduction the result pointer has to be set in a proper way.
 * The following options are possible:
 *
 *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in objective direction
 *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds
 *  - SCIP_SUCCESS    : the presolver found a reduction
 *  - SCIP_DIDNOTFIND : the presolver searched, but did not find a presolving change
 *  - SCIP_DIDNOTRUN  : the presolver was skipped
 *  - SCIP_DELAYED    : the presolver was skipped, but should be called again
 *
 * Please see also the @ref CONS_ADDITIONALPROPERTIES section to learn about the properties
 * CONSHDLR_MAXPREROUNDS and CONSHDLR_DELAYPRESOL, which influence the behaviour of SCIP
 * calling CONSPRESOL.
 *
 * @subsection CONSACTIVE
 *
 * The CONSACTIVE callback method is called each time a constraint of the constraint handler is activated.
 * For example, if a constraint is added locally to a subproblem, the CONSACTIVE callback is called whenever the
 * search enters the subtree where the constraint exists.
 *
 * @subsection CONSDEACTIVE
 *
 * The CONSDEACTIVE callback method is called each time a constraint of the constraint handler is deactivated.
 * For example, if a constraint is added locally to a subproblem, the CONSDEACTIVE callback is called whenever the
 * search leaves the subtree where the constraint exists.
 *
 * @subsection CONSENABLE
 *
 * The CONSENABLE callback method is called each time a constraint of the constraint handler is enabled.
 * Constraints might be active without being enabled. In this case, only the feasibility checks are executed,
 * but domain propagation and separation is skipped.
 *
 * @subsection CONSDISABLE
 *
 * The CONSDISABLE callback method is called each time a constraint of the constraint handler is disabled.
 *
 * @subsection CONSPRINT
 *
 * The CONSPRINT callback method is called, when the user asks SCIP to display the problem to the screen
 * or save the problem into a file. This is, however, only the case if the user requested the CIP format.
 * For more details about reading and writing with SCIP we refer to the \ref READER "file readers". In this
 * callback method the constraint handler should display the data of the constraint in an appropriate form.
 * The output format that is defined by the CONSPRINT callbacks is called CIP format.
 * In later versions of SCIP, the constraint handlers should also be able to parse (i.e., read) constraints
 * which are given in CIP format.
 *
 * @subsection CONSCOPY
 *
 * The CONSCOPY callback method is used whenever constraints should be copied from one SCIP instance into another SCIP
 * instance. This method comes with the necessary parameters to do so, most importantly with a mapping of the variables of the
 * source SCIP instance to the corresponding variables of the target SCIP instance, and a mapping for the constraints
 * in the same way. For a complete list of all arguments of this callback method see type_cons.h.
 *
 * To get the corresponding target variable of a given source variable, you can use the variable map directly:
 *
 * \code
 * targetvar = (SCIP_VAR*) (size_t) SCIPhashmapGetImage(varmap, sourcevar);
 * \endcode
 *
 * We recommend, however, to use the method SCIPgetVarCopy() which gets besides others the variable map and the constraint map as input
 * and returns the requested target variable. The advantage of using SCIPgetVarCopy() is that in the case
 * the required variable does not yet exist, it is created and added to the copy automatically:
 *
 * \code
 * SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevar, &targetvar, varmap, consmap, global) );
 * \endcode
 *
 * Finally, the result pointer <code>valid</code> has to be set to TRUE if (and only if!) the copy process was successful.
 *
 * <b>Note:</b> Be careful when setting the valid pointer. If you set the valid pointer to TRUE, but the constraint was
 * not copied one-to-one, then optimal solutions might be cut off during the search (see section
 * CONSHDLRCOPY above).
 *
 * For an example implementation we refer to cons_linear.h. Additional documentation and the complete list of all
 * parameters can be found in the file in type_cons.h.
 *
 * @subsection CONSPARSE
 *
 * This method is the counter part to CONSPRINT. The ideal idea is that a constraint handler is able to parse the output
 * which it generated via the CONSPRINT method and creates the corresponding constraint. If the parsing was successfully
 * the result pointer success should be set to TRUE. An example implementation can be found in the \ref cons_linear.h
 * "linear constraint handler".
 *
 * @subsection CONSDELVARS
 *
 * This method should iterate over the given constraints and delete all variables that were marked for deletion by SCIPdelVar().
 * Variable deletion is especially interesting for branch-cut-and-price applications. If your constraint handler allows
 * the addition of variables during the solving process (see "modifiable" attribute of constraints), then you might also want to
 * implement this callback. This would allow you to not only create variables during solving, but also remove them dynamically
 * from the problem to reduce memory consumption in case they are no longer necessary.
 * During presolving, SCIP may also find that some variables are not needed anymore and then try
 * to delete them. Thus, if you do not implement this callback, the constraint handler should capture its variables via
 * SCIPcaptureVar() to prevent SCIP from erroneously deleting them.
 *
 * Additional documentation and the complete list of all parameters can be found in the file type_cons.h.
 *
 * @subsection CONSGETVARS
 *
 * The CONSGETVARS callback of a constraint handler can be implemented to give access to the constraint variables
 * as array, independently from the internal data structure of the constraint. The buffer array
 * is already passed, together with its length. Consider implementing @ref CONSGETNVARS, too, to have
 * information about the number of variables in this constraint.
 *
 * @subsection CONSGETNVARS
 *
 * This callback can be implemented to return the number of variables involved into a particular constraint.
 * In order to have access to the variable pointers, consider implementing @ref CONSGETVARS.
 *
 * @section CONS_FURTHERINFO Further documentation
 *
 * Further documentation can be found in @ref type_cons.h for callback descriptions and a complete
 * list of all callback parameters, or in @ref scip.h
 * for globally available functions.
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
 * If the pricer finds one or more variables with negative reduced costs or negative Farkas value, it should
 * call SCIPcreateVar() and SCIPaddPricedVar() to create and add the variable to the problem. Additionally,
 * the pricer has to add the variable to all constraints in which it appears. Therefore, a pricer needs to
 * know the constraints of the model and their meaning. Note that all constraints for which additional variables
 * are generated by a pricer have to be flagged as "modifiable" in the SCIPcreateCons() call.
 *
 * We now explain how users can add their own pricers.
 * For example, look into the stable set pricer for the coloring problem (examples/Coloring/src/pricer_coloring.c) of the
 * Coloring example project.
 * The example is written in C. C++ users can easily adapt the code by using the scip::scip::ObjPricer wrapper base class and
 * implement the scip_...() virtual methods instead of the SCIP_DECL_PRICER... callback methods.
 *
 * Additional documentation for the callback methods of a pricer can be found in the file
 * type_pricer.h.
 *
 * Notice that if your pricer cannot cope with variable bounds other than 0 and infinity, you have to mark
 * all constraints containing priced variables as modifiable, and you may have to disable reduced cost
 * strengthening by setting propagating/rootredcost/freq to -1.
 *
 * Here is what you have to do to implement a pricer:
 * -# Copy the template files src/scip/pricer_xyz.c and src/scip/pricer_xyz.h into files "pricer_mypricer.c"
 *    and "pricer_mypricer.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "mypricer".
 * -# Adjust the properties of the pricer (see \ref PRICER_PROPERTIES).
 * -# Define the pricer data (see \ref PRICER_DATA). This is optional.
 * -# Implement the interface methods (see \ref PRICER_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref PRICER_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref PRICER_ADDITIONALCALLBACKS).  This is optional.
 *
 *
 * @section PRICER_PROPERTIES Properties of a Pricer
 *
 * At the top of the new file "pricer_mypricer.c" you can find the pricer properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the pricer properties by calling the constructor
 * of the abstract base class scip::ObjPricer from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par PRICER_NAME: the name of the pricer.
 * This name is used in the interactive shell to address the pricer.
 * Additionally, if you are searching for a pricer with SCIPfindPricer(), this name is looked up.
 * Names have to be unique: no two pricers may have the same name.
 *
 * \par PRICER_DESC: the description of the pricer.
 * This string is printed as a description of the pricer in the interactive shell.
 *
 * \par PRICER_PRIORITY: the priority of the pricer.
 * In each pricing round during the price-and-cut loop of the subproblem processing, the included pricers are
 * called in a predefined order, which is given by the priorities of the pricers.
 * The higher the priority, the earlier the pricer is called.
 * Usually, you will have only one pricer in your application and the priority is therefore irrelevant.
 *
 * \par PRICER_DELAY: the default for whether the pricer should be delayed, if other variables with negative reduced
 * costs have already been found in the current pricing round.
 * Variables may be declared to be "removable" in the SCIPcreateVar() call. This means that SCIP may remove the variable
 * from the LP if it was inactive (i.e., sitting at zero) for a number of LP solves. Nevertheless, after the removal of the
 * column from the LP, the variable still exists, and SCIP can calculate reduced costs and add it to the LP again if
 * necessary.
 * \n
 * If the PRICER_DELAY flag is set to TRUE (which is the common setting), all those existing variables with negative reduced costs
 * are added to the LP, and the LP is resolved before the pricer is called. Thus, the pricer can assume that all existing variables
 * have non-negative reduced costs if the \ref PRICERREDCOST method is called or non-positive Farkas value if the \ref PRICERFARKAS
 * method is called.
 * \n
 * In some applications, this inner pricing loop on the already existing variables can significantly slow down the solving process,
 * since it may lead to the addition of only very few variables in each pricing round. If this is an issue in your application,
 * you should consider setting the PRICER_DELAY flag to FALSE. You must, however, be aware of the fact that there may be already
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
 * If you are using C++, you can add pricer data, as usual, as object variables to your class.
 * \n
 * Defining pricer data is optional. You can leave the struct empty.
 *
 *
 * @section PRICER_INTERFACE Interface Methods
 *
 * At the bottom of "pricer_mypricer.c" you can find the interface method SCIPincludePricerMypricer(), which also appears in "pricer_mypricer.h".
 * It is called by the user, if (s)he wants to include the pricer, i.e., if (s)he wants to solve a model for which variables should
 * be generated by this pricer.
 *
 * This method only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the pricer. For this, you can either call SCIPincludePricer(),
 * or SCIPincludePricerBasic() since SCIP version 3.0. In the latter variant, \ref PRICER_ADDITIONALCALLBACKS "additional callbacks"
 * must be added via setter functions as, e.g., SCIPsetPricerCopy(). We recommend this latter variant because
 * it is more stable towards future SCIP versions which might have more callbacks, whereas source code using the first
 * variant must be manually adjusted with every SCIP release containing new callbacks for pricers in order to compile.
 *
 *
 * In addition, the pricer has to be activated before the solution process starts, like it is done
 * in the reader of the Coloring example (examples/Coloring/src/reader_col.c) by calling
 * \code
 * SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, "coloring")) );
 * \endcode
 *
 * If you are using pricer data, you have to allocate the memory for the data at this point.
 * You can do this by calling:
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &pricerdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_PricerData afterwards.
 *
 * You may also add user parameters for your pricer, see the method SCIPincludePricerColoring() in the pricer of the Coloring example
 * for an example of how to add user parameters.
 *
 *
 * @section PRICER_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Pricer
 *
 * The fundamental callback methods have to be implemented in order to obtain an operational algorithm.
 * They are passed together with the pricer itself to SCIP using SCIPincludePricer() or SCIPincludePricerBasic(),
 * see @ref PRICER_INTERFACE.
 *
 * In the case of a pricer, there are two fundamental callback methods, namely the @ref PRICERREDCOST and the
 * @ref PRICERFARKAS callbacks, which both search for new variables and add them to the problem.
 * These methods have to be implemented for every pricer; the other callback methods are optional.
 * In the C++ wrapper class scip::ObjPricer, the scip_redcost() method (which corresponds to the PRICERREDCOST callback)
 * is a virtual abstract member function. You have to implement it in order to be able to construct an object of your
 * pricer class.
 *
 * Additional documentation for the callback methods can be found in type_pricer.h.
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
 * methods of the constraint handlers to add the necessary variable entries to the constraints, see pub_cons.h.
 *
 * In the usual case that the pricer either adds a new variable or ensures that there are no further variables with negative dual feasibility,
 * the result pointer should be set to SCIP_SUCCESS. Only if the pricer aborts pricing without creating a new variable, but
 * there might exist additional variables with negative dual feasibility, the result pointer should be set to SCIP_DIDNOTRUN.
 * In this case, which sometimes is referred to as "early branching", the LP solution will not be used as a lower bound.
 * The pricer can, however, store a valid lower bound in the <code>lowerbound</code> pointer.
 *
 * Pricers usually need the dual LP solution as input for the pricing algorithm.
 * Since SCIP does not know the semantics of the individual constraints in the problem, the dual solution
 * has to be provided by the constraint handlers.
 * For example, the \ref cons_setppc.h "setppc constraint handler", which deals with set partitioning, packing, and covering constraints, provides
 * the method SCIPgetDualsolSetppc() to access the dual solution value for a single constraint.
 * Similarly, the dual solution of a linear constraint can be queried with the method SCIPgetDualsolLinear() of cons_linear.h.
 * The reduced costs of the existing variables can be accessed with the method SCIPgetVarRedcost().
 *
 * @subsection PRICERFARKAS
 *
 * If the current LP relaxation is infeasible, it is the task of the pricer to generate additional variables that can
 * potentially render the LP feasible again. In standard branch-and-price, these are variables with positive Farkas values,
 * and the PRICERFARKAS method should identify those variables.
 *
 * If the LP was proven to be infeasible, we have an infeasibility proof by the dual Farkas multipliers \f$y\f$.
 * With the values of \f$y\f$, an implicit inequality \f$y^T A x \ge y^T b\f$ is associated, with \f$b\f$ given
 * by the sides of the LP rows and the sign of \f$y\f$:
 *  - if \f$y_i\f$ is positive, \f$b_i\f$ is the left hand side of the row,
 *  - if \f$y_i\f$ is negative, \f$b_i\f$ is the right hand side of the row.
 *
 * \f$y\f$ is chosen in a way, such that the valid inequality  \f$y^T A x \ge y^T b\f$  is violated by all \f$x\f$,
 * especially by the (for this inequality least infeasible solution) \f$x'\f$ defined by
 *  - \f$x'_i := ub_i\f$, if \f$y^T A_i \ge 0\f$
 *  - \f$x'_i := lb_i\f$, if \f$y^T A_i < 0\f$.
 * Pricing in this case means to add variables \f$i\f$ with positive Farkas value, i.e., \f$y^T A_i x'_i > 0\f$.
 *
 * To apply Farkas pricing, the pricer needs to know the Farkas values of the constraints. Like the dual solution values for
 * feasible LP solutions, the dual Farkas values for infeasible solutions can be obtained by constraint handler interface
 * methods such as the SCIPgetDualfarkasLinear() method of the linear constraint handler.
 * The Farkas values for the bounds of the variables are just the regular reduced costs and can be accessed with SCIPgetVarRedcost().
 *
 * It is useful to note that Farkas pricing is the same as the regular pricing with a zero objective function.
 * Therefore, a typical implementation of a pricer would consist of a generic pricing algorithm that gets a dual solution and an
 * objective function vector as input and generates variables by calling SCIPcreateVar() and SCIPaddPricedVar().
 * The PRICERREDCOST callback would call this function with the regular objective function and the regular dual solution vector,
 * while the PRICERFARKAS callback would call this function with a zero objective function and the Farkas vector.
 * From a practical point of view, it is usually the simplest approach to provide just one Boolean flag to the generic pricing
 * algorithm in order to identify whether it is reduced cost or Farkas pricing. Then, the algorithm would just call the appropriate
 * methods to access the dual solution or objective function, depending on the Boolean flag.
 *
 * @section PRICER_ADDITIONALCALLBACKS Additional Callback Methods of a Pricer
 *
 * The additional callback methods do not need to be implemented in every case.
 * However, some of them have to be implemented for most applications. They can either be passed directly with
 * SCIPincludePricer() to SCIP or via specific <b>setter functions</b> after a call of SCIPincludePricerBasic(),
 * see also @ref PRICER_INTERFACE.
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
 * If you have allocated memory for fields in your pricer data, remember to free this memory
 * before freeing the pricer data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection PRICERCOPY
 *
 * The PRICERCOPY callback is executed when the SCIP instance is copied, e.g. to solve a sub-SCIP. By defining this
 * callback as <code>NULL</code> the user disables the inclusion of the pricer into all copied SCIP
 * instances. This means that primal heuristics will work on a sub-SCIP that contains only a part of the variables
 * and no variables are priced in during the solving process of the sub-SCIP. Therefore, primal solutions found in the
 * copied problem are typically still valid for the original problem and used for its solving process,
 * but dual reductions cannot be transferred to the original problem.
 *
 * <b>Note:</b> If you implement this callback, be careful when setting the valid pointer. The valid pointer should be
 * set to TRUE if (and only if!) you can make sure that all necessary data of the pricer are copied
 * correctly. If the complete problem is validly copied, i.e. if the copy methods of all problem defining plugins
 * (constraint handlers and pricers) return <code>*valid = TRUE</code>, then dual reductions found for the copied problem can be
 * transferred to the original SCIP instance. Thus, if the valid pointer is wrongly set to TRUE, it might happen that
 * optimal solutions are cut off.
 *
 * @subsection PRICERINIT
 *
 * The PRICERINIT callback is executed after the problem is transformed.
 * The pricer may, e.g., use this call to replace the original constraints stored in its pricer data by transformed
 * constraints, or to initialize other elements of its pricer data.
 *
 * @subsection PRICEREXIT
 *
 * The PRICEREXIT callback is executed before the transformed problem is freed.
 * In this method, the pricer should free all resources that have been allocated for the solving process in PRICERINIT.
 *
 * @subsection PRICERINITSOL
 *
 * The PRICERINITSOL callback is executed when the presolving is finished and the branch-and-bound process is about to begin.
 * The pricer may use this call to initialize its branch-and-bound specific data.
 *
 * @subsection PRICEREXITSOL
 *
 * The PRICEREXITSOL callback is executed before the branch-and-bound process is freed.
 * The pricer should use this call to clean up its branch-and-bound data, which was allocated in PRICERINITSOL.
 *
 * @section PRICER_REMARKS Further remarks
 *
 * If you use your own branching rule (e.g., to branch on constraints), make sure that it is able to branch on \a "pseudo solutions".
 * Otherwise, SCIP will use its default branching rules, if necessary (which all branch on variables). This
 * could disturb the pricing problem or branching might not even be possible, e.g., if all variables created thus far have already been fixed.
 *
 * Note that if the original problem is a maximization problem, SCIP will transform the problem into a minimization
 * problem by multiplying the objective function by -1. The pricer has to take care of this by multiplying
 * the original objective function value of all variables created during the solving process by -1.
 *
 * In some cases, bounds on variables are implicitly enforced by constraints of the problem and the objective function.
 * Therefore, these bounds do not need to be added to the LP explicitly, which has the advantage that the pricing routine does not need to
 * care about the corresponding dual values.
 * We call these bounds lazy bounds, they may be set by SCIPchgVarLbLazy() and SCIPchgVarUbLazy() for upper or lower bounds, respectively.
 * If the lazy bound is tighter than the local bound, the corresponding bound is not put into the LP.
 * In diving mode, lazy bounds are explicitly put into the LP, because changing the objective (which is only possible in diving)
 * might reverse the implicitly given bounds. When diving is finished, the bounds are again removed from the LP.
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
 * We now explain how users can add their own presolvers.
 * Take the dual fixing presolver (src/scip/presol_dualfix.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the scip::ObjPresol wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_PRESOL... callback methods.
 *
 * Additional documentation for the callback methods of a presolver, in particular for their input parameters,
 * can be found in the file type_presol.h.
 *
 * Here is what you have to do to implement a presolver:
 * -# Copy the template files src/scip/presol_xyz.c and src/scip/presol_xyz.h into files named "presol_mypresolver.c"
 *    and "presol_mypresolver.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "mypresolver".
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
 * of the abstract base class scip::ObjPresol from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par PRESOL_NAME: the name of the presolver.
 * This name is used in the interactive shell to address the presolver.
 * Additionally, if you are searching for a presolver with SCIPfindPresol(), this name is looked up.
 * Names have to be <b>unique</b>: no two presolvers may have the same name.
 *
 * \par PRESOL_DESC: the description of the presolver.
 * This string is printed as a description of the presolver in the interactive shell.
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
 * which also appears in "presol_mypresolver.h"
 * SCIPincludePresolMypresolver() is called by the user, if (s)he wants to include the presolver,
 * i.e., if (s)he wants to use the presolver in his/her application.
 *
 * This method only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the presolver. For this, you can either call SCIPincludePresol(),
 * or SCIPincludePresolBasic() since SCIP version 3.0. In the latter variant, \ref PRESOL_ADDITIONALCALLBACKS "additional callbacks"
 * must be added via setter functions as, e.g., SCIPsetPresolCopy(). We recommend this latter variant because
 * it is more stable towards future SCIP versions which might have more callbacks, whereas source code using the first
 * variant must be manually adjusted with every SCIP release containing new callbacks for presolvers in order to compile.
 *
 * If you are using presolver data, you have to allocate the memory for the data at this point.
 * You can do this by calling:
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &presoldata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_PresolData afterwards. For freeing the
 * presolver data, see \ref PRESOLFREE.
 *
 * You may also add user parameters for your presolver, see \ref PARAM for how to add user parameters and
 * the method SCIPincludePresolTrivial() in src/scip/presol_trivial.c for an example.
 *
 *
 * @section PRESOL_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Presolver
 *
 * The fundamental callback methods of the plugins are the ones that have to be implemented in order to obtain
 * an operational algorithm.
 * They are passed together with the presolver itself to SCIP using SCIPincludePresol() or SCIPincludePresolBasic(),
 * see @ref PRESOL_INTERFACE.
 *
 *  Presolver plugins have only one fundamental callback method, namely the @ref PRESOLEXEC method.
 * This method has to be implemented for every presolver; the other callback methods are optional.
 * In the C++ wrapper class scip::ObjPresol, the scip_exec() method (which corresponds to the PRESOLEXEC callback) is a virtual
 * abstract member function.
 * You have to implement it in order to be able to construct an object of your presolver class.
 *
 * Additional documentation for the callback methods, in particular to their input parameters,
 * can be found in type_presol.h.
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
 * The additional callback methods do not need to be implemented in every case. However, some of them have to be
 * implemented for most applications, they can be used, for example, to initialize and free private data.
 * Additional callbacks can either be passed directly with SCIPincludePresol() to SCIP or via specific
 * <b>setter functions</b> after a call of SCIPincludePresolBasic(), see also @ref PRESOL_INTERFACE.
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
 * If you have allocated memory for fields in your presolver data, remember to free this memory
 * before freeing the presolver data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection PRESOLINIT
 *
 * The PRESOLINIT callback is executed after the problem is transformed.
 * The presolver may, e.g., use this call to initialize its presolver data.
 * The difference between the original and the transformed problem is explained in
 * "What is this thing with the original and the transformed problem about?" on \ref FAQ.
 *
 * @subsection PRESOLCOPY
 *
 * The PRESOLCOPY callback is executed when a SCIP instance is copied, e.g. to
 * solve a sub-SCIP. By
 * defining this callback as
 * <code>NULL</code> the user disables the execution of the specified
 * presolver for all copied SCIP instances. This may deteriorate the performance
 * of primal heuristics using sub-SCIPs.
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
 * The PRESOLEXITPRE callback is executed after presolving finishes and before the branch-and-bound process begins.
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
 * We now explain how users can add their own separators.
 * Take the separator for the class of Gomory mixed integer inequalities (src/scip/sepa_gomory.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the scip::ObjSepa wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_SEPA... callback methods.
 *
 * Additional documentation for the callback methods of a separator, in particular for the input parameters,
 * can be found in the file type_sepa.h.
 *
 * Here is what you have to do to implement a separator:
 * -# Copy the template files src/scip/sepa_xyz.c and src/scip/sepa_xyz.h into files "sepa_myseparator.c"
 *    and "sepa_myseparator.h".
      \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "myseparator".
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
 * of the abstract base class scip::ObjSepa from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par SEPA_NAME: the name of the separator.
 * This name is used in the interactive shell to address the separator.
 * Additionally, if you are searching for a separator with SCIPfindSepa(), this name is looked up.
 * Names have to be unique: no two separators may have the same name.
 *
 * \par SEPA_DESC: the description of the separator.
 * This string is printed as a description of the separator in the interactive shell.
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
 * The frequency can be adjusted by the user. This property of the separator only defines the default value of the frequency.
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
 * \par SEPA_USESSUBSCIP: Does the separator use a secondary SCIP instance?
 * Some heuristics and separators solve MIPs or SAT problems and use a secondary SCIP instance. Examples are
 * Large Neighborhood Search heuristics such as RINS and Local Branching or the CGMIP separator. To avoid recursion,
 * these plugins usually deactivate all other plugins that solve MIPs. If a separator uses a secondary SCIP instance,
 * this parameter has to be TRUE and it is recommended to call SCIPsetSubscipsOff() for the secondary SCIP instance.
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
 * At the bottom of "sepa_myseparator.c", you can find the interface method SCIPincludeSepaMyseparator(),
 * which also appears in "sepa_myseparator.h"
 * SCIPincludeSepaMyseparator() is called by the user, if (s)he wants to include the separator,
 * i.e., if (s)he wants to use the separator in his/her application.
 *
 * This method only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the separator. For this, you can either call SCIPincludeSepa(),
 * or SCIPincludeSepaBasic() since SCIP version 3.0. In the latter variant, \ref SEPA_ADDITIONALCALLBACKS "additional callbacks"
 * must be added via setter functions as, e.g., SCIPsetSepaCopy(). We recommend this latter variant because
 * it is more stable towards future SCIP versions which might have more callbacks, whereas source code using the first
 * variant must be manually adjusted with every SCIP release containing new callbacks for separators in order to compile.
 *
 * If you are using separator data, you have to allocate the memory
 * for the data at this point. You can do this by calling:
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
 * an operational algorithm.
 * They are passed together with the separator itself to SCIP using SCIPincludeSepa() or SCIPincludeSepaBasic(),
 * see @ref SEPA_INTERFACE.
 *
 * Separator plugins have two callbacks, @ref SEPAEXECLP and @ref SEPAEXECSOL, of which at least one must be implemented.
 *
 * Additional documentation for the callback methods, in particular to their input parameters,
 * can be found in type_sepa.h.
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
 * the 'result' variable (see type_sepa.h):
 *  - detecting that the node is infeasible in the variable's bounds and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint (result SCIP_CONSADDED)
 *  - reducing a variable's domain (result SCIP_REDUCEDDOM)
 *  - adding a cutting plane to the LP (result SCIP_SEPARATED)
 *  - stating that the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the separator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the separator was skipped, but should be called again (result SCIP_DELAYED)
 *  - stating that a new separation round should be started without calling the remaining separator methods (result SCIP_NEWROUND)
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
 * the 'result' variable (see type_sepa.h):
 *  - detecting that the node is infeasible in the variable's bounds and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint (result SCIP_CONSADDED)
 *  - reducing a variable's domain (result SCIP_REDUCEDDOM)
 *  - adding a cutting plane to the LP (result SCIP_SEPARATED)
 *  - stating that the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the separator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the separator was skipped, but should be called again (result SCIP_DELAYED)
 *  - stating that a new separation round should be started without calling the remaining separator methods (result SCIP_NEWROUND)
 *
 *
 * @section SEPA_ADDITIONALCALLBACKS Additional Callback Methods of a Separator
 *
 * The additional callback methods do not need to be implemented in every case. However, some of them have to be
 * implemented for most applications, they can be used, for example, to initialize and free private data.
 * Additional callbacks can either be passed directly with SCIPincludeSepa() to SCIP or via specific
 * <b>setter functions</b> after a call of SCIPincludeSepaBasic(), see also @ref SEPA_INTERFACE.
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
 * If you have allocated memory for fields in your separator data, remember to free this memory
 * before freeing the separator data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection SEPACOPY
 *
 * The SEPACOPY callback is executed when a SCIP instance is copied, e.g. to
 * solve a sub-SCIP. By
 * defining this callback as
 * <code>NULL</code> the user disables the execution of the specified
 * separator for all copied SCIP instances. This may deteriorate the performance
 * of primal heuristics using sub-SCIPs.
 *
 * @subsection SEPAINIT
 *
 * The SEPAINIT callback is executed after the problem is transformed.
 * The separator may, e.g., use this call to initialize its separator data.
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
 * The SEPAINITSOL callback is executed when the presolving is finished and the branch-and-bound process is about to
 * begin. The separator may use this call to initialize its branch-and-bound specific data.
 *
 * @subsection SEPAEXITSOL
 *
 * The SEPAEXITSOL callback is executed before the branch-and-bound process is freed. The separator should use this call
 * to clean up its branch-and-bound data, in particular to release all LP rows that it has created or captured.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page PROP How to add propagators
 *
 * Propagators are used to tighten the domains of the variables. Like for cutting planes, there are two different types
 * of domain propagations. Constraint based (primal) domain propagation algorithms are part of the corresponding
 * constraint handlers, see \ref CONSPROP. In contrast, domain propagators usually provide dual propagations, i.e.,
 * propagations that can be applied using the objective function and the current best known primal solution. This
 * section deals with such propagators.
 *
 * A complete list of all propagators contained in this release can be found \ref PROPAGATORS "here".
 *
 * We now explain how users can add their own propagators.  Take the pseudo objective function propagator
 * (src/scip/prop_pseudoobj.c) as an example.  As all other default plugins, it is written in C. C++ users can easily
 * adapt the code by using the scip::ObjProp wrapper base class and implement the @c scip_...() virtual methods instead
 * of the @c SCIP_DECL_PROP... callback methods.
 *
 * Additional documentation for the callback methods of a propagator can be found in the file type_prop.h.
 *
 * Here is what you have to do to implement a propagator:
 * -# Copy the template files src/scip/prop_xyz.c and src/scip/prop_xyz.h into files named "prop_mypropagator.c"
 *    and "prop_mypropagator.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "mypropagator".
 * -# Adjust the properties of the propagator (see \ref PROP_PROPERTIES).
 * -# Define the propagator data (see \ref PROP_DATA). This is optional.
 * -# Implement the interface methods (see \ref PROP_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref PROP_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref PROP_ADDITIONALCALLBACKS). This is optional.
 *
 * @section PROP_PROPERTIES Properties of a Propagator
 *
 * At the top of the new file "prop_mypropagator.c" you can find the propagator properties. These are given as compiler
 * defines. The presolving-related properties are optional,
 * they only have to be defined if the propagator supports presolving routines.
 * In the C++ wrapper class, you have to provide the propagator properties by calling the constructor of the
 * abstract base class scip::ObjProp from within your constructor.  The properties you have the following meaning:
 *
 * @subsection PROP_FUNDAMENTALPROPERTIES Fundamental properties of a propagator
 *
 * \par PROP_NAME: the name of the propagator.
 * This name is used in the interactive shell to address the propagator.  Additionally, if you are searching for a
 * propagator with SCIPfindProp(), this name is searched for.  Names have to be unique: no two propagators may have the
 * same name.
 *
 * \par PROP_DESC: the description of the propagator.
 * This string is printed as a description of the propagator in the interactive shell.
 *
 * \par PROP_PRIORITY: the priority of the propagator.
 * In each propagation round, the propagators and propagation methods of the constraint handlers are called in a
 * predefined order, which is given by the priorities of the propagators and the check priorities of the constraint
 * handlers.  First, the propagators with non-negative priority are called in order of decreasing priority.  Next, the
 * propagation methods of the different constraint handlers are called in order of decreasing check priority.  Finally,
 * the propagators with negative priority are called in order of decreasing priority.  \n The priority of the
 * propagators should be set according to the complexity of the propagation algorithm and the impact of the domain
 * propagations: propagators providing fast algorithms that usually have a high impact (i.e., tighten many bounds)
 * should have a high priority.
 *
 * \par PROP_FREQ: the default frequency for propagating domains.
 * The frequency defines the depth levels at which the propagation method \ref PROPEXEC is called.  For example, a
 * frequency of 7 means, that the propagation callback is executed for subproblems that are in depth 0, 7, 14, ... of
 * the branching tree. A frequency of 0 means that propagation is only applied in preprocessing and at the root node. A
 * frequency of -1 disables the propagator.
 * \n
 * The frequency can be adjusted by the user. This property of the propagator only defines the default value of the
 * frequency.\n
 * <b>Note:</b> If you want to have a more flexible control of when to execute the propagation algorithm, you have to
 * assign a frequency of 1 and implement a check at the beginning of your propagation algorithm whether you really want
 * to execute the domain propagation or not. If you do not want to execute it, set the result code to SCIP_DIDNOTRUN.
 *
 * \par PROP_DELAY: the default for whether the propagation method should be delayed, if other propagators or constraint handlers found domain reductions.
 * If the propagator's propagation method is marked to be delayed, it is only executed after no other propagator or
 * constraint handler found a domain reduction in the current iteration of the domain propagation loop.  If the
 * propagation method of the propagator is very expensive, you may want to mark it to be delayed until all cheap
 * propagation methods have been executed.
 *
 * \par PROP_TIMING: the timing mask of the propagator.
 * SCIP calls the domain propagation routines at different places in the node processing loop.
 * This property indicates at which places the propagator is called.
 * Possible values are defined in type_timing.h and can be concatenated, e.g., as in SCIP_PROPTIMING_ALWAYS.
 *
 * @subsection PROP_ADDITIONALPROPERTIES Optional propagator properties
 *
 * The following properties are optional and only need to be defined if the propagator supports
 * presolving, that is, if the \ref PROPPRESOL "presolving callback" is implemented.
 *
 * \par PROP_PRESOL_PRIORITY: the priority of the presolving method.
 * This attribute is analogous to the PROP_PRIORITY flag, but deals with the preprocessing method of the presolver.
 *
 * \par PROP_PRESOL_MAXROUNDS: the default maximal number of presolving rounds the propagator participates in.
 * The preprocessing is executed in rounds.
 * If enough changes have been applied to the model, an additional preprocessing round is performed.
 * The MAXROUNDS parameter of a propagator denotes the maximal number of preprocessing rounds, the propagator
 * participates in.
 * A value of -1 means, that there is no limit on the number of rounds.
 * A value of 0 means, the preprocessing callback of the propagator is disabled.
 *
 * \par PROP_PRESOL_DELAY: the default for whether the presolving method should be delayed, if other propagators or constraint handlers found presolving reductions.
 * This property is analogous to the PROP_DELAY flag, but deals with the preprocessing method of the propagator.
 *
 * @section PROP_DATA Propagator Data
 *
 * Below the title "Data structures" you can find a struct called <code>struct SCIP_PropData</code>.  In this data
 * structure, you can store the data of your propagator. For example, you should store the adjustable parameters of the
 * propagator in this data structure.  If you are using C++, you can add propagator data as object variables to your
 * class as usual .
 * \n
 * Defining propagator data is optional. You can leave the struct empty.
 *
 *
 * @section PROP_INTERFACE Interface Methods
 *
 * At the bottom of "prop_mypropagator.c", you can find the interface method SCIPincludeSepaMypropagator(),
 * which also appears in "prop_mypropagator.h"
 * SCIPincludePropMypropagator() is called by the user, if (s)he wants to include the propagator,
 * i.e., if (s)he wants to use the propagator in his/her application.
 *
 * This method only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the propagator. For this, you can either call SCIPincludeProp(),
 * or SCIPincludePropBasic() since SCIP version 3.0. In the latter variant, \ref PROP_ADDITIONALCALLBACKS "additional callbacks"
 * must be added via setter functions as, e.g., SCIPsetPropCopy(). We recommend this latter variant because
 * it is more stable towards future SCIP versions which might have more callbacks, whereas source code using the first
 * variant must be manually adjusted with every SCIP release containing new callbacks for separators in order to compile.
 *
 *
 * If you are using propagator data, you have to allocate the memory for the data at this point.  You can do this by
 * calling
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &propdata) );
 * \endcode
 * You also have to initialize the fields in <code>struct SCIP_PropData</code> afterwards.
 *
 * You may also add user parameters for your propagator, see the method SCIPincludePropPseudoobj() in
 * src/scip/prop_pseudoobj.c for an example.
 *
 *
 * @section PROP_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Propagator
 *
 * The fundamental callback methods of the plugins are the ones that have to be implemented in order to obtain
 * an operational algorithm.
 * They are passed together with the propagator itself to SCIP using SCIPincludeProp() or SCIPincludePropBasic(),
 * see @ref PROP_INTERFACE.
 *
 * Propagator plugins have one fundamental callback method, namely the \ref PROPEXEC method
 * method.  This method has to be implemented for every propagator; the other callback methods are optional.  In the
 * C++ wrapper class scip::ObjProp, the scip_exec() method (which corresponds to the \ref PROPEXEC
 * callback) is a virtual abstract member function. You have to
 * implement it in order to be able to construct an object of your propagator class.
 *
 * Additional documentation for the callback methods can be found in type_prop.h.
 *
 * @subsection PROPEXEC
 *
 * The PROPEXEC callback is called during presolving and during the subproblem processing. It should perform the actual
 * domain propagation, which means that it should tighten the variables' bounds.  The technique of domain propagation,
 * which is the main workhorse of constraint programming, is called "node preprocessing" in the Integer Programming
 * community.
 *
 * The PROPEXEC callback has the following options:
 *  - detecting that the node is infeasible in the variables' bounds and can be cut off (result SCIP_CUTOFF)
 *  - reducing (i.e, tightening) the domains of some variables (result SCIP_REDUCEDDOM)
 *  - stating that the propagator searched, but did not find domain reductions, cutting planes, or cut constraints
 *    (result SCIP_DIDNOTFIND)
 *  - stating that the propagator was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the propagator was skipped, but should be called again (result SCIP_DELAYED)
 *
 *
 *
 * @section PROP_ADDITIONALCALLBACKS Additional Callback Methods of a Propagator
 *
 * The additional callback methods do not need to be implemented in every case. However, some of them have to be
 * implemented for most applications, they can be used, for example, to initialize and free private data.
 * Additional callbacks can either be passed directly with SCIPincludeProp() to SCIP or via specific
 * <b>setter functions</b> after a call of SCIPincludePropBasic(), see also @ref PROP_INTERFACE.
 *
 * @subsection PROPRESPROP
 *
 * If the propagator wants to support \ref CONF "conflict analysis", it has to supply the PROPRESPROP method.  It also should call
 * SCIPinferVarLbProp() or SCIPinferVarUbProp() in the domain propagation instead of SCIPchgVarLb() or SCIPchgVarUb() in
 * order to deduce bound changes on variables.  In the SCIPinferVarLbProp() and SCIPinferVarUbProp() calls, the
 * propagator provides a pointer to itself and an integer value "inferinfo" that can be arbitrarily chosen.
 *
 * The propagation conflict resolving method PROPRESPROP must then be implemented to provide the "reasons" for the bound
 * changes, i.e., the bounds of variables at the time of the propagation, which forced the propagator to set the
 * conflict variable's bound to its current value. It can use the "inferinfo" tag to identify its own propagation rule
 * and thus identify the "reason" bounds. The bounds that form the reason of the assignment must then be provided by
 * calls to SCIPaddConflictLb() and SCIPaddConflictUb() in the propagation conflict resolving method.
 *
 * See the description of the propagation conflict resolving method \ref CONSRESPROP of constraint handlers for
 * further details.
 *
 * Omitting the PROPRESPROP callback circumvents the implementation of the usually rather complex conflict resolving method.
 * Yet, it
 * will make the conflict analysis less effective. We suggest to first omit the conflict resolving method and check how
 * effective the propagation method is. If it produces a lot of propagations for your application, you definitely should
 * consider implementing the conflict resolving method.
 *
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
 *    SCIPfreeMemory(scip, &propdata);
 *
 *    SCIPpropSetData(prop, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you have allocated memory for fields in your propagator data, remember to free this memory
 * before freeing the propagator data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection PROPINIT
 *
 * The PROPINIT callback is executed after the problem is transformed.  The propagator may, e.g., use this call to
 * initialize its propagator data.
 *
 * @subsection PROPCOPY
 *
 * The PROPCOPY callback is executed when a SCIP instance is copied, e.g. to
 * solve a sub-SCIP. By
 * defining this callback as
 * <code>NULL</code> the user disables the execution of the specified
 * propagator for all copied SCIP instances. This may deteriorate the performance
 * of primal heuristics using sub-SCIPs.
 *
 * @subsection PROPEXIT
 *
 * The PROPEXIT callback is executed before the transformed problem is freed.
 * In this method, the propagator should free all resources that have been allocated for the solving process in PROPINIT.
 *
 * @subsection PROPINITPRE
 *
 * The PROPINITPRE callback is executed before the preprocessing is started, even if presolving is turned off.
 * The propagator may use this call to initialize its presolving data before the presolving process begins.
 *
 * @subsection PROPEXITPRE
 *
 * The PROPEXITPRE callback is executed after the preprocessing has been finished, even if presolving is turned off.
 * The propagator may use this call, e.g., to clean up its presolving data.
 * Besides clean up, no time consuming operations should be done.
 *
 * @subsection PROPINITSOL
 *
 * The PROPINITSOL callback is executed when the presolving is finished and the branch-and-bound process is about to
 * begin.
 * The propagator may use this call to initialize its branch-and-bound specific data.
 *
 * @subsection PROPEXITSOL
 *
 * The PROPEXITSOL callback is executed before the branch-and-bound process is freed.
 * The propagator should use this call to clean up its branch-and-bound data.
 *
 * @subsection PROPPRESOL
 *
 * Seaches for domain propagations, analogous to the \ref PROPEXEC callback.
 * However, this callback is called during preprocessing.
 *
 * To inform SCIP that the presolving method found a reduction the result pointer has to be set in a proper way.
 * The following options are possible:
 *
 *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in objective direction
 *  - SCIP_CUTOFF     : at least one domain reduction that renders the problem infeasible has been found
 *  - SCIP_SUCCESS    : the presolver found a domain reduction
 *  - SCIP_DIDNOTFIND : the presolver searched, but did not find a presolving change
 *  - SCIP_DIDNOTRUN  : the presolver was skipped
 *  - SCIP_DELAYED    : the presolver was skipped, but should be called again
 *
 *
 * Please see also the @ref PROP_ADDITIONALPROPERTIES section to learn about the properties
 * PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, and PROP_PRESOL_DELAY, which influence the behaviour of SCIP
 * calling PROPPRESOL.
 *
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page BRANCH How to add branching rules
 *
 * Branching rules are used to split the problem at the current node into smaller subproblems. Branching rules can be called at three
 * different occasions, which is why they have three different execution methods (see \ref
 * BRANCHRULE_ADDITIONALCALLBACKS).  Branching is performed if:
 * - the LP solution of the current problem is fractional. In this case, the integrality constraint handler calls the
 *   \ref BRANCHEXECLP methods of the branching rules.
 * - the list of external branching candidates is not empty. This will only be the case if branching candidates were added
 *   by a user's \ref RELAX "relaxation handler" or \ref CONS "constraint handler" plugin, calling SCIPaddExternBranchCand().
 *   These branching candidates should be processed by the \ref BRANCHEXECEXT method.
 * - if an integral solution violates one or more constraints and this infeasibility could not be resolved in the callback methods
 *   \ref CONSENFOLP and \ref CONSENFOPS of the corresponding constraint handlers. In this case, the \ref BRANCHEXECPS method will be called. This is the
 *   standard case, if you use SCIP as a pure CP or SAT solver. If the LP or any other type of relaxation is used, then
 *   branching on pseudo solutions works as a last resort.
 *
 * The idea of branching rules is to take a global view on the problem. In contrast, branching paradigms which are
 * specific to one type of constraint are best implemented within the enforcement callbacks of your constraint handler.
 * See, e.g., the constraint specific branching rules provided by the constraint handlers for special ordered sets
 * (src/scip/cons_sos{1,2}.c)).
 * \n
 * All branching rules that come with the default distribution of SCIP create two subproblems by splitting a single
 * variable's domain.  It is, however, fully supported to implement much more general branching schemes, for example by
 * creating more than two subproblems, or by adding additional constraints to the subproblems instead of tightening the
 * domains of the variables.
 * \n
 * A complete list of all branching rules contained in this release can be found \ref BRANCHINGRULES "here".
 *
 * We now explain how users can add their own branching rules.  Take the most infeasible LP branching rule
 * (src/scip/branch_mostinf.c) as an example.  As all other default plugins, it is written in C. C++ users can easily
 * adapt the code by using the scip::ObjBranchrule wrapper base class and implement the scip_...() virtual methods instead of
 * the SCIP_DECL_BRANCH... callback methods.
 *
 * Additional documentation for the callback methods of a branching rule can be found in the file type_branch.h.
 *
 * Here is what you have to do to implement a branching rule:
 * -# Copy the template files src/scip/branch_xyz.c and src/scip/branch_xyz.h into files named
 *    "branch_mybranchingrule.c" and "branch_mybranchingrule.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "mybranchingrule".
 * -# Adjust the properties of the branching rule (see \ref BRANCHRULE_PROPERTIES).
 * -# Define the branching rule data (see \ref BRANCHRULE_DATA). This is optional.
 * -# Implement the interface methods (see \ref BRANCHRULE_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref BRANCHRULE_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref BRANCHRULE_ADDITIONALCALLBACKS). This is optional.
 *
 *
 * @section BRANCHRULE_PROPERTIES Properties of a Branching Rule
 *
 * At the top of the new file "branch_mybranchingrule.c" you can find the branching rule properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the branching rule properties by calling the constructor
 * of the abstract base class scip::ObjBranchrule from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par BRANCHRULE_NAME: the name of the branching rule.
 * This name is used in the interactive shell to address the branching rule.
 * Additionally, if you are searching for a branching rule with SCIPfindBranchrule(), this name is looked up.
 * Names have to be unique: no two branching rules may have the same name.
 *
 * \par BRANCHRULE_DESC: the description of the branching rule.
 * This string is printed as a description of the branching rule in the interactive shell.
 *
 * \par BRANCHRULE_PRIORITY: the default value for the priority of the branching rule.
 * In the subproblem processing, the branching rules are called in decreasing order of their priority until
 * one succeeded to branch. Since most branching rules are able to generate a branching in all situations,
 * only the rule of highest priority is used. In combination with the BRANCHRULE_MAXDEPTH and
 * BRANCHRULE_MAXBOUNDDIST settings, however, interesting strategies can be easily employed. For example,
 * the user can set the priority of the "full strong branching" strategy to the highest value and assign the
 * second highest value to the "reliable pseudo cost" rule. If (s)he also sets the maximal depth for the
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
 * \par BRANCHRULE_MAXBOUNDDIST: the default value for the maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for applying branching.
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
 * At the bottom of "branch_mybranchingrule.c", you can find the interface method SCIPincludeBranchruleMybranchingrule(),
 * which also appears in "branch_mybranchingrule.h"
 * SCIPincludeBranchruleMybranchingrule() is called by the user, if (s)he wants to include the branching rule,
 * i.e., if (s)he wants to use the branching rule in his/her application.
 *
 * This method only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the branching rule. For this, you can either call
 * SCIPincludeBranchrule(),
 * or SCIPincludeBranchruleBasic() since SCIP version 3.0. In the latter variant, \ref BRANCHRULE_ADDITIONALCALLBACKS "additional callbacks"
 * must be added via setter functions as, e.g., SCIPsetBranchruleCopy(). We recommend this latter variant because
 * it is more stable towards future SCIP versions which might have more callbacks, whereas source code using the first
 * variant must be manually adjusted with every SCIP release containing new callbacks for branchrule in order to compile.
 *
 *
 * If you are using branching rule data, you have to allocate the memory for the data at this point.
 * You can do this by calling:
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
 * The additional callback methods do not need to be implemented in every case. However, some of them have to be
 * implemented for most applications, they can be used, for example, to initialize and free private data.
 * Additional callbacks can either be passed directly with SCIPincludeBranchrule() to SCIP or via specific
 * <b>setter functions</b> after a call of SCIPincludeBranchruleBasic(), see also @ref BRANCHRULE_INTERFACE.
 *
 * The most important callback methods are the \ref BRANCHEXECLP, \ref BRANCHEXECEXT,
 * and \ref BRANCHEXECPS methods, which perform the actual task of generating a branching.
 *
 * Additional documentation for the callback methods can be found in type_branch.h.
 *
 * @subsection BRANCHEXECLP
 *
 * The BRANCHEXECLP callback is executed during node processing if a fractional LP solution is available. It should
 * split the current problem into smaller subproblems. Usually, the branching is done in a way such that the current
 * fractional LP solution is no longer feasible in the relaxation of the subproblems.  It is, however, possible to
 * create a child node for which the fractional LP solution is still feasible in the relaxation, for example, by
 * branching on a variable with integral LP value.  In every case, you have to make sure that each subproblem is a
 * proper restriction of the current problem.  Otherwise, you risk to produce an infinite path in the search tree.
 *
 * The user gains access to the branching candidates, i.e., to the fractional variables, and their LP solution values by
 * calling the method SCIPgetLPBranchCands(). Furthermore, SCIP provides two methods for performing the actual
 * branching, namely SCIPbranchVar() and SCIPcreateChild().
 *
 * Given an integral variable \f$x\f$ with fractional LP solution value \f$x^*\f$, the method SCIPbranchVar() creates
 * two child nodes; one contains the bound \f$x \le \lfloor x^* \rfloor\f$ and the other one contains the bound \f$x \ge
 * \lceil x^* \rceil\f$, see the BRANCHEXECLP callback in src/scip/branch_mostinf.c for an example. In addition, if a
 * proven lower objective bound of a created child node is known, like after strong branching has been applied, the user
 * may call the method SCIPupdateNodeLowerbound() in order to update the child node's lower bound.
 *
 * Please also see the \ref BRANCHEXEC "further information for the three execution methods".
 *
 * @subsection BRANCHEXECEXT
 *
 * The BRANCHEXECEXT callback is executed during node processing if no LP solution is available and the list of
 * external branching candidates is not empty. It should split the current problem into smaller subproblems. If you
 * do not use relaxation handlers or constraints handlers that provide external branching candidates, you do not need to
 * implement this callback.
 *
 * In contrast to the LP branching candidates and the pseudo branching candidates, the list of external branching
 * candidates will not be generated automatically. The user has to add all variables to the list by calling
 * SCIPaddExternBranchCand() for each of them. Usually, this will happen in the execution method of a relaxation handler or in the
 * enforcement methods of a constraint handler.
 *
 * The user gains access to these branching candidates by calling the method SCIPgetExternBranchCands(). Furthermore,
 * SCIP provides two methods for performing the actual branching with a given solution value, namely SCIPbranchVarVal()
 * and SCIPcreateChild(). SCIPbranchVarVal() allows users to specify the branching point for a variable in contrast to
 * SCIPbranchVar(), which will always use the current LP or pseudo solution.
 *
 * This paragraph contains additional information regarding how the method SCIPbranchVarVal() works. For external branching candidates,
 * there are three principle possibilities:
 * - Given a continuous variable \f$x\f$ with solution value \f$x^*\f$, the method SCIPbranchVarVal() creates
 *   two child nodes; one contains the bound \f$x \le x^* \f$ and the other one contains the bound \f$x \ge x^* \f$.
 * - Given an integer variable \f$x\f$ with fractional solution value \f$x^*\f$, the method
 *   SCIPbranchVarVal() creates two child nodes; one contains the bound \f$x \le \lfloor x^* \rfloor\f$ and the other
 *   one contains the bound \f$x \ge \lceil x^* \rceil\f$.
 * - Given an integer variable \f$x\f$ with integral solution value \f$x^*\f$, the method SCIPbranchVarVal()
 *   creates three child nodes; one contains the bound \f$x \le x^* -1\f$, one contains the bound \f$x \ge x^* +1\f$,
 *   one contains the fixing \f$x = x^*\f$.
 *
 * See the BRANCHEXECEXT callback in src/scip/branch_random.c for an example. In addition, if a proven lower bound of a
 * created child node is known the user may call the method SCIPupdateNodeLowerbound() in order to update the child
 * node's lower bound.
 *
 * Please also see the \ref BRANCHEXEC "further information for the three execution methods".
 *
 * @subsection BRANCHEXECPS
 *
 * The BRANCHEXECPS callback is executed during node processing if no LP solution is available and at least one of the
 * integer variables is not yet fixed. It should split the current problem into smaller subproblems. PS stands for
 * pseudo solution which is the vector of all variables set to their locally best (w.r.t. the objective function)
 * bounds.
 *
 * The user gains access to the branching candidates, i.e., to the non-fixed integer variables, by calling the method
 * SCIPgetPseudoBranchCands(). Furthermore, SCIP provides two methods for performing the actual branching, namely
 * SCIPbranchVar() and SCIPcreateChild().
 *
 * Given an integer variable \f$x\f$ with bounds \f$[l,u]\f$ and not having solved the LP, the method SCIPbranchVar()
 * creates two child nodes:
 * - If both bounds are finite, then the two children will contain the domain reductions \f$x \le x^*\f$, and \f$x \ge
 *   x^*+1\f$ with \f$x^* = \lfloor \frac{l + u}{2}\rfloor\f$. The current pseudo solution will remain feasible in one
 *   of the branches, but the hope is that halving the domain's size leads to good propagations.
 * - If only one of the bounds is finite, the variable will be fixed to that bound in one of the child nodes. In the
 *   other child node, the bound will be shifted by one.
 * - If both bounds are infinite, three children will be created: \f$x \le 1\f$, \f$x \ge 1\f$, and \f$x = 0\f$.

 *
 * See the BRANCHEXECPS callback in src/scip/branch_random.c for an example. In addition, if a proven lower bound of a
 * created child node is known, the user may call the method SCIPupdateNodeLowerbound() in order to update the child
 * node's lower bound.
 *
 * Please also see the \ref BRANCHEXEC "further information for the three execution methods".
 *
 * @subsection BRANCHEXEC Further information for the three execution methods
 *
 * In order to apply more general branching schemes, one should use the method SCIPcreateChild().
 * After having created a child node, the additional restrictions of the child node have to be added with calls to
 * SCIPaddConsNode(), SCIPchgVarLbNode(), or SCIPchgVarUbNode().
 * \n
 * In the method SCIPcreateChild(), the branching rule has to assign two values to the new nodes: a node selection
 * priority for each node and an estimate for the objective value of the best feasible solution contained in the subtree
 * after applying the branching. If the method SCIPbranchVar() is used, these values are automatically assigned. For
 * variable based branching schemes, one might use the methods SCIPcalcNodeselPriority() and the method
 * SCIPcalcChildEstimate().
 *
 * In some cases, the branching rule can tighten the current subproblem instead of producing a branching. For example,
 * strong branching might have proven that rounding up a variable would lead to an infeasible LP relaxation and thus,
 * the variable must be rounded down. Therefore, the BRANCHEXECLP, BRANCHEXECPS and BRANCHEXECREL callbacks may also
 * produce domain reductions or add additional constraints to the current subproblem.
 *
 * The execution callbacks have the following options:
 *  - detecting that the node is infeasible and can be cut off (result SCIP_CUTOFF)
 *  - adding an additional constraint (e.g. a conflict constraint) (result SCIP_CONSADDED; note that this action
 *    must not be performed if the input "allowaddcons" is FALSE)
 *  - reducing the domain of a variable such that the current LP solution becomes infeasible (result SCIP_REDUCEDDOM)
 *  - applying a branching (result SCIP_BRANCHED)
 *  - stating that the branching rule was skipped (result SCIP_DIDNOTRUN).
 *
 * Only the BRANCHEXECLP callback has the possibility to add a cutting plane to the LP (result SCIP_SEPARATED).
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
 * If you have allocated memory for fields in your branching rule data, remember to free this memory
 * before freeing the branching rule data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection BRANCHINIT
 *
 * The BRANCHINIT callback is executed after the problem is transformed.
 * The branching rule may, e.g., use this call to initialize its branching rule data.
 *
 * @subsection BRANCHCOPY
 *
 * The BRANCHCOPY callback is executed when a SCIP instance is copied, e.g. to
 * solve a sub-SCIP. By
 * defining this callback as
 * <code>NULL</code> the user disables the execution of the specified
 * branching rule for all copied SCIP instances. This may deteriorate the performance
 * of primal heuristics using sub-SCIPs.
 *
 * @subsection BRANCHEXIT
 *
 * The BRANCHEXIT callback is executed before the transformed problem is freed.
 * In this method, the branching rule should free all resources that have been allocated for the solving process in
 * BRANCHINIT.
 *
 * @subsection BRANCHINITSOL
 *
 * The BRANCHINITSOL callback is executed when the presolving is finished and the branch-and-bound process is about to
 * begin.
 * The branching rule may use this call to initialize its branch-and-bound specific data.
 *
 * @subsection BRANCHEXITSOL
 *
 * The BRANCHEXITSOL callback is executed before the branch-and-bound process is freed.
 * The branching rule should use this call to clean up its branch-and-bound data.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page NODESEL How to add node selectors
 *
 * Node selectors are used to decide which of the leaves in the current branching tree is selected as next subproblem
 * to be processed. The ordering relation of the tree's leaves for storing them in the leaf priority queue is also
 * defined by the node selectors.
 * \n
 * A complete list of all node selectors contained in this release can be found \ref NODESELECTORS "here".
 *
 * We now explain how users can add their own node selectors.
 * Take the node selector for depth first search (src/scip/nodesel_dfs.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the scip::ObjNodesel wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_NODESEL... callback methods.
 *
 * Additional documentation for the callback methods of a node selector can be found in the file type_nodesel.h.
 *
 * Here is what you have to do to implement a node selector:
 * -# Copy the template files src/scip/nodesel_xyz.c and src/scip/nodesel_xyz.h into files named "nodesel_mynodeselector.c"
 *    and "nodesel_mynodeselector.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "mynodeselector".
 * -# Adjust the properties of the node selector (see \ref NODESEL_PROPERTIES).
 * -# Define the node selector data (see \ref NODESEL_DATA). This is optional.
 * -# Implement the interface methods (see \ref NODESEL_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref NODESEL_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref NODESEL_ADDITIONALCALLBACKS). This is optional.
 *
 *
 * @section NODESEL_PROPERTIES Properties of a Node Selector
 *
 * At the top of the new file "nodesel_mynodeselector.c" you can find the node selector properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the node selector properties by calling the constructor
 * of the abstract base class scip::ObjNodesel from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par NODESEL_NAME: the name of the node selector.
 * This name is used in the interactive shell to address the node selector.
 * Additionally, if you are searching for a node selector with SCIPfindNodesel(), this name is looked up.
 * Names have to be unique: no two node selectors may have the same name.
 *
 * \par NODESEL_DESC: the description of the node selector.
 * This string is printed as a description of the node selector in the interactive shell.
 *
 * \par NODESEL_STDPRIORITY: the default priority of the node selector in the standard mode.
 * The first step of each iteration of the main solving loop is the selection of the next subproblem to be processed.
 * The node selector of highest priority (the active node selector) is called to do this selection.
 * In particular, if you implemented your own node selector plugin which you want to be applied, you should choose a priority
 * which is greater then all priorities of the SCIP default node selectors.
 * Note that SCIP has two different operation modes: the standard mode and the memory saving mode. If the memory
 * limit - given as a parameter by the user - is almost reached, SCIP switches from the standard mode to the memory saving
 * mode in which different priorities for the node selectors are applied. NODESEL_STDPRIORITY is the priority of the
 * node selector used in the standard mode.
 * \n
 * Note that this property only defines the default value of the priority. The user may change this value arbitrarily by
 * adjusting the corresponding parameter setting.
 *
 * \par NODESEL_MEMSAVEPRIORITY: the default priority of the node selector in the memory saving mode.
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
 * At the bottom of "nodesel_mynodeselector.c", you can find the interface method SCIPincludeNodeselMynodeselector(),
 * which also appears in "nodesel_mynodeselector.h"
 * SCIPincludeNodeselMynodeselector() is called by the user, if (s)he wants to include the node selector,
 * i.e., if (s)he wants to use the node selector in his/her application.
 *
 * This method only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the node selector. For this, you can either call
 * SCIPincludeNodesel(),
 * or SCIPincludeNodeselBasic() since SCIP version 3.0. In the latter variant, \ref NODESEL_ADDITIONALCALLBACKS "additional callbacks"
 * must be added via setter functions as, e.g., SCIPsetNodeselCopy(). We recommend this latter variant because
 * it is more stable towards future SCIP versions which might have more callbacks, whereas source code using the first
 * variant must be manually adjusted with every SCIP release containing new callbacks for node selectors in order to compile.
 *
 *
 * If you are using node selector data, you have to allocate the memory for the data at this point.
 * You can do this by calling:
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
 * The fundamental callback methods of the plugins are the ones that have to be implemented in order to obtain
 * an operational algorithm.
 * They are passed together with the node selector itself to SCIP using SCIPincludeNodesel() or SCIPincludeNodeselBasic(),
 * see @ref NODESEL_INTERFACE.
 *
 * Node selector plugins have two fundamental callback methods, namely the NODESELSELECT method and the NODESELCOMP method.
 * These methods have to be implemented for every node selector; the other callback methods are optional.
 * They implement the two requirements every node selector has to fulfill: Selecting a node from the queue to be processed
 * next and, given two nodes, deciding which of both is favored by the node selector's selection rule. The first
 * task is implemented in the NODESELSELECT callback, the second one in the NODESELCOMP callback.
 * In the C++ wrapper class scip::ObjNodesel, the scip_select() method and the scip_comp() method (which correspond to the
 * NODESELSELECT callback and the NODESELCOMP callback, respectively) are virtual abstract member functions.
 * You have to implement them in order to be able to construct an object of your node selector class.
 *
 * Additional documentation for the callback methods can be found in type_nodesel.h.
 *
 * @subsection NODESELSELECT
 *
 * The NODESELSELECT callback is the first method called in each iteration in the main solving loop. It should decide
 * which of the leaves in the current branching tree is selected as the next subproblem to be processed.
 * It can arbitrarily decide between all leaves stored in the tree, but for performance reasons,
 * the current node's children and siblings are often treated different from the remaining leaves.
 * This is mainly due to the warm start capabilities of the simplex algorithm and the expectation that the bases of
 * neighboring vertices in the branching tree very similar.
 * The node selector's choice of the next node to process can
 * have a large impact on the solver's performance, because it influences the finding of feasible solutions and the
 * development of the global dual bound.
 *
 * Besides the ranking of the node selector, every node gets assigned a node selection priority by the branching rule
 * that created the node. See the \ref BRANCHEXECLP and \ref BRANCHEXECPS callbacks of the branching rules for details.
 * For example, the node where the branching went in the same way as the deviation from the branching variable's
 * root solution could be assigned a higher priority than the node where the branching went in the opposite direction.
 *
 * The following methods provide access to the various types of leaf nodes:
 * - SCIPgetPrioChild() returns the child of the current node with the largest node selection priority, as assigned by the
 *   branching rule.
 *   If no child is available (for example, because the current node was pruned), a NULL pointer is returned.
 * - SCIPgetBestChild() returns the best child of the current node with respect to the node selector's ordering relation as
 *   defined by the \ref NODESELCOMP callback. If no child is available, a NULL pointer is returned.
 * - SCIPgetPrioSibling() returns the sibling of the current node with the largest node selection priority.
 *   If no sibling is available (for example, because all siblings of the current node have already been processed), a NULL
 *   pointer is returned.
 *   Note that in binary branching every node has at most one sibling, but since SCIP supports arbitrary branching rules,
 *   this might not always be the case.
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
 *  - value > 0, if node 1 comes after (is worse than) node 2.
 *
 * @section NODESEL_ADDITIONALCALLBACKS Additional Callback Methods of a Node Selector
 *
 * The additional callback methods do not need to be implemented in every case. However, some of them have to be
 * implemented for most applications, they can be used, for example, to initialize and free private data.
 * Additional callbacks can either be passed directly with SCIPincludeNodesel() to SCIP or via specific
 * <b>setter functions</b> after a call of SCIPincludeNodeselBasic(), see also @ref NODESEL_INTERFACE.
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
 * If you have allocated memory for fields in your node selector data, remember to free this memory
 * before freeing the node selector data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection NODESELINIT
 *
 * The NODESELINIT callback is executed after the problem is transformed.
 * The node selector may, e.g., use this call to initialize its node selector data.
 *
 * @subsection NODESELCOPY
 *
 * The NODESELCOPY callback is executed when a SCIP instance is copied, e.g. to
 * solve a sub-SCIP. By
 * defining this callback as
 * <code>NULL</code> the user disables the execution of the specified
 * node selector for all copied SCIP instances. This may deteriorate the performance
 * of primal heuristics using sub-SCIPs.
 *
 * @subsection NODESELEXIT
 *
 * The NODESELEXIT callback is executed before the transformed problem is freed.
 * In this method, the node selector should free all resources that have been allocated for the solving process
 * in NODESELINIT.
 *
 * @subsection NODESELINITSOL
 *
 * The NODESELINITSOL callback is executed when the presolving is finished and the branch-and-bound process is about to
 * begin.
 * The node selector may use this call to initialize its branch-and-bound specific data.
 *
 * @subsection NODESELEXITSOL
 *
 * The NODESELEXITSOL callback is executed before the branch-and-bound process is freed.
 * The node selector should use this call to clean up its branch-and-bound data.
 */


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page HEUR How to add primal heuristics
 *
 * Feasible solutions can be found in two different ways during the traversal of the branch-and-bound tree. On one
 * hand, the solution of a node's relaxation may be feasible with respect to the constraints (including the integrality).
 * On the other hand, feasible solutions can be discovered by primal heuristics.
 * \n
 * A complete list of all primal heuristics contained in this release can be found \ref PRIMALHEURISTICS "here".
 *
 * We now explain how users can add their own primal heuristics.
 * Take the simple and fast LP rounding heuristic (src/scip/heur_simplerounding.c) as an example.
 * The idea of simple rounding is to iterate over all fractional variables of an LP solution and round them down,
 * if the variables appears only with nonnegative coefficients in the system Ax <= b and round them up if
 * the variables appears only with nonpositive coefficients.
 * If one of both conditions applies for each of the fractional variables, this will give a feasible solution.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the scip::ObjHeur wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_HEUR... callback methods.
 *
 * Additional documentation for the callback methods of a primal heuristic can be found in the file type_heur.h.
 *
 * Here is what you have to do to implement a primal heuristic:
 * -# Copy the template files src/scip/heur_xyz.c and src/scip/heur_xyz.h into files named "heur_myheuristic.c"
 *    and "heur_myheuristic.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "myheuristic".
 * -# Adjust the properties of the primal heuristic (see \ref HEUR_PROPERTIES).
 * -# Define the primal heuristic data (see \ref HEUR_DATA). This is optional.
 * -# Implement the interface methods (see \ref HEUR_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref HEUR_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref HEUR_ADDITIONALCALLBACKS). This is optional.
 *
 *
 * @section HEUR_PROPERTIES Properties of a Primal Heuristic
 *
 * At the top of the new file "heur_myheuristic.c" you can find the primal heuristic properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the primal heuristic properties by calling the constructor
 * of the abstract base class scip::ObjHeur from within your constructor.
 * Of course, all of them are of relevant, but the most important ones for controlling the performance
 * are usually HEUR_FREQ and HEUR_TIMING.
 * The properties you have to set have the following meaning:
 *
 * \par HEUR_NAME: the name of the primal heuristic.
 * This name is used in the interactive shell to address the primal heuristic.
 * Additionally, if you are searching for a primal heuristic with SCIPfindHeur(), this name is looked up.
 * Names have to be unique: no two primal heuristics may have the same name.
 *
 * \par HEUR_DESC: the description of the primal heuristic.
 * This string is printed as a description of the primal heuristic in the interactive shell when you call "display heuristics".
 *
 * \par HEUR_DISPCHAR: the display character of the primal heuristic.
 * In the interactive shell, this character is printed in the first column of a status information row, if the primal
 * heuristic found the feasible solution belonging to the primal bound. Note that a star '*' stands for an integral
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
 * a high priority (like simple rounding). In addition, the interaction between different types of primal heuristics should be taken into account.
 * For example, improvement heuristics, which try to generate improved solutions by inspecting one or more of the feasible
 * solutions that have already been found, should have a low priority (like Crossover which by default needs at least 3 feasible solutions).
 *
 * \par HEUR_FREQ: the default frequency for executing the primal heuristic.
 * The frequency together with the frequency offset (see HEUR_FREQOFS) defines the depth levels at which the execution
 * method of the primal heuristic \ref HEUREXEC is called. For example, a frequency of 7 together with a frequency offset
 * of 5 means, that the \ref HEUREXEC callback is executed for subproblems that are in depth 5, 12, 19, ... of the branching tree. A
 * frequency of 0 together with a frequency offset of 3 means, that the execution method is only called at those nodes that are in
 * depth level 3 (i.e., at most for \f$2^3 = 8\f$ nodes if binary branching is applied).
 * Typical cases are: A frequency of 0 and an offset of 0 which means that
 * the heuristic is only called at the root node and a frequency of -1 which disables the heuristic.
 * \n
 * The frequency can be adjusted by the user. This property of the primal heuristic only defines the default value of the
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
 * heuristic is called. Use -1 for no limit (a usual case).
 *
 * \par HEUR_TIMING: the execution timing of the primal heuristic.
 * Primal heuristics have different entry points during the solving process and the execution timing parameter defines the
 * entry point at which the primal heuristic is executed first.
 * \n
 * The primal heuristic can be called first:
 * - before the processing of the node starts (SCIP_HEURTIMING_BEFORENODE)
 * - after each LP solve during the cut-and-price loop (SCIP_HEURTIMING_DURINGLPLOOP)
 * - after the cut-and-price loop was finished (SCIP_HEURTIMING_AFTERLPLOOP)
 * - after the processing of a node <em>with solved LP</em>  was finished (SCIP_HEURTIMING_AFTERLPNODE)
 * - after the processing of a node <em>without solved LP</em> was finished (SCIP_HEURTIMING_AFTERPSEUDONODE)
 * - after the processing of the last node in the current plunge was finished, <em>and only if the LP was solved for
 *   this node</em> (SCIP_HEURTIMING_AFTERLPPLUNGE)
 * - after the processing of the last node in the current plunge was finished, <em>and only if the LP was not solved
 *   for this node</em> (SCIP_HEURTIMING_AFTERPSEUDOPLUNGE).
 * \par
 * A plunge is the successive solving of child and sibling nodes in the search tree.
 * The flags listed above can be combined to call the heuristic at multiple times by concatenating them with a bitwise OR.
 * Two useful combinations are already predefined:
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
 * Very fast primal heuristics that require an LP solution can also be called "after each LP solve during the
 * cut-and-price loop". Rounding heuristics, like the simple and fast LP rounding heuristic
 * (src/scip/heur_simplerounding.c), belong to this group of primal heuristics.
 * \n
 * Most heuristics, however, are called either after a node was completely processed
 * (e.g. expensive rounding heuristics like RENS), or even only after a full plunge was finished (e.g., diving heuristics).
 *
 * \par HEUR_USESSUBSCIP: Does the heuristic use a secondary SCIP instance?
 * Some heuristics and separators solve MIPs or SAT problems using a secondary SCIP instance. Examples are
 * Large Neighborhood Search heuristics such as RINS and Local Branching or the CGMIP separator. To avoid recursion,
 * these plugins usually deactivate all other plugins that solve MIPs. If a heuristic uses a secondary SCIP instance,
 * this parameter has to be TRUE and it is recommended to call SCIPsetSubscipsOff() for the secondary SCIP instance.
 *
 * Computational experiments indicate that for the overall performance of a MIP solver, it is important to evenly
 * spread the application of the heuristics across the branch-and-bound tree. Thus, the assignment of the parameters
 * HEUR_FREQ, HEUR_FREQOFS, and HEUR_TIMING should contribute to this aim.
 *
 * Note that all diving heuristics in the SCIP distribution (see, e.g., src/scip/heur_guideddiving.c) check whether other diving
 * heuristics have already been called at the current node. This can be done by comparing SCIPgetLastDivenode(scip) and
 * SCIPgetNNodes(scip). If the two are equal, and if the current node is not the root node (SCIPgetDepth(scip) > 0), diving
 * heuristics should be delayed by returning the result code 'SCIP_DELAYED'. This is an additional contribution to the goal of
 * not calling multiple similar heuristics at the same node.
 *
 *
 * @section HEUR_DATA Primal Heuristic Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_HeurData".
 * In this data structure, you can store the data of your primal heuristic. For example, you should store the adjustable
 * parameters of the primal heuristic or a working solution in this data structure.
 * If you are using C++, you can add primal heuristic data as usual as object variables to your class.
 * \n
 * Defining primal heuristic data is optional. You can leave the struct empty.
 *
 *
 * @section HEUR_INTERFACE Interface Methods
 *
 * At the bottom of "heur_myheuristic.c", you can find the interface method SCIPincludeHeurMyheuristic(),
 * which also appears in "heur_myheuristic.h"
 * SCIPincludeHeurMyheuristic() is called by the user, if (s)he wants to include the heuristic,
 * i.e., if (s)he wants to use the heuristic in his/her application.
 *
 * This method only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the heuristic. For this, you can either call
 * SCIPincludeHeur(),
 * or SCIPincludeHeurBasic() since SCIP version 3.0. In the latter variant, \ref HEUR_ADDITIONALCALLBACKS "additional callbacks"
 * must be added via setter functions as, e.g., SCIPsetHeurCopy(). We recommend this latter variant because
 * it is more stable towards future SCIP versions which might have more callbacks, whereas source code using the first
 * variant must be manually adjusted with every SCIP release containing new callbacks for heuristics in order to compile.
 *
 * If you are using primal heuristic data, you have to allocate the memory for the data at this point.
 * You can do this by calling:
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_HeurData afterwards.
 *
 * You may also add user parameters for your primal heuristic, see the method SCIPincludeHeurFeaspump() in
 * src/scip/heur_oneopt.c for an example where a single Boolean parameter is added.
 *
 *
 * @section HEUR_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Primal Heuristic
 *
 * The fundamental callback methods of the plugins are the ones that have to be implemented in order to obtain
 * an operational algorithm.
 * They are passed together with the primal heuristic itself to SCIP using SCIPincludeHeur() or SCIPincludeHeurBasic(),
 * see @ref HEUR_INTERFACE.
 *
 *
 * Primal heuristic plugins have only one fundamental callback method, namely the HEUREXEC method.
 * This method has to be implemented for every primal heuristic; the other callback methods are optional.
 * In the C++ wrapper class scip::ObjHeur, the scip_exec() method (which corresponds to the HEUREXEC callback) is a virtual
 * abstract member function. You have to implement it in order to be able to construct an object of your primal heuristic
 * class.
 *
 * Additional documentation for the callback methods can be found in type_heur.h.
 *
 * @subsection HEUREXEC
 *
 * The HEUREXEC callback is called at different positions during the node processing loop, see HEUR_TIMING. It should
 * search for feasible solutions and add them to the solution pool. For creating a new feasible solution, the
 * methods SCIPcreateSol() and SCIPsetSolVal() can be used. Afterwards, the solution can be added to the storage by
 * calling the method SCIPtrySolFree() (or SCIPtrySol() and SCIPfreeSol()).
 *
 * The HEUREXEC callback gets a SCIP pointer, a pointer to the heuristic itself, the current point in the
 * solve loop and a result pointer as input (see type_heur.h).
 *
 * The heuristic has to set the result pointer appropriately!
 * Therefore it has the following options:
 *  - finding at least one feasible solution (result SCIP_FOUNDSOL)
 *  - stating that the primal heuristic searched, but did not find a feasible solution (result SCIP_DIDNOTFIND)
 *  - stating that the primal heuristic was skipped (result SCIP_DIDNOTRUN)
 *  - stating that the primal heuristic was skipped, but should be called again (result SCIP_DELAYED).
 *
 *
 * @section HEUR_ADDITIONALCALLBACKS Additional Callback Methods of a Primal Heuristic
 *
 * The additional callback methods do not need to be implemented in every case. However, some of them have to be
 * implemented for most applications, they can be used, for example, to initialize and free private data.
 * Additional callbacks can either be passed directly with SCIPincludeHeur() to SCIP or via specific
 * <b>setter functions</b> after a call of SCIPincludeHeurBasic(), see also @ref HEUR_INTERFACE.
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
 * If you have allocated memory for fields in your primal heuristic data, remember to free this memory
 * before freeing the primal heuristic data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection HEURINIT
 *
 * The HEURINIT callback is executed after the problem is transformed.
 * The primal heuristic may, e.g., use this call to initialize its primal heuristic data.
 *
 * @subsection HEURCOPY
 *
 * The HEURCOPY callback is executed when a SCIP instance is copied, e.g. to
 * solve a sub-SCIP. By
 * defining this callback as
 * <code>NULL</code> the user disables the execution of the specified
 * heuristic for all copied SCIP instances. This may deteriorate the performance
 * of primal heuristics using sub-SCIPs.
 *
 * @subsection HEUREXIT
 *
 * The HEUREXIT callback is executed before the transformed problem is freed.
 * In this method, the primal heuristic should free all resources that have been allocated for the solving process in
 * HEURINIT.
 *
 * @subsection HEURINITSOL
 *
 * The HEURINITSOL callback is executed when the presolving is finished and the branch-and-bound process is about to
 * begin. The primal heuristic may use this call to initialize its branch-and-bound specific data.
 *
 * @subsection HEUREXITSOL
 *
 * The HEUREXITSOL callback is executed before the branch-and-bound process is freed. The primal heuristic should use this
 * call to clean up its branch-and-bound data, which was allocated in HEURINITSOL.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page RELAX How to add relaxation handlers
 *
 * SCIP provides specific support for LP relaxations of constraint integer programs. In addition, relaxation handlers,
 * also called relaxators, can be used to include other relaxations, e.g. Lagrange relaxations or semidefinite
 * relaxations. The relaxation handler manages the necessary data structures and calls the relaxation solver to generate dual
 * bounds and primal solution candidates.
 * \n
 * However, the data to define a single relaxation must either be extracted by the relaxation handler itself (e.g., from
 * the user defined problem data, the LP information, or the integrality conditions), or be provided by the constraint
 * handlers. In the latter case, the constraint handlers have to be extended to support this specific relaxation.
 * \n
 *
 * We now explain how users can add their own relaxation handlers using the C interface. It is very easy to
 * transfer the C explanation to C++: whenever a method should be implemented using the SCIP_DECL_RELAX... notion,
 * reimplement the corresponding virtual member function of the abstract scip::ObjRelax wrapper base class.
 * Unfortunately, SCIP does not contain a default relaxation handler plugin, which could be used as an example.
 *
 * Additional documentation for the callback methods of a relaxation handler can be found in the file type_relax.h.
 *
 * Here is what you have to do to implement a relaxation handler:
 * -# Copy the template files src/scip/relax_xyz.c and src/scip/relax_xyz.h into files named "relax_myrelaxator.c"
 *    and "relax_myrelaxator.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "myrelaxator".
 * -# Adjust the properties of the relaxation handler (see \ref RELAX_PROPERTIES).
 * -# Define the relaxation handler data (see \ref RELAX_DATA). This is optional.
 * -# Implement the interface methods (see \ref RELAX_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref RELAX_FUNDAMENTALCALLBACKS).
 * -# Implement the additional callback methods (see \ref RELAX_ADDITIONALCALLBACKS). This is optional.
 *
 *
 * @section RELAX_PROPERTIES Properties of a Relaxation Handler
 *
 * At the top of the new file "relax_myrelaxator.c" you can find the relaxation handler properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the relaxation handler properties by calling the constructor
 * of the abstract base class scip::ObjRelax from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par RELAX_NAME: the name of the relaxation handler.
 * This name is used in the interactive shell to address the relaxation handler.
 * Additionally, if you are searching for a relaxation handler with SCIPfindRelax(), this name is looked up.
 * Names have to be unique: no two relaxation handlers may have the same name.
 *
 * \par RELAX_DESC: the description of the relaxation handler.
 * This string is printed as a description of the relaxation handler in the interactive shell.
 *
 * \par RELAX_PRIORITY: the priority of the relaxation handler.
 * During each relaxation solving round, the included relaxation handlers and the
 * price-and-cut loop for solving the LP relaxation are called in a predefined order, which is given by the priorities
 * of the relaxation handlers.
 * First, the relaxation handlers with non-negative priority are called in the order of decreasing priority.
 * Next, the price-and-cut loop for solving the LP relaxation is executed.
 * Finally, the relaxation handlers with negative priority are called in the order of decreasing priority.
 * \n
 * Usually, you will have only one relaxation handler in your application and thus only have to decide whether it should
 * be called before or after solving the LP relaxation. For this decision you should consider the complexity of
 * the relaxation solving algorithm and the impact of the resulting solution: if your relaxation handler provides a fast
 * algorithm that usually has a high impact (i.e. the relaxation is a good approximation of the
 * feasible region of the subproblem and the solution severely improves the dual bound), it should have a non-negative
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
 * At the bottom of "relax_myrelaxator.c", you can find the interface method SCIPincludeRelaxMyrelaxator(),
 * which also appears in "relax_myrelaxator.h".
 * SCIPincludeRelaxMyrelaxator() is called by the user, if (s)he wants to include the relaxation handler,
 * i.e., if (s)he wants to use the relaxation handler in his/her application.
 *
 * This method only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the relaxation handler. For this, you can either call
 * SCIPincludeRelax(),
 * or SCIPincludeRelaxBasic() since SCIP version 3.0. In the latter variant, \ref RELAX_ADDITIONALCALLBACKS "additional callbacks"
 * must be added via setter functions as, e.g., SCIPsetRelaxCopy(). We recommend this latter variant because
 * it is more stable towards future SCIP versions which might have more callbacks, whereas source code using the first
 * variant must be manually adjusted with every SCIP release containing new callbacks for relaxation handlers in order to compile.
 *
 * If you are using relaxation handler data, you have to allocate the memory for the data at this point.
 * You can do this by calling:
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_RelaxData afterwards.
 *
 * You may also add user parameters for your relaxation handler, see the method SCIPincludeConshdlrKnapsack() in
 * the \ref cons_knapsack.h "knapsack constraint handler" for an example of how to add user parameters.
 *
 *
 * @section RELAX_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Relaxation Handler
 *
 * The fundamental callback methods of the plugins are the ones that have to be implemented in order to obtain
 * an operational algorithm.
 * They are passed together with the relaxation handler itself to SCIP using SCIPincludeRelax() or SCIPincludeRelaxBasic(),
 * see @ref RELAX_INTERFACE.
 *
 *
 * Relaxation handler plugins have only one fundamental callback method, namely the \ref RELAXEXEC method.
 * This method has to be implemented for every relaxation handler; the other callback methods are optional.
 * In the C++ wrapper class scip::ObjRelax, the scip_exec() method (which corresponds to the \ref RELAXEXEC callback) is a virtual
 * abstract member function.
 * You have to implement it in order to be able to construct an object of your relaxation handler class.
 *
 * Additional documentation for the callback methods can be found in type_relax.h.
 *
 * @subsection RELAXEXEC
 * The RELAXEXEC is called in each relaxation solving round. It should solve the current
 * subproblem's relaxation.
 *
 * Note that, like the LP relaxation, the relaxation handler should only operate on variables for which the corresponding
 * column exists in the transformed problem. Typical methods called by a relaxation handler are SCIPconstructLP() and SCIPflushLP() to
 * make sure that the LP of the current node is constructed and its data can be accessed via calls to SCIPgetLPRowsData()
 * and SCIPgetLPColsData(), SCIPseparateSol() to call the cutting plane separators for a given primal solution, and
 * SCIPupdateLocalLowerbound() to update the current node's dual bound after having solved the relaxation.
 * In addition, you may want to call SCIPtrySolFree() if you think that you have found a feasible primal solution.
 *
 * The primal solution of the relaxation can be stored inside the data structures of SCIP with
 * <code>SCIPsetRelaxSolVal()</code> and <code>SCIPsetRelaxSolVals()</code> and later accessed by
 * <code>SCIPgetRelaxSolVal()</code>.
 * Furthermore, there is a list of external branching candidates, that can be filled by relaxation handlers and constraint handlers,
 * allowing branching rules to take these candidates as a guide on how to split the problem into subproblems.
 * Relaxation handlers should store appropriate candidates in this list using the method <code>SCIPaddExternBranchCand()</code>.
 *
 * Usually, the RELAXEXEC callback only solves the relaxation and provides a lower (dual) bound with a call to
 * SCIPupdateLocalLowerbound().
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
 *  - interrupting the solving process to wait for additional input, e.g., cutting planes (result SCIP_SUSPENDED)
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
 * The additional callback methods do not need to be implemented in every case. However, some of them have to be
 * implemented for most applications, they can be used, for example, to initialize and free private data.
 * Additional callbacks can either be passed directly with SCIPincludeRelax() to SCIP or via specific
 * <b>setter functions</b> after a call of SCIPincludeRelaxBasic(), see also @ref RELAX_INTERFACE.
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
 * If you have allocated memory for fields in your relaxation handler data, remember to free this memory
 * before freeing the relaxation handler data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection RELAXINIT
 *
 * The RELAXINIT callback is executed after the problem is transformed.
 * The relaxation handler may, e.g., use this call to initialize its relaxation handler data.
 *
 * @subsection RELAXCOPY
 *
 * The RELAXCOPY callback is executed when a SCIP instance is copied, e.g. to
 * solve a sub-SCIP. By
 * defining this callback as
 * <code>NULL</code> the user disables the execution of the specified
 * relaxation handler for all copied SCIP instances. This may deteriorate the performance
 * of primal heuristics using sub-SCIPs.
 *
 * @subsection RELAXEXIT
 *
 * The RELAXEXIT callback is executed before the transformed problem is freed.
 * In this method, the relaxation handler should free all resources that have been allocated for the solving process in
 * RELAXINIT.
 *
 * @subsection RELAXINITSOL
 *
 * The RELAXINITSOL callback is executed when the presolving is finished and the branch-and-bound process is about to
 * begin. The relaxation handler may use this call to initialize its branch-and-bound specific data.
 *
 * @subsection REALXEXITSOL
 *
 * The RELAXEXITSOL callback is executed before the branch-and-bound process is freed.
 * The relaxation handler should use this call to clean up its branch-and-bound data.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page READER How to add file readers
 *
 * Mainly, file readers are called to parse an input file and generate a constraint integer programming model. They
 * create constraints and variables and activate variable pricers if necessary. However, they can also be called, for
 * example, to parse an input file containing information about a primal solution or fixing of variables. Besides that
 * it is possible to use some of them for writing (exporting) the problem in a specific format.  \n A complete list of
 * all file readers contained in this release can be found \ref FILEREADERS "here".
 *
 * Since a file reader is also responsible for writing a file, the user may
 * ask why the readers have not the name "filehandler". This name would
 * represent this plugin much better than the used one.
 * \n
 * The used name "readers" is historically grown. In the beginning of SCIP
 * there was no need to write/export problems. Therefore, the the plugin
 * name "readers" was best fitting for this plugin since only reading was essential.
 * It turned out, however, that it is quite nice to write/export certain subproblem during
 * the solving process mainly for debugging. Therefore, a writing callback
 * was added to the "readers" plugin.
 *
 * We now explain how users can add their own file readers.
 * Take the file reader for MIPs in IBM's Mathematical Programming System format (src/scip/reader_mps.c) as an example.
 * As all other default plugins, it is written in C. C++ users can easily adapt the code by using the scip::ObjReader wrapper
 * base class and implement the scip_...() virtual methods instead of the SCIP_DECL_READER... callback methods.
 *
 * Additional documentation for the callback methods of a file reader can be found in the file type_reader.h.
 *
 * Here is what you have to do to implement a file reader named "myreader" in C:
 * -# Copy the template files src/scip/reader_xyz.c and src/scip/reader_xyz.h into files named
 *    "reader_myreader.c" and "reader_myreader.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "myreader".
 * -# Adjust the \ref READER_PROPERTIES "properties of the file reader".
 * -# Define the \ref READER_DATA "file reader data". This is optional.
 * -# Implement the \ref READER_INTERFACE "interface methods".
 * -# Implement the \ref READER_FUNDAMENTALCALLBACKS "fundamental callback methods".
 * -# Implement the \ref READER_ADDITIONALCALLBACKS "additional callback methods". This is optional.
 *
 *
 * @section READER_PROPERTIES Properties of a File Reader
 *
 * At the top of the new file "reader_myreader.c" you can find the file reader properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the file reader properties by calling the constructor
 * of the abstract base class scip::ObjReader from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par READER_NAME: the name of the file reader.
 * This name is used in the interactive shell to address the file reader.
 * Additionally, if you are searching for a file reader with SCIPfindReader(), this name is looked up.
 * Names have to be unique: no two file readers may have the same name.
 *
 * \par READER_DESC: the description of the file reader.
 * This string is printed as a description of the file reader in the interactive shell.
 *
 * \par READER_EXTENSION: the file name extension of the file reader.
 * Each file reader is hooked to a single file name extension. It is automatically called if the user wants to read in a
 * file of corresponding name. The extensions of the different file readers have to be unique.
 * Note that the additional extension '.gz', '.z', or '.Z' (indicating a gzip compressed file) are ignored for assigning
 * an input file to a reader.
 * \n
 * It is not possible to hook up a (single) file reader with more than one file extension.
 * It is, however, not necessary to implement the same (parsing/writing) methods more than once, if you want to
 * support several file extension with the same parser. To do so look at the files reader_lp.c
 * and reader_rlp.c. Both support the LP format.
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
 * At the bottom of "reader_myreader.c", you can find the interface method SCIPincludeReaderMyreader(),
 * which also appears in "reader_myreader.h".
 * SCIPincludeReaderMyreader() is called by the user, if (s)he wants to include the reader,
 * i.e., if (s)he wants to use the reader in his/her application.
 *
 * This method only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the reader. For this, you can either call
 * SCIPincludeReader(),
 * or SCIPincludeReaderBasic() since SCIP version 3.0. In the latter variant, \ref READER_ADDITIONALCALLBACKS "additional callbacks"
 * must be added via setter functions as, e.g., SCIPsetReaderCopy(). We recommend this latter variant because
 * it is more stable towards future SCIP versions which might have more callbacks, whereas source code using the first
 * variant must be manually adjusted with every SCIP release containing new callbacks for readers in order to compile.
 *
 * If you are using file reader data, you have to allocate the memory for the data at this point.
 * You can do this by calling:
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
 * File reader plugins have no fundamental callback methods. This is due to
 * the fact that a file reader can be used for reading and/or writing a
 * file. A file reader is only useful if the reader method \ref READERREAD
 * and/or the writing method \ref READERWRITE is implemented.  One of these
 * methods should be implemented for every file reader; the other callback
 * methods \ref READERCOPY and \ref READERFREE are optional.  In the C++ wrapper class scip::ObjReader, the
 * scip_read() and scip_write() methods (which corresponds to the \ref
 * READERREAD and \ref READERWRITE callback) are virtual member
 * functions. At least one of them should be implemented.
 *
 * Additional documentation for the callback methods can be found in type_reader.h.
 *
 *
 * @section READER_ADDITIONALCALLBACKS Additional Callback Methods of a File Reader
 *
 * Additional callbacks can either be passed directly with SCIPincludeReader() to SCIP or via specific
 * <b>setter functions</b> after a call of SCIPincludeReaderBasic(), see also @ref READER_INTERFACE.
 *
 *
 * File reader plugins contain only additional callback methods, namely the methods \ref READERREAD,
 * \ref READERWRITE, \ref READERFREE, and \ref READERCOPY. Therefore, these are not needed to be implemented. However,
 * at least \ref READERREAD and/or \ref READERWRITE should be implemented (see notes
 * \ref READER_FUNDAMENTALCALLBACKS "above").
 *
 *
 * @subsection READERREAD
 *
 * The READERREAD callback is called when the user invokes SCIP to read in a file with file name extension
 * corresponding to the READER_EXTENSION property of the file reader. This is usually triggered by a call to the method
 * SCIPreadProb() or by an interactive shell command.
 * The READERREAD callback should parse the input file and perform the desired action, which usually means
 * generating a constraint integer programming model, adding a primal solution, fixing variables
 * in an existing model.
 * \n
 * Typical methods called by a file reader that is used to read/generate constraint
 * integer programming models are, for example,
 *
 * - creating an empty problem: SCIPcreateProb()
 * - creating the variables: SCIPcreateVar(), SCIPchgVarType(), SCIPchgVarLb(), SCIPchgVarUb(), SCIPaddVar(), and
 *   SCIPreleaseVar()
 * - modifying the objective function: SCIPchgVarObj() and SCIPsetObjsense().
 * - creating the constraints: SCIPcreateConsLinear(), SCIPaddCoefLinear(), SCIPchgLhsLinear(), SCIPchgRhsLinear(),
 *   SCIPaddCons(), and SCIPreleaseCons()
 *
 * Primal solutions can only be created for the transformed problem. Therefore, the user has to call SCIPtransformProb()
 * before (s)he reads in the file containing the solution and adds it to the solution pool via the method SCIPreadSol().
 *
 *
 * @subsection READERWRITE
 *
 * The READERWRITE callback is called when the user invokes SCIP to write a problem (original or transformed)
 * in the format the reader supports. This is only possible if this callback is implemented. To write the problem
 * all necessary information is given through the parameters of this callback method (see type_reader.h). This
 * information should be used to output the problem in the requested format. This callback method is usually
 * triggered by the call of the methods SCIPwriteOrigProblem(), SCIPwriteTransProblem(), SCIPprintOrigProblem(),
 * or SCIPprintTransProblem().
 * \n
 * A typical method called by a file reader which is used to write/export a constraint
 * integer programming model is SCIPinfoMessage(). This method outputs a given string into a file
 * or into stdout.
 * \n
 * For an example we refer to the writing method of the MPS reader (see reader_mps.c).
 *
 *
 * @subsection READERCOPY
 *
 * The READERCOPY callback is executed when a SCIP instance is copied, e.g. to solve a sub-SCIP. By defining this
 * callback as <code>NULL</code> the user disables the execution of the specified reader for all copied SCIP
 * instances. The question might arise why to copy that plugin. In case of debugging it is nice to be able to
 * write/display the copied instances. Since the reader is in charge of that, you might want to copy the plugin. Below
 * you see a standard implementation.
 *
 * \code
 * static
 * SCIP_DECL_READERCOPY(readerCopyMyreader)
 * {
 *    assert(scip != NULL);
 *    assert(reader != NULL);
 *    assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
 *
 *    SCIP_CALL( SCIPincludeReaderMyreader(scip) );
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
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
 * If you have allocated memory for fields in your file reader data, remember to free this memory
 * before freeing the file reader data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DIALOG How to add dialogs
 *
 * SCIP comes with a command line shell which allows the user to read in problem instances, modify the solver's
 * parameters, initiate the optimization and display certain statistics and solution information. This shell consists
 * of dialogs, which are organized as a tree in SCIP. A node of this tree which is not a leaf represents a menu in
 * the shell and the children of this node correspond to the entries of this menu (which can again be menus). All
 * different dialogs are managed by a dialog handler, which, in particular, is responsible for executing the dialog
 * corresponding to the user's command in the shell. The concept of a dialog handler is different to that
 * of a constraint handler, which is used to manage objects of the same structure, see \ref CONS. In particular, SCIP
 * features only one dialog handler (dialog_default.h), whereas there may exist different constraint handlers.
 * \n
 * A complete list of all dialogs contained in this release can be found \ref DIALOGS "here".
 *
 * We now explain how users can extend the interactive shell by adding their own dialog.
 * We give the explanation for creating your own source file for each additional dialog. Of course, you can collect
 * different dialogs in one source file. Take src/scip/dialog_default.c, where all default dialog plugins are collected, as an
 * example.
 * As all other default plugins, the default dialog plugin and the template dialog are written in C. C++ users can easily
 * adapt the code by using the scip::ObjDialog wrapper base class and implement the scip_...() virtual methods instead of the
 * SCIP_DECL_DIALOG... callback methods.
 *
 * Additional documentation for the callback methods of a dialog can be found in the file type_dialog.h.
 *
 * Here is what you have to do to add a dialog (assuming your dialog is named "mydialog"):
 * -# Copy the template files src/scip/dialog_xyz.c and src/scip/dialog_xyz.h into files named "dialog_mydialog.c"
 *    and "dialog_mydialog.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "mydialog".
 * -# Adjust the \ref DIALOG_PROPERTIES "properties of the dialog".
 * -# Define the \ref DIALOG_DATA "dialog data". This is optional.
 * -# Implement the \ref DIALOG_INTERFACE "interface methods".
 * -# Implement the \ref DIALOG_FUNDAMENTALCALLBACKS "fundamental callback methods".
 * -# Implement the \ref DIALOG_ADDITIONALCALLBACKS "additional callback methods". This is optional.
 *
 *
 * @section DIALOG_PROPERTIES Properties of a Dialog
 *
 * At the top of the new file "dialog_mydialog.c" you can find the dialog properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the dialog properties by calling the constructor
 * of the abstract base class scip::ObjDialog from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par DIALOG_NAME: the name of the dialog.
 * In the interactive shell, this name appears as the command name of the dialog in the parent dialog.
 * Additionally, if you are searching an entry in a menu with SCIPdialogFindEntry(), this name is looked up.
 * Names within one menu have to be unique: no two dialogs in the same menu may have the same name.
 *
 * \par DIALOG_DESC: the description of the dialog.
 * This string is printed as a description of the dialog in the interactive shell if the additional
 * callback method \ref DIALOGDESC is not implemented.
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
 * This method only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the dialog, which can be done by the following lines of code:
 * \code
 * if( !SCIPdialogHasEntry(parentdialog, DIALOG_NAME) )
 * {
 *    SCIP_CALL( SCIPcreateDialog(scip, &dialog, dialogExecMydialog, dialogDescMydialog, dialogFreeMydialog,
 *          DIALOG_NAME, DIALOG_DESC, DIALOG_ISSUBMENU, dialogdata) );
 *
 *    SCIP_CALL( SCIPaddDialogEntry(scip, parentdialog, dialog) );
 *
 *    SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
 * }
 * \endcode
 * Here "parentdialog" has to be an existing dialog which is defined to be a menu (see DIALOG_ISSUBMENU), e.g.,
 * the default root dialog. The method SCIPgetRootDialog() returns the root dialog.
 *
 * The interface method is called by the user, if (s)he wants to include the dialog, i.e., if (s)he wants to use the dialog in
 * his/her application.
 * Note that in order to be able to link the new dialog to an existing default dialog
 * (except the root dialog) it has to be included <b>after the
 * default dialogs plugin</b>, i.e., the SCIPincludeDialogMydialog() call has to occur after the
 * SCIPincludeDialogDefault() call. The SCIPincludeDialogDefault() method is called from within the SCIPincludeDefaultPlugins()
 * method. Therefore, it suffices to include your dialog plugins after you have called SCIPincludeDefaultPlugins().
 * In case you want to add a dialog to the <b>root dialog</b>, you just use the following
 * lines of code to get/create the root dialog.
 *
 * \code
 * SCIP_DIALOG* root;
 *
 * root = SCIPgetRootDialog(scip);
 * if( root == NULL )
 * {
 *    SCIP_CALL( SCIPcreateRootDialog(scip, &root) );
 * }
 * assert( root != NULL );
 * \endcode
 *
 * Therefore, in this case you do not have to worry about the calls of
 * SCIPincludeDialogDefault() and SCIPincludeDefaultPlugins() .
 *
 * If you are using dialog data, you have to allocate the memory for the data at this point.
 * You can do this by calling:
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &dialogdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_DialogData afterwards.
 *
 * Consider the following example. The user wants to add a "drawgraph" command to the root menu of SCIP.
 * (S)he copies the "dialog_xyz.c" and "dialog_xyz.h" files into files "dialog_drawgraph.c" and "dialog_drawgraph.h", respectively.
 * Then, (s)he puts the following code into the SCIPincludeDialogDrawgraph() method, compare SCIPincludeDialogDefault() in
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
 *       SCIP_CALL( SCIPcreateRootDialog(scip, &root) );
 *    }
 *    assert( root != NULL );
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
 * In the C++ wrapper class scip::ObjDialog, the scip_exec() method (which corresponds to the \ref DIALOGEXEC callback) is a virtual
 * abstract member function.
 * You have to implement it in order to be able to construct an object of your dialog class.
 *
 * Additional documentation for the callback methods can be found in type_dialog.h.
 *
 * @subsection DIALOGEXEC
 *
 * The DIALOGEXEC method is invoked, if the user selected the dialog's command name in the parent's menu. It should
 * execute what is stated in DIALOG_DESC, e.g., the display constraint handlers dialog should display information about
 * the constraint handlers included in SCIP, see src/scip/dialog_default.c.
 *
 * For typical methods called by the execution method, have a look at src/scip/dialog_default.c.
 *
 * The callback has to return which dialog should be processed next. This can be, for example, the root dialog
 * (SCIPdialoghdlrGetRoot()), the parent dialog (SCIPdialogGetParent()) or NULL, which stands for closing the interactive
 * shell.
 *
 *
 * @section DIALOG_ADDITIONALCALLBACKS Additional Callback Methods of a Dialog
 *
 * The additional callback methods do not need to be implemented in every case.
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
 * If you have allocated memory for fields in your dialog data, remember to free this memory
 * before freeing the dialog data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection DIALOGDESC
 *
 * This method is called when the help menu of the parent is displayed. It should output (usually a single line of)
 * information describing the meaning of the dialog.
 * \n
 * If this callback is not implemented, the description string of the dialog (DIALOG_DESC) is displayed instead.
 *
 * @subsection DIALOGCOPY
 *
 * The DIALOGCOPY callback is executed when a SCIP instance is copied, e.g. to solve a sub-SCIP. By defining this
 * callback as <code>NULL</code> the user disables the execution of this dialog for all copied SCIP instances. In
 * general there is no need to copy any dialog since it is most unlikely to start the interactive shell of the copied
 * instances.
 *
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DISP How to add display columns
 *
 * While solving a constraint integer program, SCIP displays status information in a column-like fashion. The current
 * number of processed branching tree nodes, the solving time, and the relative gap between primal and dual bound are
 * examples of such display columns. There already exists a wide variety of display columns which can be activated or
 * deactivated on demand, see src/scip/disp_default.c. Additionally, the user can implement his/her own display columns
 * in order to track problem or algorithm specific values.
 * \n
 * A complete list of all displays contained in this release can be found \ref DISPLAYS "here".
 *
 * We now explain users can add their own display columns.
 * We give the explanation for creating your own source file for each additional display column. Of course, you can collect
 * different additional display columns in one source file.
 * Take src/scip/disp_default.c, where all default display columns are collected, as an example.
 * As all other default plugins, the default display column plugins and the display column template are written in C.
 * C++ users can easily adapt the code by using the scip::ObjDisp wrapper base class and implement the scip_...() virtual methods
 * instead of the SCIP_DECL_DISP... callback methods.
 *
 *
 * Additional documentation for the callback methods of a display column can be found in the file type_disp.h.
 *
 * Here is what you have to do to implement a display column (assuming your display column is named "mydisplaycolumn"):
 * -# Copy the template files src/scip/disp_xyz.c and src/scip/disp_xyz.h into files named "disp_mydisplaycolumn.c"
 *    and "disp_mydisplaycolumn.h".
      \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "mydisplaycolumn".
 * -# Adjust the \ref DISP_PROPERTIES "properties of the display column".
 * -# Define the  \ref DISP_DATA "display column data". This is optional.
 * -# Implement the \ref DISP_INTERFACE "interface methods".
 * -# Implement the \ref DISP_FUNDAMENTALCALLBACKS "fundamental callback methods".
 * -# Implement the \ref DISP_ADDITIONALCALLBACKS "additional callback methods". This is optional.
 *
 *
 * @section DISP_PROPERTIES Properties of a Display Column
 *
 * At the top of the new file "disp_mydisplaycolumn.c" you can find the display column properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the display column properties by calling the constructor
 * of the abstract base class scip::ObjDisp from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par DISP_NAME: the name of the display column.
 * This name is used in the interactive shell to address the display column.
 * Additionally, if you are searching for a display column with SCIPfindDisp(), this name is looked up.
 * Names have to be unique: no two display columns may have the same name.
 *
 * \par DISP_DESC: the description of the display column.
 * This string is printed as a description of the display column in the interactive shell.
 *
 * \par DISP_HEADER: the header of the display column.
 * This string is printed as the header of the display column in the status information display.
 *
 * \par DISP_WIDTH: the width of the display column.
 * This parameter defines the width (number of characters) of the display column. The value of the parameter has to be
 * greater than or equal to the number of characters in the header string.
 *
 * \par DISP_PRIORITY: the priority of the display column.
 * The total width of status information lines is bounded by the parameter "display width". The display columns actually contained
 * in the status information display are selected in decreasing order of their priority. Furthermore, the user can force
 * columns to be displayed or not to be displayed in the status information display. For that, (s)he has to switch the value
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
 * This method only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the display column by calling the method
 * SCIPincludeDisp().
 *
 * The interface method is called by the user, if (s)he wants to include the display column, i.e., if (s)he wants to use the display column in his
 * application.
 *
 * If you are using display column data, you have to allocate the memory for the data at this point.
 * You can do this by calling:
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &dispdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_DispData afterwards.
 *
 * Although this is very uncommon, you may also add user parameters for your display column, see the method
 * SCIPincludeConshdlrKnapsack() in the \ref cons_knapsack.h "knapsack constraint handler" for an example.
 *
 *
 * @section DISP_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Display Column
 *
 * Display column plugins have only one fundamental callback method, namely the \ref DISPOUTPUT method.
 * This method has to be implemented for every display column; the other callback methods are optional.
 * In the C++ wrapper class scip::ObjDisp, the scip_output() method (which corresponds to the \ref DISPOUTPUT callback) is a virtual
 * abstract member function.
 * You have to implement it in order to be able to construct an object of your display column class.
 *
 * Additional documentation for the callback methods can be found in type_disp.h.
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
 * The additional callback methods do not need to be implemented in every case.
 * They can be used, for example, to initialize and free private data.
 *
 * @subsection DISPCOPY
 *
 * The DISPCOPY callback is executed when a SCIP instance is copied, e.g. to solve a sub-SCIP. By defining this callback
 * as <code>NULL</code> the user disables the execution of the specified column. In general it is probably not needed to
 * implement that callback since the output of the copied instance is usually suppressed. In the other case or for
 * debugging the callback should be implement.
 *
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
 * If you have allocated memory for fields in your display column data, remember to free this memory
 * before freeing the display column data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 * @subsection DISPINIT
 *
 * The DISPINIT callback is executed after the problem is transformed.
 * The display column may, e.g., use this call to initialize its display column data.
 *
 * @subsection DISPEXIT
 *
 * The DISPEXIT callback is executed before the transformed problem is freed.
 * In this method, the display column should free all resources that have been allocated for the solving process in
 * \ref DISPINIT.
 *
 * @subsection DISPINITSOL
 *
 * The DISPINITSOL callback is executed when the presolving is finished and the branch-and-bound process is about to
 * begin. The display column may use this call to initialize its branch-and-bound specific data.
 *
 * @subsection DISPEXITSOL
 *
 * The DISPEXITSOL callback is executed before the branch-and-bound process is freed. The display column should use this
 * call to clean up its branch-and-bound data specific data.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page EVENT How to add event handler
 *
 * While solving a constraint integer program, SCIP drops thousands of events such as SCIP_EVENTTYPE_VARFIXED (a
 * complete list of all events is given in type_event.h). These events can be caught and used to do something after a
 * certain event happens. Events can be used to speed up the solution process. For example, the set partitioning
 * constraint is only worth propagating if one of the involved variables is fixed. This can be detected by
 * catching the event SCIP_EVENTTYPE_VARFIXED. To be able to catch an event it is necessary to write an event handler
 * which defines what to do after a certain event was caught.
 *
 * We now explain how users can add their own event handlers. We give the explanation for creating your own
 * source file for each additional event handler. Of course, you can collect different event handlers in one source file
 * or you can put the event handler directly into the constraint handler.  In a \ref EVENTUSAGE "second step" we discuss
 * the usage of an event handler. This means how to catch and drop events. \ref EVENTTYPES "Finally", we give some notes on the existing
 * types of events.
 *
 * Take src/scip/cons_logior.c, where the event handler is directly included into the constraint handler. As all other
 * default plugins, the event handlers are written in C. C++ users can easily adapt the code by using the scip::ObjEventhdlr
 * wrapper base class and implement the scip_...() virtual methods instead of the SCIP_DECL_EVENT... callback methods.
 *
 * Additional documentation for the callback methods of an event handler can be found in the file type_event.h. There is
 * also an example written in C which deals with an event handler. You find this example in the directory
 * "examples/Eventhdlr/". An C++ example can be found within the TSP project (examples/TSP/src/EventhdlrNewSol.cpp).
 *
 * Here is what you have to do to implement an event handler (assuming your event handler is named "bestsol"):
 * -# Copy the template files src/scip/event_xyz.c and src/scip/event_xyz.h into files named "event_bestsol.c"
 *    and "event_bestsol.h".
      \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "bestsol".
 * -# Adjust the \ref EVENTHDLR_PROPERTIES "properties of the event handler".
 * -# Implement the \ref EVENT_INTERFACE "interface methods".
 * -# Implement the \ref EVENT_FUNDAMENTALCALLBACKS "fundamental callback methods".
 * -# Implement the \ref EVENT_ADDITIONALCALLBACKS "additional callback methods". This is optional.
 *
 *
 * @section EVENTHDLR_PROPERTIES Properties of a Event Handler
 *
 * At the top of the new file "event_bestsol.c" you can find the event handler properties.
 * These are given as compiler defines.
 * In the C++ wrapper class, you have to provide the event handler properties by calling the constructor
 * of the abstract base class scip::ObjEventhdlr from within your constructor.
 * The properties you have to set have the following meaning:
 *
 * \par EVENT_NAME: the name of the event handler.
 * This name has to be unique with respect to all other event handlers. If you are searching for an event handler with
 * SCIPfindEventhdlr(), this name is looked up.
 *
 * \par EVENT_DESC: the description of the event handler.
 * This string is printed as a description of the event handler.
 *
 * @section EVENTHDLR_DATA Event Handler Data
 *
 * Below the header "Data structures" you can find a struct which is called "struct SCIP_EventhdlrData".
 * In this data structure, you can store the data of your event handler. For example, you should store the adjustable
 * parameters of the event handler in this data structure.
 * If you are using C++, you can add event handler data as usual as object variables to your class.
 * \n
 * Defining event handler data is optional. You can leave the struct empty.
 *
 *
 * @section EVENT_INTERFACE Interface Methods
 *
 * At the bottom of "event_bestsol.c", you can find the interface method SCIPincludeEventBestsol(),
 * which also appears in "event_bestsol.h".
 * SCIPincludeEventBestsol() is called by the user, if (s)he wants to include the event handler,
 * i.e., if (s)he wants to use the event handler in his/her application.
 *
 * This method only has to be adjusted slightly.
 * It is responsible for notifying SCIP of the presence of the event handler. For this, you can either call
 * SCIPincludeEventhdlr(),
 * or SCIPincludeEventhdlrBasic() since SCIP version 3.0. In the latter variant, \ref EVENT_ADDITIONALCALLBACKS "additional callbacks"
 * must be added via setter functions as, e.g., SCIPsetReaderCopy(). We recommend this latter variant because
 * it is more stable towards future SCIP versions which might have more callbacks, whereas source code using the first
 * variant must be manually adjusted with every SCIP release containing new callbacks for event handlers in order to compile.
 *
 * If you are using event handler data, you have to allocate the memory for the data at this point.
 * You can do this by calling:
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_EventhdlrData afterwards.
 *
 * Although this is very uncommon, you may also add user parameters for your event handler, see the method
 * SCIPincludeConshdlrKnapsack() in the \ref cons_knapsack.h "knapsack constraint handler" for an example.
 *
 *
 * @section EVENT_FUNDAMENTALCALLBACKS Fundamental Callback Methods of a Event Handler
 *
 * The fundamental callback methods of the plugins are the ones that have to be implemented in order to obtain
 * an operational algorithm.
 * They are passed together with the event handler itself to SCIP using SCIPincludeEventhdlr() or SCIPincludeEventhdlrBasic(),
 * see @ref EVENT_INTERFACE.
 *
 *
 * Event handler plugins have only one fundamental callback method, namely the \ref EVENTEXEC method.  This method has
 * to be implemented for every event handler; the other callback methods are optional.  In the C++ wrapper class
 * scip::ObjEventhdlr, the scip_exec() method (which corresponds to the \ref EVENTEXEC callback) is a virtual abstract member
 * function.  You have to implement it in order to be able to construct an object of your event handler class.
 *
 * Additional documentation for the callback methods can be found in type_event.h.
 *
 * @subsection EVENTEXEC
 *
 * The EVENTEXEC callback is called after the requested event happened. Then the event handler can do some action in
 * reaction to the event.
 *
 * Typical the execution method sets a parameter to TRUE to indicate later in solving process that something happened
 * which should be analyzed further. In the \ref cons_knapsack.h "knapsack constraint handler" you find such a typical
 * example.
 *
 * @section EVENT_ADDITIONALCALLBACKS Additional Callback Methods of a Event Handler
 *
 * The additional callback methods do not need to be implemented in every case. However, some of them have to be
 * implemented for most applications, they can be used, for example, to initialize and free private data.
 * Additional callbacks can either be passed directly with SCIPincludeEventhdlr() to SCIP or via specific
 * <b>setter functions</b> after a call of SCIPincludeEventhdlrBasic(), see also @ref EVENT_INTERFACE.
 *
 * @subsection EVENTCOPY
 *
 * The EVENTCOPY callback is executed when a SCIP instance is copied, e.g. to solve a sub-SCIP. By defining this
 * callback as <code>NULL</code> the user disables the execution of the specified event handler for all copied SCIP
 * instances. Note that in most cases the event handler in the copied instance will be initialize by those objects (such
 * as constraint handlers or propagators) which need this event handler (see \ref cons_knapsack.h). In these cases the copy
 * callback can be ignored. In case of general events, such as a new best solution being found
 * (SCIP_EVENTTYPE_BESTSOLFOUND), you might want to implement that callback. The event handler example which you find
 * in the directory "examples/Eventhdlr/" uses that callback.
 *
 * \code
 * static
 * SCIP_DECL_EVENTCOPY(eventCopyBestsol)
 * {
 *    assert(scip != NULL);
 *    assert(eventhdlr != NULL);
 *    assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
 *
 *    SCIP_CALL( SCIPincludeEventHdlrBestsol(scip) );
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 *
 *
 * @subsection EVENTFREE
 *
 * If you are using event handler data, you have to implement this method in order to free the event handler data.
 * This can be done by the following procedure:
 * \code
 * static
 * SCIP_DECL_EVENTFREE(eventFreeBestsol)
 * {
 *    SCIP_EVENTHDLRDATA* eventhdlrdata;
 *
 *    eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
 *    assert(eventhdlrdata != NULL);
 *
 *    SCIPfreeMemory(scip, &eventhdlrdata);
 *
 *    SCIPeventhdlrSetData(eventhdlr, NULL);
 *
 *    return SCIP_OKAY;
 * }
 * \endcode
 * If you have allocated memory for fields in your event handler data, remember to free this memory
 * before freeing the event handler data itself.
 * If you are using the C++ wrapper class, this method is not available.
 * Instead, just use the destructor of your class to free the member variables of your class.
 *
 *
 * @subsection EVENTINIT
 *
 * The EVENTINIT callback is executed after the problem is transformed.
 * The event handler may, e.g., use this call to initialize its event handler data.
 *
 * @subsection EVENTEXIT
 *
 * The EVENTEXIT callback is executed before the transformed problem is freed.
 * In this method, the event handler should free all resources that have been allocated for the solving process in
 * \ref EVENTINIT.
 *
 * @subsection EVENTINITSOL
 *
 * The EVENTINITSOL callback is executed when the presolving is finished and the branch-and-bound process is about to
 * begin. The event handler may use this call to initialize its branch-and-bound specific data.
 *
 * @subsection EVENTEXITSOL
 *
 * The EVENTEXITSOL callback is executed before the branch-and-bound process is freed. The event handler should use this
 * call to clean up its branch-and-bound data specific data.
 *
 * @section EVENTUSAGE Catching and Dropping Events
 *
 * After you have implemented the event handler, you have to tell SCIP for which events this event handler should be
 * used. This can be a general events, such as <code>SCIP_EVENTTYPE_BESTSOLFOUND</code>, or a variable event which is the most common
 * way.
 *
 * In case of a general (not variable) event you use the function SCIPcatchEvent() to attach to an event and
 * SCIPdropEvent() to release this event later.
 *
 * \code
 * SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, NULL) );
 * \endcode
 *
 * \code
 * SCIP_CALL( SCIPdropEvent( scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, NULL) );
 * \endcode
 *
 * If you want trigger some variable event, you use the method SCIPcatchVarEvent() to attach the variable event and
 * SCIPdropVarEvent() to drop it later.
 *
 * \code
 * SCIP_CALL( SCIPcatchVarEvent( scip, var, SCIP_EVENTTYPE_VARFIXED, eventhdlr, NULL, NULL) );
 * \endcode
 *
 * \code
 * SCIP_CALL( SCIPdropVarEvent( scip, var, SCIP_EVENTTYPE_VARFIXED, eventhdlr, NULL, NULL) );
 * \endcode
 *
 * @section EVENTTYPES Event types
 *
 * All available events are listed in type_event.h. There are atomic events such as <code>SCIP_EVENTTYPE_VARFIXED</code>
 * and combined events such as <code>SCIP_EVENTTYPE_VARCHANGED</code>. The events are encoded via bit masks. Each atomic
 * event has a unique power of two. This enables combination of the atomic events.
 *
 * SCIP only throws atomic events. However, an event handler might be interested in bunch of events. Through the
 * underlying bit masks it is possible to combine the atomic events. For example, <code>SCIP_EVENTTYPE_VARCHANGED</code>
 * is an event which combines the events <code>SCIP_EVENTTYPE_VARFIXED</code>, <code>SCIP_EVENTTYPE_VARUNLOCKED</code>,
 * <code>SCIP_EVENTTYPE_OBJCHANGED</code>, <code>SCIP_EVENTTYPE_GBDCHANGED</code>,
 * <code>SCIP_EVENTTYPE_DOMCHANGED</code>, and <code>SCIP_EVENTTYPE_IMPLADDED</code>.
 *
 * \code
 * #define SCIP_EVENTTYPE_VARCHANGED     (SCIP_EVENTTYPE_VARFIXED | SCIP_EVENTTYPE_VARUNLOCKED | SCIP_EVENTTYPE_OBJCHANGED
 *                                    | SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_DOMCHANGED | SCIP_EVENTTYPE_IMPLADDED)
 * \endcode
 *
 * Depending on the event type, the event offers different information. The methods which can be used to gain
 * access to this information are given in pub_event.h.
 *
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page NLPI How to add interfaces to nonlinear programming solvers
 *
 * NLPIs are used to interface a solver for nonlinear programs (NLP).
 * It is used, e.g., to solve convex relaxations of the problem or to find locally optimal solutions of
 * nonlinear relaxations or subproblems.
 * The NLPI has been designed such that it can be used independently from SCIP.
 *
 * While the NLPI itself corresponds to the solver interface, the NLPIPROBLEM corresponds to the
 * (solver specific) representation of a concrete nonlinear program.
 * An NLP is specified as a set of indexed variables with variable bounds, an objective function,
 * and a set of constraints, where each constraint is specified as a function which is restricted to lie
 * between given left and right hand sides (possibly infinite).
 * A function consists of a linear, quadratic, and general nonlinear part.
 * The linear and quadratic parts are specified via variable indices and coefficients, while the
 * general nonlinear part is specified via an expression tree.
 * That is, the user of the NLPI does not provide function evaluation callbacks but an algebraic representation of the NLP.
 * Interfaces for solvers that require function evaluations can make use of the NLPIORACLE, which
 * provides a set of methods to compute functions values, gradients, Jacobians, and Hessians for a given NLP.
 * See the interface to Ipopt for an example on how to use the NLPIORACLE.
 *
 * A complete list of all NLPIs contained in this release can be found \ref NLPIS "here".
 *
 * We now explain how users can add their own NLP solver interface.
 * Take the interface to Ipopt (src/nlpi/nlpi_ipopt.cpp) as an example.
 * Unlike most other plugins, it is written in C++.
 * Additional documentation for the callback methods of an NLPI, in particular for their input parameters,
 * can be found in the file type_nlpi.h.
 *
 * Here is what you have to do to implement an NLPI:
 * -# Copy the template files src/nlpi/nlpi_xyz.c and src/nlpi/nlpi_xyz.h into files named "nlpi_mynlpi.c"
 *    and "nlpi_mynlpi.h".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor and replace all occurrences of "xyz" by "mynlpi".
 * -# Adjust the properties of the nlpi (see \ref NLPI_PROPERTIES).
 * -# Define the NLPI and NLPIPROBLEM data (see \ref NLPI_DATA).
 * -# Implement the interface methods (see \ref NLPI_INTERFACE).
 * -# Implement the fundamental callback methods (see \ref NLPI_FUNDAMENTALCALLBACKS).
 *
 *
 * @section NLPI_PROPERTIES Properties of an NLPI
 *
 * At the top of the new file "nlpi_mynlpi.c", you can find the NLPI properties.
 * These are given as compiler defines.
 * The properties you have to set have the following meaning:
 *
 * \par NLPI_NAME: the name of the NLP solver interface.
 * This name is used in the interactive shell to address the NLPI.
 * Additionally, if you are searching for an NLPI with SCIPfindNLPI(), this name is looked up.
 * Names have to be unique: no two NLPIs may have the same name.
 *
 * \par NLPI_DESC: the description of the NLPI.
 * This string is printed as a description of the NLPI in the interactive shell.
 *
 * \par NLPI_PRIORITY: the priority of the NLPI.
 * If an NLP has to be solved, an NLP solver has to be selected.
 * By default, the solver with the NLPI with highest priority is selected.
 * The priority of an NLPI should be set according to performance of the solver:
 * solvers that provide fast algorithms that are usually successful on a wide range of problems should have a high priority.
 * An easy way to list the priorities of all NLPIs is to type "display nlpis" in the interactive shell of SCIP.
 *
 * @section NLPI_DATA NLPI Data
 *
 * Below the header "Data structures" you can find structs which are called "struct SCIP_NlpiData" and "struct SCIP_NlpiProblem".
 * In this data structure, you can store the data of your solver interface and of a specific NLP problem.
 * For example, you could store a pointer to the block memory data structure in the SCIP_NlpiData data structure
 * and store a pointer to an NLPIoracle in the SCIP_NlpiProblem data structure.
 *
 * @section NLPI_INTERFACE Interface Methods
 *
 * At the bottom of "nlpi_mynlpi.c", you can find the interface method SCIPcreateNlpSolverXyz(),
 * which also appears in "nlpi_mynlpi.h".
 * \n
 * This method only has to be adjusted slightly.
 * It is responsible for creating an NLPI that contains all properties and callback methods of your
 * solver interface by calling the method SCIPnlpiCreate().
 * SCIPcreateNlpSolverXyz() is called by the user (e.g., SCIP), if (s)he wants to use this solver interface in his/her application.
 *
 * If you are using NLPI data, you have to allocate the memory for the data at this point.
 * You can do this by calling:
 * \code
 * SCIP_CALL( SCIPallocMemory(scip, &nlpidata) );
 * \endcode
 * You also have to initialize the fields in struct SCIP_NlpiData afterwards. For freeing the
 * NLPI data, see \ref NLPIFREE.
 *
 *
 * @section NLPI_FUNDAMENTALCALLBACKS Fundamental Callback Methods of an NLPI
 *
 * The fundamental callback methods of the plugins are the ones that have to be implemented in order to obtain
 * an operational algorithm. Currently, all NLPI callbacks are fundamental.
 *
 * Additional documentation of the callback methods, in particular to their input parameters,
 * can be found in type_nlpi.h.
 *
 * @subsection NLPICOPY
 *
 * The NLPICOPY callback is executed if the plugin should be copied, e.g., when a SCIP instance is copied.
 *
 * @subsection NLPIFREE
 *
 * The NLPIFREE callback is executed if the NLP solver interface data structure should be freed, e.g., when a SCIP instance is freed.
 *
 * @subsection NLPIGETSOLVERPOINTER
 *
 * The NLPIGETSOLVERPOINTER callback can be used to pass a pointer to a solver specific data structure to the user.
 *
 * @subsection NLPICREATEPROBLEM
 *
 * The NLPICREATEPROBLEM callback is executed if a particular NLP problem is to be created.
 * The callback method should initialize a SCIP_NlpiProblem struct here that corresponds to an empty NLP.
 *
 * @subsection NLPIFREEPROBLEM
 *
 * The NLPIFREEPROBLEMPOINTER callback is executed if a particular NLP problem is to be freed.
 * The callback method should free a SCIP_NlpiProblem struct here.
 *
 * @subsection NLPIGETPROBLEMPOINTER
 *
 * The NLPIGETPROBLEMPOINTER callback can be used to pass a pointer to a solver specific data structure of the NLP to the user.
 *
 * @subsection NLPIADDVARS
 *
 * The NLPIADDVARS callback is executed if a set of variables with lower and upper bounds and names should be added to a particular NLP.
 * The callback method must add the new variables behind the previously added variables, if any.
 * If NULL is given for the lower bounds arguments, -infinity is assumed as lower bound for each new variable.
 * If NULL is given for the upper bounds arguments, +infinity is assumed as upper bound for each new variable.
 * It is also permitted to use NULL for the names argument.
 *
 * @subsection NLPIADDCONSTRAINTS
 *
 * The NLPIADDCONSTRAINTS callback is executed if a set of constraints should be added to a particular NLP.
 * Constraints are specified by providing left and right hand sides, linear and quadratic coefficients, expression trees, and constraint names.
 * All of these arguments are optional, giving NULL for left hand sides corresponds to -infinity, giving NULL for right hand sides corresponds to +infinity.
 *
 * @subsection NLPISETOBJECTIVE
 *
 * The NLPISETOBJECTIVE callback is executed to set the objective function of a particular NLP.
 *
 * @subsection NLPICHGVARBOUNDS
 *
 * The NLPICHGVARBOUNDS callback is executed to change the bounds on a set of variables of an NLP.
 *
 * @subsection NLPICHGCONSSIDES
 *
 * The NLPICHGCONSSIDES callback is executed to change the sides on a set of constraints of an NLP.
 *
 * @subsection NLPIDELVARSET
 *
 * The NLPIDELVARSET callback is executed to delete a set of variables from an NLP.
 * The caller provides an array in which for each variable it is marked whether it should be deleted.
 * In the same array, the method should return the new position of each variable in the NLP, or -1 if it was deleted.
 *
 * @subsection NLPIDELCONSSET
 *
 * The NLPIDELCONSSET callback is executed to delete a set of constraints from an NLP.
 * The caller provides an array in which for each constraint it is marked whether it should be deleted.
 * In the same array, the method should return the new position of each constraint in the NLP, or -1 if it was deleted.
 *
 * @subsection NLPICHGLINEARCOEFS
 *
 * The NLPICHGLINEARCOEFS callback is executed to change the coefficients in the linear part of the objective function or a constraint of an NLP.
 *
 * @subsection NLPICHGQUADCOEFS
 *
 * The NLPICHGQUADCOEFS callback is executed to change the coefficients in the quadratic part of the objective function or a constraint of an NLP.
 *
 * @subsection NLPICHGEXPRTREE
 *
 * The NLPICHGEXPRTREE callback is executed to replace the expression tree of the objective function or a constraint of an NLP.
 *
 * @subsection NLPICHGNONLINCOEF
 *
 * The NLPICHGNONLINCOEF callback is executed to change a single parameter in the (parametrized) expression tree of the objective function or a constraint of an NLP.
 *
 * @subsection NLPICHGOBJCONSTANT
 *
 * The NLPICHGOBJCONSTANT callback is executed to change the constant offset of the objective function of an NLP.
 *
 * @subsection NLPISETINITIALGUESS
 *
 * The NLPISETINITIALGUESS callback is executed to provide primal and dual initial values for the variables and constraints of an NLP.
 * For a local solver, these values can be used as a starting point for the search.
 * It is possible to pass a NULL pointer for any of the arguments (primal values of variables, dual values of variable bounds, dual values of constraints).
 * In this case, the solver should clear previously set starting values and setup its own starting point.
 *
 * @subsection NLPISOLVE
 *
 * The NLPISOLVE callback is executed when an NLP should be solved.
 * The solver may use the initial guess provided by \ref NLPISETINITIALGUESS as starting point.
 * The status of the solving process and solution can be requested by
 * \ref NLPIGETSOLSTAT, \ref NLPIGETTERMSTAT, \ref NLPIGETSOLUTION, and \ref NLPIGETSTATISTICS.
 *
 * @subsection NLPIGETSOLSTAT
 *
 * The NLPIGETSOLSTAT callback can be used to request the solution status (solved, infeasible, ...) after an NLP has been solved.
 *
 * @subsection NLPIGETTERMSTAT
 *
 * The NLPIGETTERMSTAT callback can be used to request the termination reason (normal, iteration limit, ...) after an NLP has been solved.
 *
 * @subsection NLPIGETSOLUTION
 *
 * The NLPIGETSOLUTION callback can be used to request the primal and dual solution values after an NLP solve.
 * The method should pass pointers to arrays of variable values to the caller.
 * It is possible to return only primal values for the variables, but no values for the dual variables, e.g., if a solver does not compute such values.
 *
 * @subsection NLPIGETSTATISTICS
 *
 * The NLPIGETSTATISTICS callback can be used to request the statistical values (number of iterations, time, ...) after an NLP solve.
 * The method should fill the provided NLPSTATISTICS data structure.
 *
 * <!-- NLPIGETWARMSTARTSIZE, NLPIGETWARMSTARTMEMO, NLPISETWARMSTARTMEMO are not documented,
      since they are currently not used, not implemented, and likely to change with a next version. -->
 *
 * @subsection NLPIGETINTPAR
 *
 * The NLPIGETINTPAR callback can be used to request the value of an integer valued NLP parameter.
 *
 * @subsection NLPISETINTPAR
 *
 * The NLPISETINTPAR callback is executed to set the value of an integer valued NLP parameter.
 *
 * @subsection NLPIGETREALPAR
 *
 * The NLPIGETREALPAR callback can be used to request the value of a real valued NLP parameter.
 *
 * @subsection NLPISETREALPAR
 *
 * The NLPISETREALPAR callback is executed to set the value of a real valued NLP parameter.
 *
 * @subsection NLPIGETSTRINGPAR
 *
 * The NLPIGETSTRINGPAR callback can be used to request the value of a string valued NLP parameter.
 *
 * @subsection NLPISETSTRINGPAR
 *
 * The NLPISETSTRINGPAR callback is executed to set the value of a string valued NLP parameter.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page EXPRINT How to add interfaces to expression interpreters
 *
 * An expression interpreter is a tool to compute point-wise and interval-wise the function values, gradients, and
 * derivatives of algebraic expressions which are given in the form of an expression tree.
 * It is used, e.g., by an NLP solver interface to compute Jacobians and Hessians for the solver.
 *
 * The expression interpreter interface in SCIP has been implemented similar to those of the LP solver interface (LPI).
 * For one binary, exactly one expression interpreter has to be linked.
 * The expression interpreter API has been designed such that it can be used independently from SCIP.
 *
 * A complete list of all expression interpreters contained in this release can be found \ref EXPRINTS "here".
 *
 * We now explain how users can add their own expression interpreters.
 * Take the interface to CppAD (\ref exprinterpret_cppad.cpp) as an example.
 * Unlike most other plugins, it is written in C++.
 *
 * Additional documentation for the callback methods of an expression interpreter, in particular for their input parameters,
 * can be found in the file \ref exprinterpret.h
 *
 * Note that the expression interpreter API has <b>BETA status</b> and thus may change in the next version.
 *
 * Here is what you have to do to implement an expression interpreter:
 * -# Copy the file \ref exprinterpret_none.c into a file named "exprinterpreti_myexprinterpret.c".
 *    \n
 *    Make sure to adjust your Makefile such that these files are compiled and linked to your project.
 * -# Open the new files with a text editor.
 * -# Define the expression interpreter data (see \ref EXPRINT_DATA).
 * -# Implement the interface methods (see \ref EXPRINT_INTERFACE).
 *
 *
 * @section EXPRINT_DATA Expression Interpreter Data
 *
 * In "struct SCIP_ExprInt", you can store the general data of your expression interpreter.
 * For example, you could store a pointer to the block memory data structure.
 *
 * @section EXPRINT_INTERFACE Interface Methods
 *
 * The expression interpreter has to implement a set of interface method.
 * In your "exprinterpret_myexprinterpret.c", these methods are mostly dummy methods that return error codes.
 *
 * @subsection SCIPexprintGetName
 *
 * The SCIPexprintGetName method should return the name of the expression interpreter.
 *
 * @subsection SCIPexprintGetDesc
 *
 * The SCIPexprintGetDesc method should return a short description of the expression interpreter, e.g., the name of the developer of the code.
 *
 * @subsection SCIPexprintGetCapability
 *
 * The SCIPexprintGetCapability method should return a bitmask that indicates the capabilities of the expression interpreter,
 * i.e., whether it can evaluate gradients, Hessians, or do interval arithmetic.
 *
 * @subsection SCIPexprintCreate
 *
 * The SCIPexprintCreate method is called to create an expression interpreter data structure.
 * The method should initialize a "struct SCIP_ExprInt" here.
 *
 * @subsection SCIPexprintFree
 *
 * The SCIPexprintFree method is called to free an expression interpreter data structure.
 * The method should free a "struct SCIP_ExprInt" here.
 *
 * @subsection SCIPexprintCompile
 *
 * The SCIPexprintCompile method is called to initialize the data structures that are required to evaluate
 * a particular expression tree.
 * The expression interpreter can store data that is particular to a given expression tree in the tree by using
 * SCIPexprtreeSetInterpreterData().
 *
 * @subsection SCIPexprintFreeData
 *
 * The SCIPexprintFreeData method is called when an expression tree is freed.
 * The expression interpreter should free the given data structure.
 *
 * @subsection SCIPexprintNewParametrization
 *
 * The SCIPexprintNewParametrization method is called when the values of the parameters in a parametrized expression tree have changed.
 *
 * @subsection SCIPexprintEval
 *
 * The SCIPexprintEval method is called when the value of an expression represented by an expression tree should be computed for a point.
 *
 * @subsection SCIPexprintEvalInt
 *
 * The SCIPexprintEvalInt method is called when an interval that contains the range of an expression represented by an expression tree with respect to intervals for the variables should be computed.
 *
 * @subsection SCIPexprintGrad
 *
 * The SCIPexprintGrad method is called when the gradient of an expression represented by an expression tree should be computed for a point.
 *
 * @subsection SCIPexprintGradInt
 *
 * The SCIPexprintGradInt method is called when an interval vector that contains the range of the gradients of an expression represented by an expression tree with respect to intervals for the variables should be computed.
 *
 * @subsection SCIPexprintHessianSparsityDense
 *
 * The SCIPexprintHessianSparsityDense method is called when the sparsity structure of the Hessian matrix should be computed and returned in dense form.
 *
 * @subsection SCIPexprintHessianDense
 *
 * The SCIPexprintHessianDense method is called when the Hessian of an expression represented by an expression tree should be computed for a point.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CONF How to use conflict analysis
 *
 * Conflict analysis is a way to automatically use the information obtained from infeasible nodes
 * in the branch-and-bound tree.
 *
 * Once a node is declared infeasible, SCIP automatically tries to infer a constraint that explains the reason for the
 * infeasibility, in order to avoid similar situations later in the search.  This explanation essentially consists of a
 * constraint stating that at least one of its variables should have a bound different from the current infeasible node,
 * because the current setting led to infeasibility. Clearly, all variables that are fixed in the current infeasible
 * node would yield such a constraint (since this leads to infeasibility). The key point rather is to infer a "small"
 * constraint that does the same job. SCIP handles this by several heuristics. For this, SCIP sets up a
 * so-called (directed) conflict graph. The nodes in this graph correspond to bound changes of variables and an arc (@a
 * u, @a v) means that the bound change corresponding to @a v is based on the bound change of @a u. In general, a node
 * will have several ingoing arcs which represent all bound changes that have been used to infer (propagate) the bound
 * change in question. The graph also contains source nodes for each bound that has been changed during branching and an
 * artificial target node representing the conflict, i.e., the infeasibility. Essentially, SCIP heuristically constructs
 * a cut in this graph that involves few "branching nodes". For details on the techniques that SCIP uses,
 * we refer to the paper @par
 * Tobias Achterberg, Conflict Analysis in Mixed Integer Programming@n
 * Discrete Optimization, 4, 4-20 (2007)
 *
 * For conflict analysis to work well, the author of a \ref CONS "Constraint Handler" or a
 * \ref PROP "Propagator" has to implement three kinds of functionality:
 *
 * -# If one detects infeasibility, one should initiate conflict analysis, see \ref INITCONFS "below".
 * -# During propagation, one should call the right functions to fix variables.
 * -# One should implement the <em>so-called reverse propagation</em>.
 *
 * If this functionality is not implemented, SCIP will still work correctly, but cannot use the information of the constraint
 * handler or the propagator for conflict analysis. In this case, each bound reduction performed by the constraint
 * handler/propagator will be treated as if it had been a branching decision.
 *
 * @section INITCONFS Initiating Conflict Analysis
 *
 * If one detects infeasibility within propagation, one should do the following:
 * -# Call SCIPinitConflictAnalysis().
 * -# Inform SCIP about the variable bounds that are the reason for the detection of infeasibility
 * via the functions SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(), or
 * SCIPaddConflictBinvar(). If there is more than one valid explanation of infeasibility, either one can be used.
 * Typically, smaller explanations tend to be better.
 * -# Call SCIPanalyzeConflict() from a propagator or SCIPanalyzeConflictCons() from a constraint
 * handler.
 *
 * This functionality allows SCIP to set up the conflict graph and perform a conflict analysis.
 *
 * @section Propagation
 *
 * When propagating variable domains, SCIP needs to be informed that the deduced variable bounds should be
 * used in conflict analysis. This can be done by the functions SCIPinferVarLbCons(),
 * SCIPinferVarUbCons(), and SCIPinferBinvarCons() for constraint handlers and SCIPinferVarLbProp(),
 * SCIPinferVarUbProp(), and SCIPinferBinvarProp() for propagators. You can pass one integer of
 * information that should indicate the reason of the propagation and can be used in reverse
 * propagation, see the next section.
 *
 * @section RESPROP Reverse Propagation
 *
 * Reverse Propagation is used to build up the conflict graph. Essentially, it provides an algorithm to detect the arcs
 * leading to a node in the conflict graph, i.e., the bound changes responsible for the new bound change deduced during
 * propagation. Reverse Propagation needs to be implemented in the RESPROP callback functions of
 * \ref CONSRESPROP "constraint handlers" or \ref PROPRESPROP "propagators".
 * These callbacks receive the following information: the variable which is under investigation (@p
 * infervar), the corresponding bound change (@p bdchgidx, @p boundtype), and the integer (@p inferinfo) that has been
 * supplied during propagation.
 *
 * One can use SCIPvarGetUbAtIndex() or SCIPvarGetLbAtIndex() to detect the bounds before or after the propagation that
 * should be investigated. Then the bounds that were involved should be passed to SCIP via SCIPaddConflictLb() and
 * SCIPaddConflictUb().  If there is more than one valid explanation of infeasibility, either one can be used.
 * Typically, smaller explanations tend to be better.
 *
 * Details and (more) examples are given in Sections @ref CONSRESPROP and @ref PROPRESPROP.
 *
 *
 * @section Example
 *
 * Consider the constraint handler @p cons_linearordering.c in the
 * <a href="http://scip.zib.de/doc/examples/LOP/index.shtml"><b>linear ordering example</b></a>
 * (see @p example/LOP directory). This constraint handler propagates the equations \f$x_{ij} + x_{ji} =
 * 1\f$ and triangle inequalities \f$x_{ij} + x_{jk} + x_{ki} \leq 2\f$.
 *
 * When propagating the equation and <code>vars[i][j]</code> is fixed to 1, the constraint handler uses
 * \code
 *    SCIP_CALL( SCIPinferBinvarCons(scip, vars[j][i], FALSE, cons, i*n + j, &infeasible, &tightened) );
 * \endcode
 * Thus, variable <code>vars[j][i]</code> is fixed to 0 (@p FALSE), and it passes <code>i*n + j </code> as @p inferinfo.
 *
 * When it propagates the triangle inequality and both <code>vars[i][j]</code> and <code>vars[j][k]</code>
 * are fixed to 1, the constraint handler uses
 * \code
 *    SCIP_CALL( SCIPinferBinvarCons(scip, vars[k][i], FALSE, cons, n*n + i*n*n + j*n + k, &infeasible, &tightened) );
 * \endcode
 * Thus, in this case, variable  <code>vars[k][i]</code> is fixed to 0 and  <code>n*n + i*n*n +  j*n + k</code> is
 * passed as <code>inferinfo</code>.
 *
 * In reverse propagation, the two cases can be distinguished by @p inferinfo: if it is less than @p n*n,
 * we deal with an equation, otherwise with a triangle inequality. The constraint handler can then extract the
 * indices @p i, @p j (and @p k in the second case) from inferinfo.
 *
 * In the first case, it has to distinguish whether <code>vars[i][j]</code> is fixed to 0 or 1 &ndash;
 * by calling SCIPaddConflictLb()
 * or SCIPaddConflictUb(), respectively, with variable <code>vars[j][i]</code>. In the second case, it is clear that the only
 * possible propagation is to fix <code>vars[i][j]</code> to 0 when both <code>vars[k][i]</code> and <code>vars[j][k]</code>
 * are fixed to 1. It then calls
 * SCIPaddConflictLb() for both <code>vars[k][i]</code> and <code>vars[j][k]</code>.
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
 *  doesn't need the object anymore, (s)he has to call the object's release method.
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
/**@page PARAM How to add additional user parameters
 *
 *  Users may add their own parameters to SCIP by calling SCIPaddXyzParam(). Using
 *  this method, there are two possibilities for where to store the actual parameter value:
 *   - If the given valueptr is NULL, SCIP stores the parameter value internally, and
 *     the user can only access the value with the SCIPgetXyzParam() and
 *     SCIPsetXyzParam() calls.
 *   - If the given valueptr is not NULL, SCIP stores the parameter value at the given
 *     address, and the user can directly manipulate the value at this address.
 *     (S)he has to be careful with memory management in string parameters: when the
 *     SCIPaddStringParam() method is called, the given address must hold a char*
 *     pointer with value NULL. The default value is then copied into this pointer,
 *     allocating memory with BMSallocMemoryArray(). If the parameter is changed, the
 *     old string is freed with BMSfreeMemoryArray() and the new one is copied to a new
 *     memory area allocated with BMSallocMemoryArray(). When the parameter is freed,
 *     the memory is freed with BMSfreeMemoryArray().
 *     The user should not interfere with this internal memory management. Accessing
 *     the string parameter through the given valueptr is okay as long as it does not
 *     involve reallocating memory for the string.
 *
 *  In some cases, it is necessary to keep track of changes in a parameter.
 *  If this is the case, the user can define a method by the PARAMCHGD callback and use this method as
 *  the @c paramchgd parameter of the @c SCIPaddXyzParam() method, also giving a pointer to the data, which is
 *  needed in this method, as @c paramdata. If this method is not NULL, it is called every time
 *  the value of the parameter is changed.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page DEBUG Debugging
 *
 *  If you need to debug your own code that uses SCIP, here are some tips and tricks:

 *  - Use <b>asserts</b> in your code to show preconditions for the parameters, invariants and postconditions.
 *    Assertions are boolean expressions which inevitably have to evaluate to <code>TRUE</code>. Consider the
 *    following example, taken from the file src/scip/cons_linear.c:
 * \verbatim
SCIP_RETCODE consdataCatchEvent(
   SCIP*                 scip,               /**< SCIP data structure *\/
   SCIP_CONSDATA*        consdata,           /**< linear constraint data *\/
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing *\/
   int                   pos                 /**< array position of variable to catch bound change events for *\/
   )
   {
      assert(scip != NULL);
      assert(consdata != NULL);
      assert(eventhdlr != NULL);
      assert(0 <= pos && pos < consdata->nvars);
   ...
   }
 * \endverbatim
 *    @n
 *    As you can see, both pointers and integers are checked for valid values at the beginning of the
 *    function <code>consdataCatchEvent()</code>. This is particularly important for, e.g., array indices like
 *    the variable <code>pos</code> in this example, where using the <code>consdata->nvars[pos]</code>
 *    pointer could result in unexspected behaviour
 *    if the asserted precondition on <code>pos</code> were not matched and \<pos\> were an arbitrary index
 *    outside the array range.
 *
 *  - In order to activate assertions, use the <b>Debug mode</b> by compiling SCIP via
 *   \code
 *    make OPT=dbg
 *   \endcode and run the code. See \ref MAKE for further information about compiler options for SCIP.
 *
 *  - Spending only little extra time on
 *    asserting preconditions saves most of the time spent on debugging!
 *
 *  - Turn on <b>additional debug output</b> by adding the line
 *    \code
 *    #define SCIP_DEBUG
 *    \endcode
 *    at the top of SCIP files you want to analyze. This will output messages included in the code using
 *    <code>SCIPdebugMessage()</code> (see \ref EXAMPLE_1).
 *    We recommend to also use <code>SCIPdebugMessage()</code> in your own code for being able to activate
 *    debug output in the same way.
 *  - If available on your system, we recommend to use a debugger like <code>gdb</code>
 *    to trace all function calls on the stack,
 *    display values of certain expressions, manually break the running code, and so forth.
 *  - If available on your system, you can use software like <a href="http://valgrind.org">valgrind</a> to check for uninitialized
 *    values or segmentation faults.
 *  - For checking the usage of SCIP memory, you can use
 *    <code>SCIPprintMemoryDiagnostic()</code>. This outputs memory that is currently in use,
 *    which can be useful after a <code>SCIPfree()</code> call.
 *  - If your code cuts off a feasible solution, but you do not know which component is responsible,
 *    you can define <code>SCIP_DEBUG_SOLUTION</code> in the file <code>debug.h</code> to be a filename
 *    containing a solution in SCIP format (see \ref EXAMPLE_2).
 *    This solution is then read and it is checked for every cut, whether the solution violates the cut.
 *
 *  @section EXAMPLE_1 How to activate debug messages
 *    For example, if we include a <code>\#define SCIP_DEBUG</code> at the top of \ref heur_oneopt.h, recompile SCIP
 *    in DBG mode, and run the SCIP interactive shell to solve p0033.mps from the
 *     <a href="http://miplib.zib.de/miplib3/miplib.html">MIPLIB 3.0</a> , we get some output like:
 * \code
 * SCIP version 1.1.0 [precision: 8 byte] [memory: block] [mode: debug] [LP solver: SoPlex 1.4.0]
 * Copyright (c) 2002-2013 Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
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
 * @section EXAMPLE_2 How to add a debug solution
 *
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
 * If we afterwards use
 * <code>\#define SCIP_DEBUG_SOLUTION "check/p0033.sol"</code> in debug.h, recompile and run SCIP,
 * it will output:
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
 *  SCIP comes along with a set of useful tools that allow to perform automated tests. The
 *  following is a step-by-step guide from setting up the test environment for evaluation and
 *  customization of test runs.
 *
 *
 *  @section SETUP Setting up the test environment
 *
 *  At first you should create a file listing all problem instances that should be part of the test.
 *  This file has to be located in the the directory <code>scip/check/testset/</code>
 *  and has to have the file extension <code>.test</code>, e.g., <code>testrun.test</code>,
 *  in order to be found by the <code>scip/check/check.sh</code> script.
 *  \n
 *  All test problems can be listed in the <code>test</code>-file by a relative path,
 *  e.g., <code>../../problems/instance1.lp</code> or absolute path, e.g., <code>/home/problems/instance2.mps</code>
 *  in this file. Only one problem should be listed each on line (since the command <code>cat</code> is used to parse this file).
 *  Note that these problems have to be readable for SCIP in order to solve them.
 *  However, you can use different file formats.
 *
 *  Optionally, you can provide a solution file in the <code>scip/check/testset/</code> directory containing
 *  known information about the feasibility and the best known objective values for the test instances.
 *  SCIP can use these values to verify the results. The file has to have the same basename as the
 *  <code>.test</code>-file, i.e., in our case <code>testrun.solu</code>. One line can only contain
 *  information about one test instance. A line has to start with the type of information given:
 *
 *  - <code>=opt=</code> stating that a problem name with an optimal objective value follows
 *  - <code>=best=</code> stating that a problem name with a best know objective value follows
 *  - <code>=inf=</code> stating that a problem name follows which is infeasible
 *
 *  With these information types you can encode for an instance named <code>instance1.lp</code> the following
 *  information:
 *  - The instance has a known optimal (objective) value of 10.
 *   \code
 *   =opt=  instance1 10
 *   \endcode
 *  - The instance has a best known solution with objective value 15.
 *   \code
 *   =best=  instance1 15
 *   \endcode
 *  - The instance is feasible (but has no objective function or we don't know a solution value)
 *   \code
 *   =feas=  instance1
 *   \endcode
 *  - The instance is infeasible.
 *   \code
 *   =inf=  instance1
 *   \endcode
 *
 *  If you don't know whether the instance is feasible or not (so the status is unknown),
 *  you can omit the instance in the <code>solu</code>-file or write
 *   \code
 *   =unkn=  instance1
 *   \endcode
 *
 * <b>Note that in all lines the file extension of the file name is omitted.</b>
 *  \n
 *  See the files <code>scip/check/testset/short.test</code> and <code>scip/check/testset/short.solu</code>
 *  for an example of a <code>test</code>-file and its corresponding <code>solu</code>-file.
 *
 *
 *
 *  @section STARTING Starting a test run
 *
 *
 *  \code
 *  make TEST=testrun test
 *  \endcode
 *
 *  in the SCIP root directory. Note that <code>testrun</code> is exactly the basename of our
 *  <code>test</code>-file (<code>testrun.test</code>). This will cause SCIP to solve our test instances
 *  one after another and to create various output files (see \ref EVAL).
 *
 *
 *  @section EVAL Evaluating a test run
 *
 *  During computation, SCIP automatically creates the directory <code>scip/check/results/</code>
 *  (if it does not already exist) and stores the following output files there.
 *
 *  \arg <code>*.out</code> - output of <code>stdout</code>
 *  \arg <code>*.err</code> - output of <code>stderr</code>
 *  \arg <code>*.set</code> - copy of the used settings file
 *
 *  \arg <code>*.res</code> - ASCII table containing a summary of the computational results
 *  \arg <code>*.tex</code> - TeX table containing a summary of the computational results
 *  \arg <code>*.pav</code> - <a href="http://www.gamsworld.org/performance/paver/">PAVER</a> output
 *
 *  The last three files in the above list, i.e., the files containing a summary of the computational results,
 *  can also be generated manually. Therefore the user has to call the <code>evalcheck.sh</code> script in the
 *  @c check directory with the corresponding @c out file as argument. For example, this may be useful if the user stopped the
 *  test before it was finished, in which case the last three files will not be automatically generated by SCIP.
 *
 *  The last column of the ASCII summary table contains the solver status. We distinguish the following statuses: (in order of priority)
 *
 *  \arg abort: solver broke before returning solution
 *  \arg fail: solver cut off a known feasible solution (value of the <code>solu</code>-file is beyond the dual bound;
 *  especially of problem is claimed to be solved but solution is not the optimal solution)
 *  \arg ok: solver solved problem with the value in solu-file
 *  \arg solved: solver solved problem which has no (optimal) value in solu-file (since we here cannot detect the direction
 *  of optimization, it is possible that a solver claims an optimal solution which contradicts a known feasible solution)
 *  \arg better: solver found solution better than known best solution (or no solution was noted in the <code>solu</code>-file so far)
 *  \arg gaplimit, sollimit: solver reached gaplimit or limit of number of solutions (at present: only in SCIP)
 *  \arg timeout: solver reached any other limit (like time or nodes)
 *  \arg unknown: otherwise
 *
 *  Additionally the <code>evalcheck.sh</code> script can generate a <code>solu</code>-file by calling
 *  \code
 *  ./evalcheck.sh writesolufile=1 NEWSOLUFILE=<solu-file> <out-file>
 *  \endcode
 *  where <code><solu-file></code> denotes the filename of the new file where the solutions shall be
 *  (and <code><out-file></code> denotes the output (<code>.out</code>) files to evaluate).
 *
 *  Another feature can be enabled by calling:
 *  \code
 *  ./evalcheck.sh printsoltimes=1 ...
 *  \endcode
 *  The output has two additional columns containing the solving time until the first and the best solution was found.
 *
 *
 *  @b Note: The @em basename of all these files is the same and has the following structure
 *  which allows us to reconstruct the test run:
 *
 *  \code
 *  check.<test name>.<binary>.<machine name>.<setting name>
 *  \endcode
 *
 *  \arg <<code>test name</code>> indicates the name of the the test file, e.g., <code>testrun</code>
 *  \arg <<code>binary</code>> defines the used binary, e.g., <code>scip-1.1.0.linux.x86.gnu.opt.spx</code>
 *  \arg <<code>machine name</code>> tells the name of the machine, e.g., <code>mycomputer</code>
 *  \arg <<code>setting name</code>> denotes the name of the used settings, e.g., <code>default</code>
 *    means the (SCIP) default settings were used
 *
 *  Using the examples out of the previous listing the six file names would have the name:
 *
 *  \code
 *  check.testrun.scip-1.1.0.linux.x86.gnu.opt.spx.mycomputer.default.<out,err,set,res,tex,pav>
 *  \endcode
 *
 *
 *  @section USING Using customized setting files
 *
 *  It is possible to use customized settings files for the test run instead of testing SCIP with default settings.
 *  These have to be placed in the directory <code>scip/settings/</code>.
 *
 *  @b Note: Accessing setting files in subfolders of the @c settings directory is currently not supported.
 *
 *  To run SCIP with a custom settings file, say for example <code>fast.set</code>, we call
 *
 *  \code
 *  make TEST=testrun SETTINGS=fast test
 *  \endcode
 *
 *  in the SCIP root directory.
 *
 *
 *  @section ADVANCED Advanced options
 *
 *  We can further customize the test run by specifying the following options in the <code>make</code> call:
 *
 *  \arg <code>TIME</code>  - time limit for each test instance in seconds [default: 3600]
 *  \arg <code>NODES</code> - node limit [default: 2100000000]
 *  \arg <code>MEM</code>   -  memory limit in MB [default: 1536]
 *  \arg <code>DISPFREQ</code> - display frequency of the output [default: 10000]
 *  \arg <code>FEASTOL</code> - LP feasibility tolerance for constraints [default: "default"]
 *  \arg <code>LOCK</code> - should the test run be locked to prevent other machines from performing the same test run [default: "false"]
 *  \arg <code>CONTINUE</code> - continue the test run if it was previously aborted [default: "false"]
 *  \arg <code>VALGRIND</code> - run valgrind on the SCIP binary; errors and memory leaks found by valgrind are reported as fails [default: "false"]
 *
 *
 *  @section COMPARE Comparing test runs for different settings
 *
 *  Often test runs are performed on the basis of different settings. In this case, it is useful to
 *  have a performance comparison. For this purpose, we can use the <code>allcmpres.sh</code> script in
 *  the @c check directory.
 *
 *  Suppose, we performed our test run with two different settings, say <code>fast.set</code> and
 *  <code>slow.set</code>. Assuming that all other parameters (including the SCIP binary), were the same,
 *  we may have the following <code>res</code>-files in the directory <code>scip/check/results/</code>
 *
 *  \code
 *  check.testrun.scip-1.1.0.linux.x86.gnu.opt.spx.mycomputer.fast.res
 *  check.testrun.scip-1.1.0.linux.x86.gnu.opt.spx.mycomputer.slow.res
 *  \endcode
 *
 *  For a comparison of both computations, we simply call
 *
 *  \code
 *  allcmpres.sh results/check.testrun.scip-1.1.0.linux.x86.gnu.opt.spx.mycomputer.fast.res \
 *               results/check.testrun.scip-1.1.0.linux.x86.gnu.opt.spx.mycomputer.slow.res
 *  \endcode
 *
 *  in the @c check directory. This produces an ASCII table on the console that provide a detailed
 *  performance comparison of both test runs. Note that the first <code>res</code>-file serves as reference
 *  computation. The following list explains the output.
 *  (The term "solver" can be considered as the combination of SCIP with a specific setting file.)
 *
 *  \arg <code>Nodes</code> - Number of processed branch-and-bound nodes.
 *  \arg <code>Time</code>  - Computation time in seconds.
 *  \arg <code>F</code>     - If no feasible solution was found, then '#', empty otherwise.
 *  \arg <code>NodQ</code>  - Equals Nodes(i) / Nodes(0), where 'i' denotes the current solver and '0' stands for the reference solver.
 *  \arg <code>TimQ</code>  - Equals Time(i) / Time(0).
 *  \arg <code>bounds check</code> - Status of the primal and dual bound check.
 *
 *  \arg <code>proc</code> - Number of instances processed.
 *  \arg <code>eval</code> - Number of instances evaluated (bounds check = "ok", i.e., solved to optimality
 *      within the time and memory limit and result is correct). Only these instances are used in the calculation
 *      of the mean values.
 *  \arg <code>fail</code> - Number of instances with bounds check = "fail".
 *  \arg <code>time</code> - Number of instances with timeout.
 *  \arg <code>solv</code> - Number of instances correctly solved within the time limit.
 *  \arg <code>wins</code> - Number of instances on which the solver won (i.e., the
 *      solver was at most 10% slower than the fastest solver OR had the best
 * 	primal bound in case the instance was not solved by any solver within
 *	the time limit).
 *  \arg <code>bett</code>    - Number of instances on which the solver was better than the
 *	reference solver (i.e., more than 10% faster).
 *  \arg <code>wors</code>    - Number of instances on which the solver was worse than the
 *	reference solver (i.e., more than 10% slower).
 *  \arg <code>bobj</code>    - Number of instances on which the solver had a better primal
 *	bound than the reference solver (i.e., a difference larger than 10%).
 *  \arg <code>wobj</code>    - Number of instances on which the solver had a worse primal
 *	bound than the reference solver (i.e., a difference larger than 10%).
 *  \arg <code>feas</code>    - Number of instances for which a feasible solution was found.
 *  \arg <code>gnodes</code>   - Geometric mean of the processed nodes over all evaluated instances.
 *  \arg <code>shnodes</code> - Shifted geometric mean of the processed nodes over all evaluated instances.
 *  \arg <code>gnodesQ</code>  - Equals nodes(i) / nodes(0), where 'i' denotes the current
 *	solver and '0' stands for the reference solver.
 *  \arg <code>shnodesQ</code> - Equals shnodes(i) / shnodes(0).
 *  \arg <code>gtime</code>    - Geometric mean of the computation time over all evaluated instances.
 *  \arg <code>shtime</code>  - Shifted geometric mean of the computation time over all evaluated instances.
 *  \arg <code>gtimeQ</code>   - Equals time(i) / time(0).
 *  \arg <code>shtimeQ</code> - Equals shtime(i) / shtime(0).
 *  \arg <code>score</code>   - N/A
 *
 *  \arg <code>all</code>   - All solvers.
 *  \arg <code>optimal auto settings</code> - Theoretical result for a solver that performed 'best of all' for every instance.
 *  \arg <code>diff</code>  - Solvers with instances that differ from the reference solver in the number of
 *       processed nodes or in the total number of simplex iterations.
 *  \arg <code>equal</code> - Solvers with instances whose number of processed nodes and total number of
 *       simplex iterations is equal to the reference solver (including a 10% tolerance) and where no timeout
 *       occured.
 *  \arg <code>all optimal</code> - Solvers with instances that could be solved to optimality by
 *       <em>all</em> solvers; in particular, no timeout occurred.
 *
 *  Since this large amount of information is not always needed, one can generate a narrower table by calling:
 *  \code
 *  allcmpres.sh short=1 ...
 *  \endcode
 *  where <code>NodQ</code>, <code>TimQ</code> and the additional comparison tables are omitted.
 *
 *  If the <code>res</code>-files were generated with the parameter <code>printsoltimes=1</code>
 *  we can enable the same feature here as well by calling:
 *  \code
 *  allcmpres.sh printsoltimes=1 ...
 *  \endcode
 *  As in the evaluation, the output contains the two additional columns of the solving time until the first and the best solution was found.
 *
 *  @section SOLVER Testing and Evaluating for other solvers
 *
 *  Analogously to the target <code>test</code> there are further targets to run automated tests with other MIP solvers.
 *  These are:
 *  \arg for <a href="http://www-01.ibm.com/software/integration/optimization/cplex-optimizer/">cplex</a>
 *  \code
 *  make testcplex
 *  \endcode
 *  \arg for <a href="http://www.gurobi.com/">gurobi</a>
 *  \code
 *  make testgurobi
 *  \endcode
 *  \arg for <a href="https://projects.coin-or.org/Cbc">cbc</a>
 *  \code
 *  make testcbc
 *  \endcode
 *  \arg for <a href="http://www.mosek.com/">mosek</a>
 *  \code
 *  make testmosek
 *  \endcode
 *  \arg for <a href="http://www.gnu.org/software/glpk/">glpk</a>
 *  \code
 *  make testglpk
 *  \endcode
 *  \arg for <a href="https://projects.coin-or.org/SYMPHONY">symphony</a>
 *  \code
 *  make testsymphony
 *  \endcode
 *  \arg for <a href="https://projects.coin-or.org/CHiPPS">blis</a>
 *  \code
 *  make testblis
 *  \endcode
 *  \arg for <a href="http://www.gams.com/">gams</a>
 *  \code
 *  make testgams GAMSSOLVER=xyz
 *  \endcode
 *  For this target, the option GAMSSOLVER has to be given to specify the name of a GAMS solver to run, e.g. GAMSSOLVER=SCIP.
 *  Additional advanced options specific to this target are:
 *    GAMS to specify the GAMS executable (default: gams),
 *    GAP to specify a gap limit (default: 0.0),
 *    CLIENTTMPDIR to specify a directory where GAMS should put its scratch files (default: /tmp),
 *    CONVERTSCIP to specify a SCIP which can be used to convert non-gams files into gams format (default: bin/scip, if existing; set to "no" to disable conversion).
 *  The following options are NOT supported (and ignored): MEM, DISPFREQ, FEASTOL, LOCK.
 *
 *  Note: This works only if the referred programs are installed globally on your machine.
 *
 *  The above options like <code>TIME</code> are also available for the other solvers.
 *
 *  For cbc, cplex, gams, and gurobi another advanced option is available:
 *  \arg <code>THREADS</code> - number of threads used in the solution process
 *
 *  After the testrun there should be an <code>.out</code>, an <code>.err</code> and a <code>.res</code> file
 *  with the same basename as described above.
 *
 *  Furthermore you can also use the script <code>allcmpres.sh</code> for comparing results of different solvers.
 *
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CHG1 Interface changes between SCIP 0.9 and SCIP 1.0
 *
 *  @section CHGPARAM New parameters
 *
 * - All functions SCIP<datatype>Param() got a new parameter "isadvanced".
 *   \n
 *   This does not influence the performance of SCIP, but the position of the parameter in the settings menu.
 *   Hence, if you do not care about this, you can assign any value to it.
 *   You should add the corresponding flag to the SCIP<datatype>Param() calls in your own source code.
 *
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CHG2 Interface changes between SCIP 1.0 and SCIP 1.1
 *
 * - SCIPcreateChild() has a new last parameter giving an estimate for value of best feasible solution in the subtree to
 *   be created. One possibility is to use SCIPgetLocalOrigEstimate() for this value.
 *
 * - The callback \ref CONSCHECK in the constraint handlers now has a new parameter <code>printreason</code> that tells
 *   a constraint handler to output the reason for a possible infeasibility of the solution to be checked using
 *   SCIPinfoMessage(). Have a look at one of the constraint handlers implemented in SCIP to see how it works. This
 *   methodology makes it possible to output the reason of a violation in human readable form, for instance, for the check
 *   at the end of a SCIP run, where the obtained best solution is checked against the original formulation.\n This change
 *   often has little effect on C-implementations, since this parameter can be safely ignored with respect to the
 *   correctness of the code. The corresponding C++ method scip::ObjConshdlr::scip_check(), however, has to be extended
 *   and will not compile otherwise.
 *
 * - SCIPcheckSolOrig() is restructured. The last two parameters have changed. They are now bools indicating
 *   whether the reason for the violation should be printed to the standard output and whether all violations should be
 *   printed. This reflects the changes in the constraint handlers above, which allow the automation of the feasibility
 *   test. The pointers to store the constraint handler or constraint are not needed anymore.
 *
 * - New parameters "extension" and "genericnames" in SCIPprintTransProblem(), SCIPprintOrigProblem(),
 *   SCIPwriteOrigProblem(), and SCIPwriteTransProblem() defining the requested format or NULL for default CIP format
 *   and using generic names for the variables and constraints. Examples are
 *   - <code>SCIPprintTransProblem(scip, NULL, NULL, TRUE)</code> displays the transformed problem in CIP format with
 *     generic variables and constraint names
 *   - <code>SCIPprintOrigProblem(scip, NULL, "lp", FALSE)</code> displays the original problem in LP format with
 *     original variables and constraint names.
 *
 * - New callback method SCIP_DECL_READERWRITE(x) in type_reader.h; this method is called to write a problem to file
 *   stream in the format the reader stands for; useful for writing the transformed problem in LP or MPS format. Hence,
 *   also SCIPincludeReader() has changed.
 *
 * - New parameter "conshdlrname" in SCIPincludeLinconsUpgrade().
 *
 * - Added user pointer to callback methods of hash table, see pub_misc.h.
 *
 * - New parameter "extension" in SCIPreadProb(),    defining a desired file format or NULL if file extension should be used.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CHG3 Interface changes between SCIP 1.1 and SCIP 1.2
 *
 *
 * @section CHGCALLBACKS New and changed callbacks
 *
 * - The callback SCIP_DECL_PRICERREDCOST(x) in the \ref PRICER "pricers" has two new parameters:
 *    - A <code>result</code> pointer determines whether the pricer guarantees that there exist no more variables. This allows for early branching.
 *    - A pointer for providing a lower bound.
 *
 * - The \ref CONS "constraint handlers" have two new callback methods (see type_cons.h for more details).
 *    - SCIP_DECL_CONSCOPY(x) - this method can be used to copy a constraint.
 *    - SCIP_DECL_CONSPARSE(x) - this method can be used to parse a constraint in CIP format.
 *
 *  @section CHGINTERFUNC New parameters in interface methods
 *
 * - SCIPcalcMIR() in scip.h has two new parameter "mksetcoefsvalid" and "sol". The parameter "mksetcoefsvalid" stores
 *   whether the coefficients of the mixed knapsack set ("mksetcoefs") computed in SCIPlpCalcMIR() are valid. If the mixed knapsack constraint obtained after aggregating LP rows
 *   is empty or contains too many nonzero elements the generation of the <b>c-MIR cut</b> is aborted in SCIPlpCalcMIR() and "mksetcoefs" is not valid.
 *   The input parameter "sol" can be used to separate a solution different from the LP solution.
 *
 * - SCIPgetVarClosestVlb() and SCIPgetVarClosestVub() in scip.h have a new parameter "sol". It can be used to obtain the <b>closest variable bound</b> w.r.t. a solution different from the LP solution.
 *
 *  @section MISCELLANEOUS Miscellaneous
 *
 * - A significant change for <b>C++ users</b> is that all include files of SCIP
 *   automatically detect C++ mode, i.e., no <code>extern "C"</code> is needed anymore.
 *
 * For further release notes we refer the \ref RELEASENOTES "Release notes".
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CHG4 Interface changes between SCIP 1.2 and SCIP 2.0
 *
 *
 * @section CHGCALLBACKS4 New and changed callbacks
 *
 *
 * - <b>Copying a SCIP instance</b>:
 *      <br>
 *      <br>
 *    - All plugins, like \ref BRANCH "branching rules" and \ref HEUR "primal heuristics", have a new callback method (see, e.g.,
 *      type_branch.h and type_heur.h for more details):
 *       - SCIP_DECL_BRANCHCOPY(x), SCIP_DECL_HEURCOPY(x) etc.
 *       - When copying a SCIP instance, these methods are called to copy the plugins.
 *      <br>
 *      <br>
 *    - Constraint handlers have two new callback methods. One for copying the constraint handler plugins
 *      SCIP_DECL_CONSHDLRCOPY() and the other for copying a constraint itself, SCIP_DECL_CONSCOPY().
 *      <br>
 *      <br>
 *    - Variables have a new callback method (see type_var.h for more details):
 *       - SCIP_DECL_VARCOPY(x) - When copying a SCIP instance, this method is called to copy the variables' data.
 *      <br>
 *      <br>
 *    - The main problem has a new callback method (see type_prob.h for more details):
 *       - SCIP_DECL_PROBCOPY(x) - When copying a SCIP instance, this method is called to copy the problem's data.
 *      <br>
 *      <br>
 *    - The argument success in SCIP_DECL_CONSCOPY has been renamed to valid.
 *
 * - <b>Branching on externally given candidates</b>:
 *      <br>
 *      <br>
 *    - The \ref BRANCH "branching rules" have a second new callback method (see type_branch.h for more details):
 *       - SCIP_DECL_BRANCHEXECEXT(x) - This method can be used to branch on external branching candidates,
 *         which can be added by a user's "relaxation handler" or "constraint handler" plugin, calling <code>SCIPaddExternBranchCand()</code>.
 *
 * - <b>Restarts</b>:
 *      <br>
 *      <br>
 *    - The callback SCIP_DECL_PROBEXITSOL(x) in the main problem has one new parameter (see type_prob.h for more details):
 *       - The parameter <code>restart</code> is <code>TRUE</code> if the callback method was triggered by a restart.
 *
 *
 * <br>
 * @section CHGINTERFUNC4 Changed interface methods
 *
 * - <b>Copying a SCIP instance</b>:
 *      <br>
 *      <br>
 *    - Every new callback method resulted in a new parameter of the include function for the corresponding plugin,
 *      e.g., SCIPincludeBranchrule() has two new parameters <code>SCIP_DECL_BRANCHCOPY((*branchcopy))</code> and
 *      <code>SCIP_DECL_BRANCHEXECREL((*branchexecrel))</code>.  In the same fashion, the new callbacks
 *      SCIP_DECL_VARCOPY and SCIP_DECL_PROBCOPY led to new parameters in SCIPcreateVar() and SCIPcreateProb() in
 *      scip.c, respectively.
 *      <br><br>
 *    - SCIPincludeHeur() and SCIPincludeSepa() in \ref scip.h, as well as scip::ObjSepa() and scip::ObjHeur(), have a new parameter:
 *       - <code>usessubscip</code> - It can be used to inform SCIP that the heuristic/separator to be included uses a secondary SCIP instance.
 *      <br><br>
 *    - SCIPapplyRens() in \ref heur_rens.h has a new parameter <code>uselprows</code>. It can be used to switch from LP rows
 *      to constraints as basis of the sub-SCIP constructed in the RENS heuristic.
 *      <br>
 *      <br>
 *    - W.r.t. to copy and the C++ wrapper classes there are two new classes. These are <code>ObjCloneable</code> and
 *      <code>ObjProbCloneable</code>. The constraint handlers and variables pricers are derived from
 *      <code>ObjProbCloneable</code> and all other plugin are derived from <code>ObjCloneable</code>. Both
 *      classes implement the function <code>iscloneable()</code> which return whether a plugin is clone
 *      able or not. Besides that
 *      each class has a function named <code>clone()</code> which differ in their signature.
 *      See objcloneable.h, objprobcloneable.h, and the TSP example for more details.
 *
 * - <b>Branching</b>:
 *      <br><br>
 *    - The method SCIPgetVarStrongbranch() has been replaced by two methods SCIPgetVarStrongbranchFrac() and
 *      SCIPgetVarStrongbranchInt().
 *      <br><br>
 *    - The methods SCIPgetVarPseudocost() and SCIPgetVarPseudocostCurrentRun() in \ref scip.h now return the pseudocost value of
 *      one branching direction, scaled to a unit interval. The former versions of SCIPgetVarPseudocost() and
 *      SCIPgetVarPseudocostCurrentRun() are now called SCIPgetVarPseudocostVal() and SCIPgetVarPseudocostValCurrentRun(), respectively.
 *      <br>
 *      <br>
 *    - The methods SCIPgetVarConflictScore() and SCIPgetVarConflictScoreCurrentRun() in \ref scip.h are now called
 *      SCIPgetVarVSIDS() and SCIPgetVarVSIDSCurrentRun(), respectively.
 *      <br><br>
 *    - The methods SCIPvarGetNInferences(), SCIPvarGetNInferencesCurrentRun(), SCIPvarGetNCutoffs(), and
 *      SCIPvarGetNCutoffsCurrentRun() are now called SCIPvarGetInferenceSum(), SCIPvarGetInferenceSumCurrentRun(),
 *      SCIPvarGetCutoffSum(), and SCIPvarGetCutoffSumCurrentRun(), respectively. Furthermore, they now return
 *      <code>SCIP_Real</code> instead of <code>SCIP_Longint</code> values.
 *
 * - <b>Others</b>:
 *      <br><br>
 *    - SCIPcutGenerationHeuristicCmir() in \ref sepa_cmir.h has three new parameters:
 *        - <code>maxmksetcoefs</code> - If the mixed knapsack constraint obtained after aggregating LP rows contains more
 *          than <code>maxmksetcoefs</code> nonzero coefficients the generation of the <b>c-MIR cut</b> is aborted.
 *        - <code>delta</code> - It can be used to obtain the scaling factor which leads to the best c-MIR cut found within
 *          the cut generation heuristic. If a <code>NULL</code> pointer is passed, the corresponding c-MIR cut will already be
 *          added to SCIP by SCIPcutGenerationHeuristicCmir(). Otherwise, the user can generate the cut and add it to SCIP
 *          on demand afterwards.
 *        - <code>deltavalid</code> - In case, the user wants to know the best scaling factor, i.e., <code>delta</code> passed is not <code>NULL</code>,
 *          <code>deltavalid</code> will be <code>TRUE</code> if the stored scaling factor <code>delta</code> will lead to a violated c-MIR cut.
 *      <br>
 *      <br>
 *    - All functions for setting <b>user parameters</b> of different types like SCIPparamSetBool(), SCIPparamSetChar(),
 *      SCIPparamSetInt(), SCIPparamSetLongint(), and SCIPparamSetString() in pub_paramset.h have a new parameter:
 *        - <code>quiet</code> - It prevents any output during the assign to a new value.
 *
 * <br>
 * @section MISCELLANEOUS4 Miscellaneous
 *
 * - The NLPI library is now a separate library that is required when linking against the SCIP library.
 *   This requires changes to Makefiles that use SCIP, see the \ref RELEASENOTES "Release notes" for more details.
 *
 * - We do not distinguish between <b>block memory</b> for the original and the transformed problem anymore. The same
 *   block memory is now used in both problem stages.
 *
 * - The usage of <b>strong branching</b> changed. Now, SCIPstartStrongbranch() and SCIPendStrongbranch() must be
 *   called before and after strong branching, respectively.
 *
 * - All <b>C++</b> objects and constructors have a SCIP pointer, now.
 *
 * - The <b>predefined setting files</b> like "settings/cuts/off.set,aggressive.set,fast.set" have been replaced by
 *   interface methods like SCIPsetHeuristics(), SCIPsetPresolving(), SCIPsetSeparating(), and SCIPsetEmphasis() in
 *   \ref scip.h and by user dialogs in the interactive shell like
 *   <br>
 *   <br>
 *   <code>SCIP&gt; set {heuristics|presolving|separating} emphasis {aggressive|fast|off}</code>
 *   <br>
 *   <br>
 *   or
 *   <br>
 *   <br>
 *   <code>SCIP&gt; set emphasis {counter|cpsolver|easycip|feasibility|hardlp|optimality}</code>
 *
 *
 * <br>
 * For further release notes we refer the \ref RELEASENOTES "Release notes".
 */

/* - SCIP now has "lazy bounds", which are useful for column generation - see @ref PRICER_REMARKS "pricer remarks" for an explanation.
 *
 * - SCIP has rudimentary support to solve quadratic nonlinear integer programs - see \ref cons_quadratic.h.
 *
 * - There are LP-interfaces to QSopt and Gurobi (rudimentary).
 *
 * - SCIP can now handle indicator constraints (reading (from LP, ZIMPL), writing, solving, ...) - see \ref cons_indicator.h.
 *
 * - One can now do "early branching" useful for column generation.
 *
 * - Can now run a black-box lexicographic dual simplex algorithm.
 */

 /*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
 /**@page CHG5 Interface changes between SCIP 2.0 and SCIP 2.1
  *
  *
  * @section CHGCALLBACKS5 New and changed callbacks
  *
  * - <b>Presolving</b>:
  *      <br>
  *      <br>
  *    - The new parameters "nnewaddconss" and "naddconss" were added to the constraint handler callback method SCIP_DECL_CONSPRESOL()
  *      and the presolver callback method SCIP_DECL_PRESOLEXEC(). These parameters were also added to corresponding C++ wrapper class methods.
  *    - Propagators are now also called in during presolving, this is supported by the new callback methods SCIP_DECL_PROPINITPRE(),
  *      SCIP_DECL_PROPEXITPRE(), and SCIP_DECL_PROPPRESOL().
  *    - New parameters "isunbounded" and "isinfeasible" for presolving initialization (SCIP_DECL_CONSINITPRE(),
  *      SCIP_DECL_PRESOLINITPRE(), SCIP_DECL_PROPINITPRE()) and presolving deinitialization (SCIP_DECL_CONSEXITPRE(),
  *      SCIP_DECL_PRESOLEXITPRE(), SCIP_DECL_PROPEXITPRE()) callbacks of presolvers,
  *      constraint handlers and propagators, telling the callback whether the problem was already declared to be
  *      unbounded or infeasible.  This allows to avoid expensive steps in these methods in case the problem is already
  *      solved, anyway.
  *      <br>
  *      <br>
  *      <DIV class="note">
  *      Note, that the C++ methods
  *      - scip::ObjConshdlr::scip_presol() corresponding to SCIP_DECL_CONSPRESOL()
  *      - scip::ObjConshdlr::scip_initpre() corresponding to  SCIP_DECL_CONSINITPRE()
  *      - scip::ObjPresol::scip_initpre() corresponding to SCIP_DECL_PRESOLINITPRE()
  *      - scip::ObjProp::scip_initpre() corresponding to SCIP_DECL_PROPINITPRE()
  *      - scip::ObjConshdlr::scip_exitpre() corresponding to SCIP_DECL_CONSEXITPRE()
  *      - scip::ObjPresol::scip_exitpre() corresponding to SCIP_DECL_PRESOLEXITPRE()
  *      -  scip::ObjProp::scip_exitpre() corresponding to  and SCIP_DECL_PROPEXITPRE()
  *      .
  *      are virtual functions. That means, if you are not adding the new parameters, your code will still compile, but these methods are not executed.
  *      </DIV>
  *
  * - <b>Constraint Handler</b>:
  *     <br>
  *     <br>
  *   - The new constraint handler callback SCIP_DECL_CONSDELVARS() is called after variables were marked for deletion.
  *     This method is optional and only of interest if you are using SCIP as a branch-and-price framework. That means,
  *     you are generating new variables during the search. If you are not doing that just define the function pointer
  *     to be NULL.
  *     <br>
  *     If this method gets implemented you should iterate over all constraints of the constraint handler and delete all
  *     variables that were marked for deletion by SCIPdelVar().
  *
  * - <b>Problem Data</b>:
  *     <br>
  *     <br>
  *   - The method SCIPcopyProb() and the callback SCIP_DECL_PROBCOPY() got a new parameter "global" to indicate whether the global problem or a local version is copied.
  *
  * - <b>Conflict Analysis</b>:
  *     <br>
  *     <br>
  *   - Added parameter "separate" to conflict handler callback method SCIP_DECL_CONFLICTEXEC() that defines whether the conflict constraint should be separated or not.
  *
  * - <b>NLP Solver Interface</b>:
  *     <br>
  *     <br>
  *   - The callbacks SCIP_DECL_NLPIGETSOLUTION() and SCIP_DECL_NLPISETINITIALGUESS() got new parameters to get/set values of dual variables.
  *   - The callback SCIP_DECL_NLPICOPY() now passes the block memory of the target SCIP as an additional parameter.
  *
  * <br>
  * @section CHGINTERFUNC5 Changed interface methods
  *
  * - <b>Writing and Parsing constraints</b>:
  *      <br>
  *      <br>
  *    - The methods SCIPwriteVarName(), SCIPwriteVarsList(), and SCIPwriteVarsLinearsum() got a new boolean parameter "type"
  *      that indicates whether the variable type should be written or not.
  *    - The method SCIPwriteVarsList() got additionally a new parameter "delimiter" that defines the character which is used for delimitation.
  *    - The methods SCIPparseVarName() and SCIPparseVarsList() got a new output parameter "endptr" that is filled with the position where the parsing stopped.
  *
  * - <b>Plugin management</b>:
  *      <br>
  *      <br>
  *    - SCIPincludeProp() got additional parameters to set the timing mask of the propagator and the new callbacks and parameters related to calling the propagator in presolving.
  *    - SCIPincludeConshdlr() got additional parameters to set the variable deletion callback function and the timing mask for propagation.
  *
  * - <b>Constraint Handlers</b>:
  *      <br>
  *      <br>
  *    - Method SCIPseparateRelaxedKnapsack() in knapsack constraint handler got new parameter "cutoff", which is a pointer to store whether a cutoff was found.
  *    - Method SCIPincludeQuadconsUpgrade() of quadratic constraint handler got new parameter "active" to indicate whether the upgrading method is active by default.
  *
  * - <b>Nonlinear expressions, relaxation, and solver interface</b>:
  *      <br>
  *      <br>
  *    - The methods SCIPexprtreeEvalSol(), SCIPexprtreeEvalIntLocalBounds(), and SCIPexprtreeEvalIntGlobalBounds() have been renamed to SCIPevalExprtreeSol(),
  *      SCIPevalExprtreeLocalBounds(), and SCIPevalExprtreeGlobalBounds() and are now located in scip.h.
  *    - Various types and functions dealing with polynomial expressions have been renamed to use the proper terms "monomial" and "polynomial".
  *    - The methods SCIPnlpGetObjective(), SCIPnlpGetSolVals(), and SCIPnlpGetVarSolVal() have been removed, use SCIPgetNLPObjval(), SCIPvarGetNLPSol()
  *      and SCIPcreateNLPSol() to retrieve NLP solution values instead.
  *    - Removed methods SCIPmarkRequireNLP() and SCIPisNLPRequired(), because the NLP is now always constructed if nonlinearities are present.
  *    - SCIPgetNLP() has been removed and NLP-methods from pub_nlp.h have been moved to scip.h, which resulted in some renamings, too.
  *    - The functions SCIPnlpiGetSolution() and SCIPnlpiSetInitialGuess() got additional arguments to get/set dual values.
  *    - The method SCIPgetNLPI() got a new parameter "nlpiproblem", which is a pointer to store the NLP solver interface problem.
  *
  * - <b>Others</b>:
  *      <br>
  *      <br>
  *    - SCIPgetVarCopy() got a new parameter "success" that will be FALSE if method is called after problem creation stage and no hash map is given or no image for
  *      the given variable is contained in the given hash map.
  *    - Removed method SCIPreadSol(); call solution reading via SCIPreadProb() which calls the solution reader for .sol files.
  *    - SCIPchgVarType() got an extra boolean parameter to store if infeasibility is recognized while upgrading a variable from continuous type to an integer type.
  *    - SCIPdelVar() got a new parameter "deleted", which stores whether the variable was successfully marked to be deleted.
  *    - SCIPcalcNodeselPriority() got a new parameter "branchdir", which defines the type of branching that was performed: upwards, downwards, or fixed.
  *    - The parameters "timelimit" and "memorylimit" were removed from SCIPapplyRens().
  *
  * <br>
  * @section MISCELLANEOUS5 Miscellaneous
  *
  *  - The result value SCIP_NEWROUND has been added, it allows a separator/constraint handler to start a new separation round
  *    (without previous calls to other separators/conshdlrs).
  *  - All timing flags are now defined type_timing.h.
  *  - The variable deletion event is now a variable specific event and not global, anymore.
  *  - The emphasis setting types now distinguish between plugin-type specific parameter settings (default, aggressive, fast, off), which are changed by
  *    SCIPsetHeuristics/Presolving/Separating(), and global emphasis settings (default, cpsolver, easycip, feasibility, hardlp, optimality, counter),
  *    which can be set using SCIPsetEmphasis().
  *
  * <br>
  * For further release notes we refer the \ref RELEASENOTES "Release notes".
  */

 /*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
 /**@page CHG6 Interface changes between SCIP 2.1 and SCIP 3.0
  *
  *
  * @section CHGCALLBACKS6 New and changed callbacks
  *
  * - <b>Conflict Analysis</b>:
  *     <br>
  *     <br>
  *   - Added parameter "relaxedbds" to conflict handler callback method SCIP_DECL_CONFLICTEXEC(). This array contains
  *     bounds which are sufficient to create a valid conflict
  *
  * - <b>Constraint Handler</b>:
  *     <br>
  *     <br>
  *   - New optional callback methods in constraint handlers: SCIP_DECL_CONSGETVARS and SCIP_DECL_CONSGETNVARS.
  *     These callbacks, if implemented, should return an array of all variables and the number of all variables used
  *     by the given constraint, respectively. (This method might, e.g., be called by a presolver)
  *   - Added a propagation timing parameter "proptiming" to SCIP_DECL_CONSPROP(), giving the current timing at which
  *     this method is called
  *   - Added a parameter 'restart' to the SCIP_DECL_CONSEXITSOL() callback method, indicating whether this call was
  *     triggered by a restart.
  *   - Added a parameter 'relaxedbd' to SCIP_DECL_CONSRESPROP() callback method. If explaining a given bound change
  *     (index), it is sufficient to explain the reason for reaching the 'relaxedbd' value, see above
  *   - Removed parameters "isunbounded", "isinfeasible" and "result" from SCIP_DECL_CONSINITPRE() and
  *     SCIP_DECL_CONSEXITPRE() callback methods. It is not allowed to determine unboundedness or infeasibility in
  *     these callbacks, anymore.
  *
  * - <b>Message Handler</b>:
  *      <br>
  *      <br>
  *   - New callback method SCIP_DECL_MESSAGEHDLRFREE() which is called when the message handler is freed.
  *   - The old callback method SCIP_DECL_MESSAGEERROR() was replaced by the callback method SCIP_DECL_ERRORPRINTING().
  *
  * - <b>Presolving</b>:
  *      <br>
  *      <br>
  *   - Removed parameters "isunbounded", "isinfeasible" and "result" from SCIP_DECL_PRESOLINITPRE() and
  *     SCIP_DECL_PRESOLSEXITPRE(). It is not allowed to determine unboundedness or infeasibility in these
  *     callbacks, anymore.
  *
  * - <b>Propagator</b>:
  *     <br>
  *     <br>
  *   - Added a propagation timing parameter "proptiming" to SCIP_DECL_PROPEXEC(), giving the
  *     current timing at which this method is called.
  *   - Added a parameter 'restart' to SCIP_DECL_PROPEXITSOL() callback method, indicating whether this call was
  *     triggered by a restart.
  *   - Added a parameter 'relaxedbd' to SCIP_DECL_PROPRESPROP() callback method. If explaining a given bound change
  *     (index), it is sufficient to explain the reason for reaching the 'relaxedbd' value.
  *   - Removed parameters "isunbounded", "isinfeasible" and "result" from SCIP_DECL_PROPINITPRE() and
  *     SCIP_DECL_PROPEXITPRE() callback methods. It is not allowed to determined unboundedness or infeasibility in
  *     these callbacks, anymore.
  *
  * - <b>NLP Solver Interface</b>:
  *     <br>
  *     <br>
  *   - New NLPI callback SCIP_DECL_NLPISETMESSAGEHDLR() to set message handler in NLP solver interfaces.
  *
  * <br>
  * @section CHGINTERFUNC6 Changed interface methods
  *
  * - <b>Plugin management</b>:
  *      <br>
  *      <br>
  *   - Added basic include methods for almost all plugin types, e.g., SCIPincludeConshdlrBasic();
  *     these methods should make the usage easier, sparing out optional callbacks and parameters.
  *   - To extend the basic functionalities, there are setter method to add
  *     optional callbacks. For example SCIPsetConshdlrParse(), SCIPsetPropCopy() or SCIPsetHeurInitsol().
  *
  * - <b>Constraint Handlers</b>:
  *      <br>
  *      <br>
  *   - Added basic creation methods for all constraints types, e.g., SCIPcreateConsBasicLinear(); these methods should make the usage easier,
  *      sparing out optional callbacks and parameters.
  *   - New methods SCIPgetConsVars() and SCIPgetConsNVars() (corresponding callbacks need to be implemented, see above)
  *
  * - <b>Problem</b>:
  *      <br>
  *      <br>
  *   - Added basic creation methods SCIPcreateVarBasic() and SCIPcreateProbBasic() and setter functions
  *   - Added method SCIPisPresolveFinished() which returns whether the presolving process would be stopped after the
  *     current presolving round, given no further reductions will be found.
  *   - Forbid problem modifications in SCIP_STAGE_{INIT,EXIT}PRESOLVE (see pre-conditions for corresponding methods in scip.h).
  *
  * - <b>Variable usage</b>:
  *      <br>
  *      <br>
  *   - Renamed SCIPvarGetBestBound() to SCIPvarGetBestBoundLocal(), SCIPvarGetWorstBound() to
  *     SCIPvarGetWorstBoundLocal() and added new methods SCIPvarGetBestBoundGlobal() and SCIPvarGetWorstBoundGlobal().
  *   - Method SCIPvarGetProbvarSum() is not public anymore, use SCIPgetProbvarSum() instead.
  *   - Replaced method SCIPvarGetRootRedcost() by SCIPvarGetBestRootRedcost().
  *
  * - <b>Message Handler</b>:
  *      <br>
  *      <br>
  *   - Changed the message handler system within SCIP heavily such that it is thread-safe. SCIPcreateMessagehdlr() in
  *     scip.{c,h} was replaced by SCIPmessagehdlrCreate() in pub_message.h/message.c with a changed parameter list.
  *   - Error messages (SCIPerrorMessage()) are not handled via the message handler anymore; per default the error
  *     message is written to stderr.
  *
  * - <b>Separation</b>:
  *      <br>
  *      <br>
  *   - New functions SCIPcreateEmptyRowCons(), SCIPcreateEmptyRowSepa(), SCIPcreateRowCons(), and SCIPcreateRowSepa()
  *     that allow to set the originating constraint handler or separator of a row respectively; this is, for instance,
  *     needed for statistics on the number of applied cuts. If rows are created outside a constraint handler or separator
  *     use SCIPcreateRowUnspec() and SCIPcreateEmptyRowUnspec(). The use of SCIPcreateEmptyRow() and SCIPcreateRow() is
  *     deprecated.
  *   - New functions SCIProwGetOrigintype(), SCIProwGetOriginCons(), and SCIProwGetOriginSepa() to obtain the originator
  *     that created a row.
  *
  * - <b>LP interface</b>:
  *      <br>
  *      <br>
  *   - SCIPlpiCreate() got a new parameter 'messagehdlr'.
  *   - SoPlex LPI supports setting of SCIP_LPPAR_DUALFEASTOL when using SoPlex version 1.6.0.5 and higher.
  *
  * - <b>Nonlinear expressions, relaxation, and solver interface</b>:
  *      <br>
  *      <br>
  *   - Renamed SCIPmarkNonlinearitiesPresent() to SCIPenableNLP() and SCIPhasNonlinearitiesPresent() to
  *     SCIPisNLPEnabled().
  *   - Method SCIPexprtreeRemoveFixedVars() is not public anymore.
  *
  * - <b>Counting</b>:
  *      <br>
  *      <br>
  *   - Changed the counting system within SCIP heavily. SPARSESOLUTION was renamed to SCIP_SPARSESOL. New method for
  *     SCIP_SPARSESOL usage, SCIPsparseSolCreate(), SCIPsparseSolFree(), SCIPsparseSolGetVars(),
  *     SCIPsparseSolGetNVars(), SCIPsparseSolGetLbs(), SCIPsparseSolGetUbs() in (pub_)misc.{c,h}.
  *   - Renamed SCIPgetCountedSparseSolutions() to SCIPgetCountedSparseSols() in cons_countsols.{c,h}.
  *
  * <br>
  * @section MISCELLANEOUS6 Miscellaneous
  *
  *   - Replaced SCIPparamSet*() by SCIPchg*Param() (where * is either Bool, Int, Longint, Real, Char, or String).
  *   - New parameter in SCIPcopy() and SCIPcopyPlugins() to indicate whether the message handler from the source SCIP
  *     should be passed to the target SCIP (only the pointer is copied and the usage counter of the message handler is
  *     increased).
  *   - SCIPprintCons() does not print termination symbol ";\n" anymore; if wanted, use SCIPinfoMessage() to print
  *     ";\n" manually
  *   - All objscip *.h file now use the default SCIP interface macros.
  *   - The methods SCIPsortedvecInsert*() have an additional parameter which can be used to receive the position where
  *     the new element was inserted.
  *   - New macro SCIPdebugPrintCons() to print constraint only if SCIP_DEBUG flag is set.
  *
  * <br>
  * For further information we refer to the \ref RELEASENOTES "Release notes" and the \ref CHANGELOG "Changelog".
  */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page COUNTER How to use SCIP to count/enumerate feasible solutions
 *
 * SCIP is capable of computing (count or enumerate) the number of feasible solutions of a given constraint integer
 * program. In case continuous variables are present, the number of feasible solutions for the projection to the
 * integral variables is counted/enumerated. This means, an assignment to the integer variables is counted if the
 * remaining problem (this is the one after fixing the integer variables w.r.t. to this assignment) is feasible.
 *
 * As a first step you have to load or create your problem in the usual way. In case of using the
 * interactive shell, you use the <code>read</code> command:
 *
 * <code>SCIP&gt; read &lt;file name&gt;</code>
 *
 * Afterwards you can count the number of feasible solution with the command <code>count</code>.
 *
 * <code>SCIP&gt; count</code>
 *
 * That means SCIP will count the number of solution but does not store (enumerate) them. If you are interested in that see
 * \ref COLLECTALLFEASEBLES.
 *
 * @note Since SCIP version 2.0.0 you do not have to worry about <tt>dual</tt> reductions anymore. These are
 * automatically turned off. The only thing you should switch off are restarts. These restarts can lead to a wrong
 * counting process. We recommend using the counting settings which can be set in the interactive shell as follows:
 *
 * <code>SCIP&gt; set emphasis counter</code>
 *
 * The SCIP callable library provides an interface method SCIPcount() which allows users to count the number of feasible
 * solutions to their problem. The method SCIPsetParamsCountsols(), which is also located in cons_countsols.h, loads the
 * predefined counting settings to ensure a safe count. The complete list of all methods that can be used for counting
 * via the callable library can be found in cons_countsols.h.
 *
 *
 * @section COUNTLIMIT Limit the number of solutions which should be counted
 *
 * It is possible to give a (soft) upper bound on the number solutions that should be counted. If this upper bound is
 * exceeded, SCIP will be stopped. The name of this parameter is <code>constraints/countsols/sollimit</code>. In
 * the interactive shell this parameter can be set as follows:
 *
 * <code>SCIP&gt; set constraints countsols sollimit 1000</code>
 *
 * In case you are using the callable library, this upper bound can be assigned by calling SCIPsetLongintParam() as follows:
 *
 * \code
 * SCIP_CALL( SCIPsetLongintParam( scip, "constraints/countsols/sollimit", 1000) );
 * \endcode
 *
 *
 * The reason why this upper bound is soft comes from the fact that, by default, SCIP uses a technique called unrestricted
 * subtree detection. Using this technique it is possible to detect several solutions at once. Therefore, it can happen
 * that the solution limit is exceeded before SCIP is stopped.
 *
 * @section COLLECTALLFEASEBLES Collect all feasible solution
 *
 * Per default SCIP only counts all feasible solutions. This means, these solutions are not stored. If you switch the
 * parameter <code>constraints/countsols/collect</code> to TRUE (the default value is FALSE) the detected solutions are
 * stored. Changing this parameter can be done in the interactive shell
 *
 * <code>SCIP&gt; set constraints countsols collect TRUE</code>
 *
 * as well as via the callable library
 *
 * \code
 * SCIP_CALL( SCIPsetBoolParam( scip, "constraints/countsols/collect", TRUE) );
 * \endcode
 *
 * @note The solution which are collected are stored w.r.t. the active variables. These are the variables which got not
 *       removed during presolving.
 *
 * In case you are using the interactive shell you can write all collected solutions to a file as follows
 *
 * <code>SCIP&gt; write allsolutions &lt;file name&gt;</code>
 *
 * In that case the sparse solutions are unrolled and lifted back into the original variable space.
 *
 * The callable library provides a method which gives access to all collected sparse solutions. That is,
 * SCIPgetCountedSparseSolutions(). The sparse solutions you get are defined w.r.t. active variables. To get solutions
 * w.r.t. to the original variables. You have to do two things:
 *
 * -# unroll each sparse solution
 * -# lift each solution into original variable space by extending the solution by those variable which got removed
 *    during presolving
 *
 * The get the variables which got removed during presolving, you can use the methods SCIPgetFixedVars() and
 * SCIPgetNFixedVars(). The method SCIPgetProbvarLinearSum() transforms given variables, scalars and constant to the
 * corresponding active variables, scalars and constant. Using this method for a single variable gives a representation
 * for that variable w.r.t. the active variables which can be used to compute the value for the considered solution (which
 * is defined w.r.t. active variables).
 *
 * For that complete procedure you can also check the source code of
 * \ref SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteAllsolutions) "SCIPdialogExecWriteAllsolutions()" cons_countsols.c which
 * does exactly that.
 *
 *
 * @section COUNTOPTIMAL Count number of optimal solutions
 *
 * If you are interested in counting the number of optimal solutions, this can be done with SCIP using the
 * <code>count</code> command by applying the following steps:
 *
 *  -# Solve the original problem to optimality and let \f$c^*\f$ be the optimal value
 *  -# Add the objective function as constraint with left and right hand side equal to \f$c^*\f$
 *  -# load the adjusted problem into SCIP
 *  -# use the predefined counting settings
 *  -# start counting the number of feasible solutions
 *
 * If you do this, SCIP will collect all optimal solutions of the original problem.
 *
 */

/**@page LICENSE License
 *
 * \verbinclude COPYING
 */

/**@page FAQ Frequently Asked Questions (FAQ)
 * \htmlinclude faqcss.inc
 * \htmlinclude faq.inc
 */


/**@page AUTHORS SCIP Authors
 * \htmlinclude authors.inc
 */

/**@page INSTALL Installation information
 * \verbinclude INSTALL
 */

/**@page RELEASENOTES Release notes
 *
 * \verbinclude SCIP-release-notes-3.0.1
 *
 * \verbinclude SCIP-release-notes-3.0
 *
 * \verbinclude SCIP-release-notes-2.1.1
 *
 * \verbinclude SCIP-release-notes-2.1
 *
 * \verbinclude SCIP-release-notes-2.0.2
 *
 * \verbinclude SCIP-release-notes-2.0.1
 *
 * \verbinclude SCIP-release-notes-2.0
 *
 * \verbinclude SCIP-release-notes-1.2
 *
 * \verbinclude SCIP-release-notes-1.1
 */

/**@page CHANGELOG CHANGELOG
 *
 * \verbinclude CHANGELOG
 *
 */

/**@defgroup PUBLICMETHODS Public Methods
 *
 * This page lists headers containing methods provided by the core of SCIP that can be used via the
 * callable library. If you are in the <a href="../html/index.shtml">User's Manual</a> you only find methods that are
 * public and, therefore, allowed to be used. The <a href="../html_devel/index.shtml">Developer's Manual</a> includes
 * all methods.
 *
 * All of the headers listed below include functions that are allowed to be called by external users. Besides those
 * functions it is also valid to call methods that are listed in one of the headers of the (default) plugins, e.g.,
 * cons_linear.h.
 *
 * If you are looking for information about a particular object of SCIP, such as a variable or a constraint, you should
 * first search the corresponding "pub_<...>.h" header. E.g., for constraints, look in pub_cons.h. If you need some
 * information about the overall problem, you should start searching in scip.h.
 *
 * Since there is a huge number of methods in scip.h, these methods are grouped into different categories. These
 * categories are:
 *
 * - Memory Management
 * - Miscellaneous Methods
 * - General SCIP Methods
 * - Message Output Methods
 * - Parameter Methods
 * - SCIP User Functionality Methods: Managing Plugins
 * - User Interactive Dialog Methods
 * - Global Problem Methods
 * - Local Subproblem Methods
 * - Solve Methods
 * - Variable Methods
 * - Conflict Analysis Methods
 * - Constraint Methods
 * - LP Methods
 * - LP Column Methods
 * - LP Row Methods
 * - Cutting Plane Methods
 * - LP Diving Methods
 * - Probing Methods
 * - Branching Methods
 * - Primal Solution Methods
 * - Event Methods
 * - Tree Methods
 * - Statistic Methods
 * - Timing Methods
 * - Numerical Methods
 * - Dynamic Arrays
 *
 */

/**@defgroup TYPEDEFINITIONS Type Definitions
 * This page lists headers which contain type definitions of callback methods.
 *
 * All headers below include the descriptions of callback methods of
 * certain plugins. For more detail see the corresponding header.
 */

/**@defgroup BRANCHINGRULES Branching Rules
 * @brief This page contains a list of all branching rule which are currently available.
 *
 * A detailed description what a branching rule does and how to add a branching rule to SCIP can be found
 * \ref BRANCH "here".
 */

/**@defgroup CONSHDLRS  Constraint Handler
 * @brief This page contains a list of all constraint handlers which are currently available.
 *
 * A detailed description what a constraint handler does and how to add a constraint handler to SCIP can be found
 * \ref CONS "here".
 */

/**@defgroup DIALOGS Dialogs
 * @brief This page contains a list of all dialogs which are currently available.
 *
 * A detailed description what a dialog does and how to add a dialog to SCIP can be found
 * \ref DIALOG "here".
 */

/**@defgroup DISPLAYS Displays
 * @brief This page contains a list of all displays (output columns)  which are currently available.
 *
 * A detailed description what a display does and how to add a display to SCIP can be found
 * \ref DISP "here".
 *
 */

/**@defgroup EXPRINTS Expression Interpreter
 * @brief This page contains a list of all expression interpreter which are currently available.
 *
 * A detailed description what a expression interpreter does and how to add a expression interpreter to SCIP can be found
 * \ref EXPRINT "here".
 */

/**@defgroup FILEREADERS File Readers
 * @brief This page contains a list of all file readers which are currently available.
 *
 * @section AVAILABLEFORMATS List of readable file formats
 *
 * The \ref SHELL "interactive shell" and the callable library are capable of reading/parsing several different file
 * formats.
 *
 * <table>
 * <tr><td>\ref reader_cip.h "CIP format"</td> <td>for SCIP's constraint integer programming format</td></tr>
 * <tr><td>\ref reader_cnf.h "CNF format"</td> <td>DIMACS CNF (conjunctive normal form) file format used for example for SAT problems</td></tr>
 * <tr><td>\ref reader_fzn.h "FZN format"</td> <td>FlatZinc is a low-level solver input language that is the target language for MiniZinc.</td></tr>
 * <tr><td>\ref reader_gms.h "GMS format"</td> <td>for mixed-integer nonlinear programs (<a href="http://www.gams.com/docs/document.htm">GAMS</a>) [write only]</td></tr>
 * <tr><td>\ref reader_lp.h  "LP format"</td>  <td>for mixed-integer (quadratically constrained quadratic) programs (CPLEX)</td></tr>
 * <tr><td>\ref reader_mps.h "MPS format"</td> <td>for mixed-integer (quadratically constrained quadratic) programs</td></tr>
 * <tr><td>\ref reader_opb.h "OPB format"</td> <td>for pseudo-Boolean optimization instances</td></tr>
 * <tr><td>\ref reader_osil.h "OSiL format"</td> <td>for mixed-integer nonlinear programs</td></tr>
 * <tr><td>\ref reader_pip.h "PIP format"</td> <td>for <a href="http://polip.zib.de/pipformat.php">mixed-integer polynomial programming problems</a></td></tr>
 * <tr><td>\ref reader_sol.h "SOL format"</td> <td>for solutions; XML-format (read-only) or raw SCIP format</td></tr>
 * <tr><td>\ref reader_wbo.h "WBO format"</td> <td>for weighted pseudo-Boolean optimization instances</td></tr>
 * <tr><td>\ref reader_zpl.h "ZPL format"</td> <td>for <a href="http://zimpl.zib.de">ZIMPL</a> models, i.e., mixed-integer linear and nonlinear
 *                                                 programming problems [read only]</td></tr>
 * </table>
 *
 * @section ADDREADER How to add a file reader
 *
 * A detailed description what a file reader does and how to add a file reader to SCIP can be found
 * \ref READER "here".
 *
 */

/**@defgroup LPIS LP Solver Interfaces
 * @brief This page contains a list of all LP solver interfaces which are currently available.
 */

/**@defgroup NODESELECTORS Node Selectors
 * @brief This page contains a list of all node selectors which are currently available.
 *
 * A detailed description what a node selector does and how to add a node selector to SCIP can be found
 * \ref NODESEL "here".
 */

/**@defgroup NLPIS NLP Solver Interfaces
 * @brief This page contains a list of all NLP solver interfaces which are currently available.
 *
 * A detailed description what a NLP solver interface does and how to add a NLP solver interface to SCIP can be found
 * \ref NLPI "here".
 */

/**@defgroup PRESOLVERS Presolvers
 * @brief This page contains a list of all presolvers which are currently available.
 *
 * A detailed description what a presolver does and how to add a presolver to SCIP can be found
 * \ref PRESOL "here".
 */

/**@defgroup PRICERS Pricers
 * @brief This page contains a list of all pricers which are currently available.
 *
 * Per default there exist no variable pricer. A detailed description what a variable pricer does and how to add a
 * variable pricer to SCIP can be found \ref PRICER "here".
 */

/**@defgroup PRIMALHEURISTICS Primal Heuristics
 * @brief This page contains a list of all primal heuristics which are currently available.
 *
 * A detailed description what a primal heuristic does and how to add a primal heuristic to SCIP can be found
 * \ref HEUR "here".
 */

/**@defgroup PROPAGATORS Propagators
 * @brief This page contains a list of all propagators which are currently available.
 *
 * A detailed description what a propagator does and how to add a propagator to SCIP can be found
 * \ref PROP "here".
 */

/**@defgroup RELAXATORS Relaxation Handlers
 * @brief This page contains a list of all relaxation handlers which are currently available.
 *
 * Note that the linear programming relaxation is not implemented via the relaxation handler plugin. Per default there
 * exist no relaxation handler. A detailed description what a variable pricer does and how to add a A detailed
 * description what a relaxation handler does and how to add a relaxation handler to SCIP can be found \ref RELAX
 * "here".
 */

/**@defgroup SEPARATORS Separators
 * @brief This page contains a list of all separators  which are currently available.
 *
 * A detailed description what a separator does and how to add a separator to SCIP can be found
 * \ref SEPA "here".
 */

/**@page PARAMETERS List of all SCIP parameters
 *
 * This page list all parameters of the current SCIP version. This list can
 * easily be generated by SCIP via the interactive shell using the following command:
 *
 * <code>SCIP&gt; set save &lt;file name&gt;</code>
 *
 * or via the function call:
 *
 * <code>SCIP_CALL( SCIPwriteParams(scip, &lt;file name&gt;, TRUE, FALSE) );</code>
 *
 * \verbinclude parameters.set
 */

