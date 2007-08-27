<div id="faq">

    <h4><a name="Section0">Sections</a></h4>
    <ul>
     <li><a href="#Section1">General Questions about SCIP</a></li>
     <li><a href="#Section2">Using SCIP as a standalone MIP-Solver</a></li>
     <li><a href="#Section3">Using SCIP included in another source code</a></li>
     <li><a href="#Section4">Using SCIP as a Branch-Cut-And-Price-Framework</a></li>
    </ul>

    <h4><a name="Section1">General Questions about SCIP</a></h4>
    <p class="q">
      <a name="Q1">What is SCIP?</a>
    </p>
    <p class="a">
      SCIP is a standalone MIP-Solver as well as a branch-cut-and-price-framework.<br> 
      You can use the precompiled binaries with which you can solve MIPs in *.mps or *.lp file format,
      you can SCIP as a subroutine for solving MIPs from your own source code, 
      or you can use it as a framework in which you implement your own plugins&mdash;or for a combination
      of these three purposes.<br>
      In this FAQ there are separate sections for the three usages of SCIP.
    <p class="q">
      <a name="Q2">When should I use SCIP?</a>
    </p>
    <p class="a">
      If you are either looking for a fast non-commercial MIP-Solver or for a 
      branch-cut-and-price-framework in which you can directly implement your own methods 
      &mdash; which can even be for more general purposes than MIP (see <a href="#Q18">here</a>).
    </p>
    <p class="q">
      <a name="Q3">I heard something about licenses. Do I have to pay for using SCIP?</a>
    </p>
    <p class="a"> 
      As long as you use it for academical, non-commercial purposes: No. 
      This will not change. For the other cases, check the explanation of the 
      <a href="http://zibopt.zib.de/ZIBLicence.shtml">ZIB academic license</a>.
    </p>
    <p class="q">
      <a name="Q4">Do I need any extra software?</a>
    </p>
    <p class="a">
      Except if you want to use SCIP as a pure CP-Solver (see <a href="#Q18">here</a>), 
      you need an underlying LP-Solver which has first to be installed and linked to the libraries 
      (see the &quot;INSTALL&quot; file in the SCIP root directory). 
      LP-solvers currently supported by SCIP are: 
    </p>
    <ul>
      <li><a href="http://www.zib.de/berthold/soplex/">SoPlex</a></li>
      <li><a href="http://www.ilog.com/products/cplex/">ILOG CPLEX</a></li>
      <li><a href="http://www.dashoptimization.com/home/products/products_optimizer.html">Dash XPress-MP</a></li>
      <li><a href="http://www.coin-or.org/Clp/">CLP</a></li>
      <li><a href="http://www.mosek.com/">Mosek</a></li>
    </ul>
    <p class="a">
      You also might need a modeling language like <a href="http://www.zib.de/koch/zimpl/">ZIMPL</a> to generate *.mps 
    or *.lp files. ZIMPL files can also directly be read by SCIP.
    </p>
    <p class="q">
      <a name="Q5">How do I get started?</a>
    </p>
    <p class="a">
      An easy way is to use the SCIP-binaries and call SCIP from a shell. 
      For that, you just have to download one of the precompiled binaries from the 
      <a href="http://scip.zib.de/download.shtml">download section</a> or the zipped source code and compile it 
      with your favorite settings.
    </p>
    <p class="a">
      Another way is to use SCIP as a solver integrated into your own program source code. 
      See the directory &quot;examples/MIPsolver/&quot; for a simple example and 
      <a href="#Q12">this point</a>.
    </p>
    <p class="a">
      A third way is to implement your own plugins into SCIP. 
      This is explained in the HowTos of all plugin types, which you can find at the 
      <a href="http://scip.zib.de/doc/html/index.html">doxygen documentation</a>.
    </p>
    <p class="q">
      <a name="Q6"> I would like to check whether some functionality is implemented into SCIP. 
	How does the naming of the methods work? Where do I find the most common methods?</a>
    </p>
    <p class="a">
      For an explanation of the naming see the 
      <a href="http://scip.zib.de/doc/html/CODE.html">coding style guidelines</a>. <br> 
      The methods being of interest for you, should normally be found either in scip.h or in pub_*.h, 
      see also <a href="http://scip.zib.de/doc/html/DOC.html">here</a>.
      In the doxygen documentation of <a href="http://scip.zib.de/doc/html/scip_8h.html">scip.h</a>, 
      the methods are ordered by topics.		      
    </p>
    <p class="q">
      <a name="Q7">The methods SCIPgetVarSol() and SCIPvarGetSol() seem to have the same
	functionality. Which one should I use?</a>
    </p>
    <p class="a">
      In fact, there is a slight difference: SCIPvarGetSol() is also able to return pseudo solution values.
      If you do not have an idea, what pseudo solutions are, SCIPgetVarSol() should be just fine.<br>
      This should be the only case of 'duplicate methods'. 
      However, if you find another one, please contact us.
    </p>
     <p class="q">
       <a name="Q27">I have installation problems. What can I do?</a>
    </p>
    <p class="a">
      If you want to use <a href="http://soplex.zib.de/">SoPlex</a> as the underlying LP-solver,
      try the following:                   
      First, download the <a href="http://zibopt.zib.de">ZIB Optimization Suite</a>.
      Then, extract the file, change into the ziboptsuite directory, and enter 'make'.
      As long as you have all the necessary libraries installed in your system, it should generate a SCIP
      binary linked to ZIMPL and Soplex. The necessary system libraries are:
    </p>
    <ol>
        <li>ZLIB (libz.a)</li>
        <li>GMP (libgmp.a)</li>
        <li>Readline (libreadline.a)</li>
    </ol>
    <p class="a">
      If you do not have all of these libraries, read the INSTALL file in the ZIB Optimization
      Suite directory.<br>
      As a summary, the call of 'make ZIMPL=false ZLIB=false READLINE=false' should work on most systems.
    </p>
    <p class="a">
      If you have any problems while use an LP-solver different from SoPlex,
      please read the SCIP INSTALL file first.                                                     
    </p>
    <p class="q">
      <a name="Q8">I changed to version 1.0 and now compiling breaks with some error like:<br>
	&quot;scip/src/scip/scip.h:343: error: too few arguments to function 'SCIP_RETCODE 
	SCIPaddIntParam(SCIP*, const char*, const char*, int*, unsigned int, int, int, int, 
	SCIP_RETCODE (*)(SCIP*, SCIP_PARAM*), SCIP_PARAMDATA*)'&quot;<br>
	What's going on?</a>
    </p>
    <p class="a">
      The SCIPadd&lt;datatype&gt;Param() functions got a new parameter flag 
      &quot;isadvanced&quot;. This has no influence to the performance, but only to the 
      position of the parameter in the settings menu. Hence, if you do not care about this,
      you can assign any value to it.<br>
      You should add the corresponding flag to the SCIPadd&lt;datatype&gt;Param() calls in your
      own source code. 
    </p>
    <p class="q">
      <a name="Q18">Can I use SCIP as a pure CP-Solver?</a>
    </p>
    <p class="a">
      Yes. There is a setting-file settings/emphasis/cpsolver.set. 
      Furthermore, you can compile SCIP without any LP-Solver by &quot;make LPS=none&quot;.
    </p>
<h4><a name="Section2">Using SCIP as a standalone MIP-Solver</a></h4>
    <p class="q">
      <a name="Q22">The output is too wide for my terminal window. What can I do?</a>
    </p>
    <p class="a">
      Type &quot;set display width&quot; followed by an appropriate number into the interactive shell.
      See also the <a href="#Q9">next question</a>.		      
    </p>
    <p class="q">
      <a name="Q9">What do the cryptic abbreviations for the columns mean which are displayed 
	during the solving process of SCIP?</a>
    </p>
    <p class="a">
      Type "display display" in the interactive shell to get an explanation of them. 
      By the way: If a letter appears in front of a display row, it indicates, 
      which heuristic found the new primal bound, a star representing an integral LP-relaxation.<br>
      Typing "display statistics" after finishing or interrupting the solving process gives 
      you plenty of extra information about the solving process.
    </p>
    <p class="q">
      <a name="Q10">How do I change the behavior of SCIP?</a>
    </p>
    <p class="a">
      You can find a variety of different general settings for SCIP in the directory "settings". 
      They can be loaded via "set load settings/&#42;/&#42;.set". 
      You can combine the general settings for cuts, presolving, and heuristics arbitrarily.<br>
      "display parameters" shows you which settings currently differ from their default, 
      "set default" resets them all. 
      Furthermore, there are complete settings in settings/emphasis/, analogously to 
      <a href="http://www.ilog.com/products/cplex/">ILOG CPLEX</a>.		    
    </p>
    <p class="q">
      <a name="Q11"> I recognized that one special plugin works very bad / very well for my problem and 
	I want to disable it / weaken its influence / intensify its influence. How do I do this?</a>
    </p>
    <p class="a">
      For using a <b class="grey">non-default branching rule or node selection strategy</b> as standard, 
      you just have to give it the highest priority, using, e.g., 
      &quot;set branching mostinf priority 9999999&quot;.<br>
      If you want to completely <b class="grey">disable a heuristic or a separator</b> you have to set its 
      frequency to -1 and the sepafreq to -1 for separation by constraint handlers, respectively.<br>
      For <b class="grey">disabling a presolver</b>, you have to set its maxrounds parameter to 0.<br>
      If you want to <b class="grey">intensify</b> the usage of a <b class="grey">heuristic</b>, 
      you can reduce its frequency to some smaller, positive value, and/or raise the quotient and 
      offset values (maxlpiter for diving heuristics, nodes for LNS heuristics).<br>
      For <b class="grey">intensifying</b> the usage of a <b class="grey">separator</b>, 
      you can raise its maxroundsroot and maxsepacutsroot values. 
      If you also want to use this separator locally, you have to set its frequency to a positive value 
      and possibly raise maxrounds and maxsepacuts.
      Compare the parameters of the heuristic/separator in the appropriate aggressive setting 
      (see <a href="#Q10">previous question</a>).<br>
      For <b class="grey">weakening</b>, you should just do the opposite operation, i.e. reducing the values 
      you would raise for intensification and vice versa.
    </p>
  <h4><a name="Section3">Using SCIP included in another source code</a></h4>
    <p class="q">
      <a name="Q12">How Do I construct a problem instance in SCIP?</a>
    </p>
    <p class="a">
      First you have to create a SCIP object via SCIPcreate(), 
      then you start to build the problem via SCIPcreateProb().<br>
      Then you create variables via SCIPcreateVar() and add them to the problem via SCIPaddVar().<br>
      The same has to be done for the constraints. 
      For example, if you want to fill in the rows of a general MIP, 
      you have to call SCIPcreateConsLinear(), SCIPaddConsLinear() and additionally SCIPreleaseCons() 
      after finishing.
      If all variables and constraints are present, you can initiate the solution process via SCIPsolve().<br>
      Make sure to also call SCIPreleaseVar() if you do not need the variable pointer anymore.
      For an explanation of creating and releasing objects, please see the 
      <a href="http://scip.zib.de/doc/html/OBJ.html">doxygen documentation.</a>.
    </p>
    <p class="q">
      <a name="Q13">I already know a solution in advance, which I want to pass to SCIP. 
	How do I do this?</a>
    </p>
    <p class="a">
      First you have to build your problem, then you have to transform your problem 
      (SCIP only accepts solutions if it is at least in the transformed stage, 
      see <a href="#Q15">here</a>) via calling SCIPtransformProb(). 
      Next, you create a new SCIP primal solution by calling SCIPcreateSol() and set all nonzero values 
      by calling SCIPsetSolVal(). <br>
      After that, you add this solution by calling SCIPtrySol() 
      (success should be true afterwards, if your solution was correct) 
      and then release it by calling SCIPsolFree().
    </p>
    <p class="q">
      <a name="Q14">What operational stages of SCIP are there and are they important for me?</a>
    </p>
    <p class="a">
      There are ten different stages during a run of SCIP.
      <!--, a flow chart can be seen here: ????????. -->
      There are some methods which cannot be called in all stages, consider for example SCIPtrySol()
      (see <a href="#Q13">previous question</a>).
    </p>
    <p class="q">
      <a name="Q15">What is this thing with the original and the transformed problem about?</a>
    </p>
    <p class="a">
      Before the solving process starts, the original problem is copied to a different memory area.
      This copy is called &quot;transformed problem&quot;, and all modifications during the presolving
      and solving process are only applied to the transformed problem.<br>
      This has two main advantages: first, the user can also modify the problem after partially solving it.
      All modifications done by SCIP (presolving, cuts, variable fixings) during the partial solving process will
      be deleted together with the transformed problem, the user can modify the original problem and restart solving.<br>
      Second, the feasibility of solutions is always tested on the original problem!
    </p>
    <p class="q">
      <a name="Q16">Why do the names, e.g. in debug messages<!--DBG FAQ????--> often differ from the ones I
	defined?</a>
    </p>
    <p class="a">
      This can have several reasons. Especially names of binary variables can get different 
      prefixes and suffixes.
      Each transformed variable and constraint (see <a href="#Q15">here</a>) gets a &quot;t_&quot;
      as prefix.
      Apart from that, the meaning of original and transformed variables and constraints
      is identical.<br>
      General integers with bounds that differ just by one will be aggregated to binary variables which 
      get the same name with the suffix &quot;_bin&quot; .
      E.g. an integer variable t_x with lower bound 4 and upper bound 5 will be aggregated to a binary 
      variable t_x_bin = t_x - 4.<br>
      Variables can have negated counterparts, e.g. for a binary t_x its (also binary) negated would be 
      t_x_neg = 1 - t_x.<br>
      The knapsack constraint handler is able to disaggregate its constraints to cliques, which
      are set packing constraints and create names that consist of the knapsack's name and 
      a suffix &quot;_clq_&lt;int&gt;&quot;.
      E.g., a knapsack constraint 
      knap:&nbsp;x_1&nbsp;+&nbsp;x2&nbsp;+2&middot;x_3&nbsp;&le;&nbsp;2
      could be disaggregated to the set packing constraints 
      knap_clq_1:&nbsp;x_1&nbsp;+&nbsp;x_3&nbsp;&le;&nbsp;1 and 
      knap_clq_2:&nbsp;x_2&nbsp;+&nbsp;x_3&nbsp;&le;&nbsp;1.
    </p>
    <p class="q">
      <a name="Q17">What is SCIP_CALL()? Do I need this?</a>
    </p>
    <p class="a">
      Yes, you do. SCIP_CALL() is a global define, which handles the return codes of all methods 
      which return a SCIP_RETCODE and should therefore parenthesize each such method.
      SCIP_OKAY is the code which is returned if everything worked well; 
      there are 18 different error codes, see type_retcode.h.
      Each method that calls methods which return a SCIP_RETCODE should itself return a SCIP_RETCODE.
      If this is not possible, use SCIP_CALL_ABORT() to catch the return codes of the methods.
      If you do not want to use this either, you have to do the exception handling
      (i.e. the case that the return code is not SCIP_OKAY) on your own. 
    </p>		   
<h4><a name="Section4">Using SCIP as a Branch-Cut-And-Price-Framework</a></h4>
    <p class="q">
      <a name="Q19">What types of plugins can I add and how do I do this?</a>
    </p>
    <p class="a">
      See the <a href="http://scip.zib.de/doc/html/OBJ.html">doxygen documentation</a> for a list of 
      plugin types. There is a HowTo for each of them.
    </p>
    <p class="q">
      <a name="Q20">When should I implement a constraint handler, when should I implement a separator?</a>
    </p>
     <p  class="a">
     This depends on whether you want to add constraints or only cutting planes. The main
     difference is that constraints can be &quot;model constraints&quot;, while cutting planes are only
     additional LP rows that strengthen the LP relaxation.
     A model constraint is a constraint that is important for the feasibility of the integral
     solutions. If you deleted a model constraint, some infeasible integral vectors would
     suddenly become feasible in the reduced model.
     A cutting plane is redundant w.r.t. integral solutions. The set of feasible integral
     vectors does not change if the cutting plane is removed. You can, however, relax this
     condition slightly and add cutting planes that do cut off feasible solutions, as long as
     at least one of the optimal solutions remains feasible.
<br><br>
      You want to use a constraint handler in the following cases:
     </p>
     <ol> 
      <li>Some of your feasibility conditions can not be expressed by existing constraint types
      (e.g., linear constraints), or you would need too many of them. For example, the
      &quot;nosubtour&quot; constraint in the TSP is equivalent to exponentially many linear constraints.
      Therefore, it is better to implement a &quot;nosubtour&quot; constraint handler that can inspect
      solutions for subtours and generate subtour elimination cuts and others (e.g., comb
      inequalities) to strengthen the LP relaxation.</li>
      <li>Although you can express your feasibility condition by a reasonable number of existing
      constraint types, you can represent and process the condition in a more efficient way. For
      example, it may be that you can, due to your structural knowledge, implement a stronger or
      faster domain propagation or find tighter cutting planes than what one could do with the
      sum of the individual &quot;simple&quot; constraints that model the feasibility condition.</li>
    </ol>
    <p class="a">
      You want to use a cutting plane separator in the following cases:
    </p>
    <ol>
      <li>You have a general purpose cutting plane procedure that can be applied to any MIP. It
      does not use problem specific knowledge. It only looks at the LP, the integrality
      conditions, and other deduced information like the implication graph.</li>
      <li>You can describe your feasibility condition by a set C of constraints of existing type
      (e.g., linear constraints). The cuts you want to separate are model specific, but apart
      from these cuts, there is nothing you can gain by substituting the set C of constraints
      with a special purpose constraint. For example, the preprocessing and the domain
      propagation methods for the special purpose constraint would do basically the same as what
      the existing constraint handler does with the set C of constraints. In this case, you
      don't need to implement the more complex constraint handler. You add constraints of
      existing type to your problem instance in order to produce a valid model, and you enrich
      the model by your problem specific cutting plane separator to make the solving process
      faster. You can easily evaluate the performance impact of your cutting planes by enabling
      and disabling the separator.
      </li>
    </ol>
    <p class="a">
      Note that a constraint handler is defined by the type of constraints that it manages. For
      constraint handlers, always think in terms of constraint programming. For example, the
      &quot;nosubtour&quot; constraint handler in the TSP manages &quot;nosubtour&quot; constraints, which demand
      that in a given graph no feasible solution can contain a tour that does not contain all
      cities. In the usual TSP problem, there is only one &quot;nosubtour&quot; constraint, because there
      is only one graph for which subtours have to be ruled out.
      The &quot;nosubtour&quot; constraint handler has various ways of enforcing the &quot;nosubtour&quot; property
      of the solutions. A simple way is to just check each integral solution candidate (in the
      CHECK, ENFOLP, and ENFOPS callback methods) for subtours. If there is a subtour, the
      solution is rejected. A more elaborate way includes the generation of &quot;subtour elimination
      cuts&quot; in the SEPALP callback method of the constraint handler. Additionally, the
      constraint handler may want to separate other types of cutting planes like comb
      inequalities in its SEPALP callback.
      </p>  
    <p class="q">
      <a name="Q21">Can I remove unnecessary display columns or&mdash;even better&mdash;add my own ones? 
	Can I change the statistics displayed at the end of solving?</a>
    </p>
    <p class="a">
      Setting the status of a display column to 0 turns it off. 
      E.g., type "set display memused status 0" in the interactive shell to disable the 
      memory information column, or include the line 
      SCIPsetIntParam(scip, "display/memused/status", 0) into your source code.
      Adding an own display column can be realized via the SCIPincludeDisp() method, see the 
      <a href="http://scip.zib.de/doc/html/DISP.html">doxygen documentation</a>.<br>
      The statistic display, which is shown by "display statistics" and SCIPprintStatistics(), 
      respectively, cannot be changed.
    </p>
    <p class="q">
      <a name="Q22">How do LP-rows look like in SCIP?</a>
    </p>
    <p class="a">
      Each row is of the form <i>lhs</i> &le; &Sigma;(<i>val</i>[<i>j</i>]&middot;<i>col</i>[<i>j</i>])
      + <i>const</i> &le; <i>rhs</i>.
      For now, <i>val</i>[<i>j</i>]&middot;<i>col</i>[<i>j</i>] can be interpreted as 
      <i>a</i><sub><i>ij</i></sub>&middot;<i>x</i><sub><i>j</i></sub> 
      (for the difference between columns and variables see <a href="#Q24">here</a>).
      The constant is essentially needed for collecting the influence of presolving reductions like 
      variable fixings and aggregations.<br>
      The <i>lhs</i> and <i>rhs</i> may take infinite values: a less-than inequality would
      have <i>lhs</i> = -&#8734;, and a greater-than inequality would have <i>rhs</i> = +&#8734;.
      For equations <i>lhs</i> is equal to <i>rhs</i>.
      An infinite left hand side can be recognized by SCIPisInfinity(scip, -lhs),
      an infinite right hand side can be recognized by SCIPisInfinity(scip, rhs).
    </p>
    <p class="q">
      <a name="Q23">How do I get the data of the current LP-relaxation?</a>
    </p>
    <p class="a">
      You can get all rows in the current LP-relaxation by calling SCIPgetRowsData(). 
      The methods SCIProwGetConstant(), SCIProwGetLhs(), SCIProwGetRhs(), SCIProwGetVals(),
      SCIProwGetNNonz(),SCIProwGetCols() then give you information about each row, 
      see <a href="#Q22">previous question</a>.<br>
      You get a columnwise representation by calling SCIPgetColsData(). 
      The methods SCIPcolGetLb() and SCIPcolGetUb() give you the locally valid bounds of the LP relaxation
      in the current branch-and-bound-node.<br>
      If you are interested in global information, you have to call SCIPcolGetVar() to get the variable 
      associated to a column (see <a href="#Q24">next question</a>), which you can ask for global bounds
      via SCIPvarGetLbGlobal() and SCIPvarGetUbGlobal() as well as the type of the variable 
      (binary, general integer, implicit integer, or continuous) by calling SCIPvarGetType().
      For more information, also see <a href="#Q25">this question</a>.
    </p>
    <p class="q">
      <a name="Q24">What is the difference between columns and variables, rows and constraints?</a>
    </p>
    <p class="a">
      The terms columns and rows always refer to the representation in the current LP-relaxation, 
      variables and constraints to your global Constraint Integer Program.<br>
      Each column has an associated variable, which it represents, 
      but not every variable must be part of the current LP-relaxation. 
      E.g., it could be already fixed, aggregated to another variable, 
      or be priced out if a column generation approach was implemented.<br>
      Each row has either been added to the LP by a constraint handler or by a cutting plane separator.
      A constraint handler is able to, but does not need to, 
      add one or more rows to the LP as a linear relaxation of each of its constraints.
      E.g., in the usual case (i.e. without using dynamic rows) the linear constraint handler 
      adds one row to the LP for each linear constraint.
    </p>
    <p class="q">
      <a name="Q25">Are the variables and rows sorted in any particular order?</a>
    </p>
    <p class="a">The variable array which you get by SCIPgetVars() is internally sorted by variable types.
      The binary variables are stored at position [0,...,nbinvars-1], 
      the general integers at [nbinvars,...,nbinvars+nintvars-1]. 
      It holds that nvars = nbinvars + ninitvars + nimplvars + ncontvars.
      There is no further sorting within these sections, as well as there is no sorting for the rows.
      But each column and each row has a unique index, 
      which can be obtained by SCIPcolGetIndex() and SCIProwGetIndex(), respectively.
    </p>
    <p class="q">
      <a name="Q26">When should I use which of the numerical comparison functions?</a>
    </p>
    <p class="a">
      There are various numerical comparison functions available, each of them using a different
      epsilon in its comparisons. Let's take the equality comparison as an example. There are
      the following methods available: SCIPisEQ(), SCIPisSumEQ(), SCIPisFeasEQ(), SCIPisRelEQ(),
      SCIPisSumRelEQ().
<br><br>
      SCIPisEQ() should be used to compare two single values that are either results of a simple
      calculation or are input data. The comparison is done w.r.t. the &quot;numerics/epsilon&quot;
      parameter, which is 1e-9 in the default settings.<br>
      SCIPisSumEQ() should be used to compare the results of two scalar products or other &quot;long&quot;
      sums of values. In these sums, numerical inaccuracy can occur due to cancellation of
      digits in the addition of values with opposite sign. Therefore, SCIPisSumEQ() uses a
      relaxed equality tolerance of &quot;numerics/sumepsilon&quot;, which is 1e-6 in the default settings.<br>
      SCIPisFeasEQ() should be used to check the feasibility of some result, for example after
      you have calculated the activity of a constraint and compare it with the left and right
      hand sides. The feasibility is checked w.r.t. the &quot;numerics/feastol&quot; parameter, and 
      equality is defined in a relative fashion in contrast to absolute differences. That means,
      two values are considered to be equal if their difference divided by the larger of their
      absolute values is smaller than &quot;numerics/feastol&quot;. This parameter is 1e-6 in the default
      settings.<br>
      SCIPisRelEQ() can be used to check the relative difference between two values, just like
      what SCIPisFeasEQ() is doing. In contrast to SCIPisFeasEQ() it uses &quot;numerics/epsilon&quot; as
      tolerance.<br>
      SCIPisSumRelEQ() is the same as SCIPisRelEQ() but uses &quot;numerics/sumepsilon&quot; as tolerance.
      It should be used to compare two results of scalar products or other &quot;long&quot; sums.
    </p>
      <!-- FAQ Punkt: Verhaeltnis zu CPLEX callable library, concert, Coin... ????????? -->
</div>