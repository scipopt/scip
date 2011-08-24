<div id="faq">

    <h4><a name="Section0">Sections</a></h4>
    <ul>
     <li><a href="#Section1">General Questions about SCIP</a></li>
     <li><a href="#Section2">Using SCIP as a standalone MIP-Solver</a></li>
     <li><a href="#Section3">Using SCIP included in another source code</a></li>
     <li><a href="#Section4">Using SCIP as a Branch-Cut-And-Price-Framework</a></li>
     <li><a href="#Section5">Specific questions about Column Generation and Branch-And-Price with SCIP</a></li>
     <li><a href="#Section6">Specific questions about the copy functionality in SCIP</a></li>
    </ul>

    <h4><a name="Section1">General Questions about SCIP</a></h4>
    <p class="q">
      <a name="Q1">What is SCIP?</a>
    </p>

    <p class="a">
      SCIP is a constraint integer program solver. It can be used as a framework for <b>branch-cut-and-price</b> 
      and contains all necessary plugins to serve as a <b>standalone solver</b> for MIP, MIQCP, and PBO.<br> 
    <ul>
      <li>You can use the <a href="http://scip.zib.de/download.shtml">precompiled binaries</a>
          with which you can solve MIPs, MIQCPS, or PBO problems in *.mps, *.lp, *.rlp, *.zpl, *.opb, *.wbo and other file formats</li>
      <li>You can use SCIP as a subroutine for solving MIPs and more general constraint integer programs from your own source code.</li>
      <li>You can use SCIP as a framework in which you implement your own plugins.</li> 
      <li>You can use SCIP in any combination of the three purposes above.</li>
   </ul>
      This FAQ contains separate sections treating each of these usages of SCIP. It further considers specific questions for some features.
    <p class="q">
      <a name="Q2">When should I use SCIP?</a>
    </p>
    <p class="a">
      If you are either looking for a fast non-commercial MIP-solver or for a 
      branch-cut-and-price-framework in which you can directly implement your own methods 
      &mdash; which can even be for more general purposes than MIP (see, e.g., <a href="#Q18">here</a>).
    </p>
    <p class="q">
      <a name="Q3">I heard something about licenses. Do I have to pay for using SCIP?</a>
    </p>
    <p class="a"> 
      As long as you use it for academic, non-commercial purposes: No. 
      This will not change. For the other cases, check the explanation of the 
      <a href="http://zibopt.zib.de/ZIBLicence.shtml">ZIB academic license</a> and always feel free to ask us.
      If you want to use SCIP commercially, please write an e-mail to koch@zib.de .
    </p>
    <p class="q">
      <a name="Q4">Do I need any extra software?</a>
    </p>
    <p class="a">
      Except if you want to use SCIP as a pure CP-Solver (see <a href="#Q18">here</a>), 
      you need an underlying LP-Solver which has first to be installed and linked to the libraries 
      (see the <a href="INSTALL.html"><tt>INSTALL</tt></a> file in the SCIP root directory). 
      LP-solvers currently supported by SCIP are: 
    </p>
    <ul>
      <li><a href="http://www.coin-or.org/Clp/">CLP</a></li>
      <li><a href="http://www.fico.com/en/Products/DMTools/Pages/FICO-Xpress-Optimization-Suite.aspx">FICO XPress</a></li>
    <li><a href="http://www.gurobi.com/">Gurobi</a> (interface is in beta stage)</li>
      <li><a href="http://www-01.ibm.com/software/integration/optimization/cplex-optimizer/">ILOG CPLEX</a></li>
      <li><a href="http://www.mosek.com/">Mosek</a></li>
      <li><a href="http://www2.isye.gatech.edu/~wcook/qsopt/">QSopt</a></li>
      <li><a href="http://soplex.zib.de/">SoPlex</a></li>
    </ul>
    <p class="a">
      We also provide some <a href="http://scip.zib.de/download.shtml">precompiled binaries</a>. Besides that,
      you might need a modeling language like <a href="http://www.zib.de/koch/zimpl/">ZIMPL</a> to generate *.mps 
      or *.lp files. ZIMPL files can also directly be read by SCIP.
      You can download a package which includes SCIP, SoPlex and ZIMPL <a href="http://scip.zib.de/download.shtml">here</a>.
    </p>
    <p class="a">
    If you want to use SCIP for mixed integer nonlinear programming, you might want to use an underlying NLP solver (e.g., 
    <a href="https://projects.coin-or.org/Ipopt">Ipopt</a>) and an expression interpreter (e.g., <a href="http://www.coin-or.org/CppAD/">CppAD</a>).
    </p>
    <p class="q">
      <a name="Q5">How do I get started?</a>
    </p>
    <p class="a">
      An easy way is to use the SCIP-binaries and call SCIP from a shell, see <a href="SHELL.html">here</a> for a tutorial. 
      For that, you just have to download one of the precompiled binaries from the 
      <a href="http://scip.zib.de/download.shtml">download section</a>, or the zipped source code and compile it 
      with your favorite settings. This is described in detail in the <a href="INSTALL.html"><tt>INSTALL</tt></a> file in the SCIP main directory.
    </p>
    <p class="a">
      Another way is to use SCIP as a solver integrated into your own program source code.
      See the directories &quot;examples/MIPsolver/&quot and &quot;examples/Queens/&quot;
      for simple examples and <a href="#Q12">this point</a>.
    </p>
    <p class="a">
      A third way is to implement your own plugins into SCIP. 
      This is explained in the HowTos for all plugin types, which you can find in the 
      <a href="http://scip.zib.de/doc/html/index.html">doxygen documentation</a>.
      See also <a href="http://scip.zib.de/doc/html/START.html">How to start a new project</a>.
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
      In fact, there is a slight difference: <code>SCIPvarGetSol()</code> is also able to return pseudo solution values.
      If you do not have an idea, what pseudo solutions are, <code>SCIPgetVarSol()</code> should be just fine.<br> This
      should be the only case of 'duplicate methods'. If you find, however, another one, please contact us.
    </p>
     <p class="q">
       <a name="Q27">I have installation problems. What can I do?</a>
    </p>
    <p class="a">
      Read the <a href="INSTALL.html"><tt>INSTALL</tt></a> file in the SCIP root directory. It contains hints of how to get around problems. You can also try
      the binaries available on the <a href="http://scip.zib.de">SCIP page</a>.
    </p>
    <p class="a">
      If you want to use <a href="http://soplex.zib.de/">SoPlex</a> as the underlying LP-solver, you can try the
      following:
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
      If you do not have all of these libraries, read the <tt>INSTALL</tt> file in the ZIB Optimization
      Suite directory.<br>
      As a summary, the call of <code>make ZIMPL=false ZLIB=false READLINE=false</code> should work on most systems.
    </p>
    <p class="a">
      If you have any problems while using an LP-solver different from SoPlex,
      please read the SCIP <a href="INSTALL.html"><tt>INSTALL</tt></a> file first.                                                     
    </p>
    <p class="q">
      <a name="Q8">I changed to a new version of SCIP and now compiling breaks with some error messages which I don't understand.
	Do you have a general hint on that?</a>
    </p>
    <p class="a">
      Maybe the parameters of a function in SCIP changed. Relevant changes between version are listed below.
    </p>
    <ul>
        <li> <a href="http://scip.zib.de/doc/html/CHG1.html">Interface changes between SCIP 0.9 and SCIP 1.0</a></li>
        <li> <a href="http://scip.zib.de/doc/html/CHG2.html">Interface changes between SCIP 1.0 and SCIP 1.1</a></li>
        <li> <a href="http://scip.zib.de/doc/html/CHG3.html">Interface changes between SCIP 1.1 and SCIP 1.2</a></li>
        <li> <a href="http://scip.zib.de/doc/html/CHG4.html">Interface changes between SCIP 1.2 and SCIP 2.0</a></li>
    </ul>
    <p class="q">
      <a name="Q18">Can I use SCIP as a pure CP-Solver?</a>
    </p>
    <p class="a">
      Yes. There is a setting-file "settings/emphasis/cpsolver.set". 
      Furthermore, you can compile SCIP without any LP-Solver by <code>make LPS=none</code>.
    </p>    
    <p class="q">
       <a name="Q19">How can I use the SCIP makefiles?</a>
    </p>
    <p class="a">
      See <a href="http://scip.zib.de/doc/html/MAKE.html">Makefiles</a>.
    </p>
    <p class="q">
       <a name="Q19">What is this business with .a and .so libraries in the directory <tt>lib/</tt>?</a>
    </p>
    <p class="a">
       When SCIP builds the binary it needs to link with the corresponding libraries of the LP-solver. There are usually
       two ways to distribute a library. In the first (.a) the library is linked <em>statically</em> to SCIP; this means
       that all of its information is packed into the binary. In the second way (.so) the library is a <em>shared</em>
       library. In this case, the code of the library is not inserted into the binary itself, but is loaded at
       runtime. This has the advantage that the binaries are smaller, but it comes at the cost that you have to make
       sure that the library is found at runtime (for most systems it suffices to put the path of the library into the
       <tt>LD_LIBRARY_PATH</tt> variable).<br/><br/> <em>Note:</em> Depending on the system shared or static libraries
       are standard. If both links are present in the <tt>lib</tt> directory, the linker chooses which version to take
       (on newer Linux systems this is usually the shared version). If you do not want this to happen you have to delete
       the link that is not intended.
    </p>
    <p class="q">
       <a name="Q20">Can I compile SCIP as a shared library?</a>
    </p>
    <p class="a">
       You can use the <code>SHARED=true</code> option
       when making SCIP. This will generate the libraries of SCIP in shared format. The binary then
       also uses this form. Note that the path to the <code>lib</code> directory of SCIP is used to
       locate the libraries. If you want to move the libraries, you might to set the
       <code>LD_LIBRARY_PATH</code> environment variable to include the new path.

    <!--
       If you have a different system or compiler (which supports shared libraries) you might try to
       copy the file <code>make.linux.x86.gnu.opt-shared</code> and modify it according to your
       needs (see also <a href="http://scip.zib.de/doc/html/MAKE.html">Makefiles</a>).
     -->
       If you are using your own building system: The "magic" changes are the <code>-FPIC</code> 
       compiler/linker option and the <code>-Wl,-rpath</code> option.
    </p>
<h4><a name="Section2">Using SCIP as a standalone MIP-Solver</a></h4>
    <p class="q">
      <a name="Q28">The output is too wide for my terminal window. What can I do?</a>
    </p>
    <p class="a">
      In the interactive shell you can set the width of the output with the following command
      &quot;set display width&quot; followed by an appropriate number.
      See also the <a href="#Q9">next question</a>.		      
    </p>
    <p class="q">
      <a name="Q9">What do the cryptic abbreviations for the columns mean which are displayed 
	during the solving process of SCIP?</a>
    </p>
    <p class="a">
      Type "display display" in the interactive shell to get an explanation of them.<br> 
      By the way: If a letter appears in front of a display row, it indicates, 
      which heuristic found the new primal bound, a star representing an integral LP-relaxation.<br>
      Typing "display statistics" after finishing or interrupting the solving process gives 
      you plenty of extra information about the solving process.<br>
      (Typing "display heuristics" gives you a list of the heuristics including their letters.)
    </p>
    <p class="q">
      <a name="Q10">How do I change the behavior of SCIP?</a>
    </p>
    <p class="a">
      You can switch the settings for all presolving, heuristics, and separation plugins to three different modes
      via the "set {presolving, heuristics, separation} emphasis" parameters in the interactive shell.
      "off" turns off the respective type of plugins, "fast" chooses settings that lead to less time spent
      in this type of plugins, decreasing their impact, and "aggressive" increases the impact
      of this type of plugins.
      You can combine these general settings for cuts, presolving, and heuristics arbitrarily.<br>
      "display parameters" shows you which settings currently differ from their default, 
      "set default" resets them all. 
      Furthermore, there are complete settings that can be set by "set emphasis", i.e. settings for pure
      feasibility problems, solution counting, and CP like search.
    </p>
    <p class="q">
      <a name="Q11"> I recognized that one special plugin works very bad / very well for my problem and 
	I want to disable it / weaken its influence / intensify its influence. How do I do this?</a>
    </p>
    <p class="a">
    <ul>
      <li>For using a <b class="grey">non-default branching rule or node selection strategy</b> as
          standard, you just have to give it the highest priority, using
      <ul>
         <li><code> SCIP&gt; set branching &lt;name of a branching rule&gt; priority 9999999</code></li>
         <li><code> SCIP&gt; set nodeselectors &lt;name of a node selector&gt; priority 9999999</code></li>
      </ul> 
          With the commands
      <ul>
         <li><code> SCIP&gt; display branching</code></li>
         <li><code> SCIP&gt; display nodeselectors</code></li>
      </ul> 
          you get a list of all branching rules and node selectors, respectively.
          These lists give information about the different priorities.</li>
      <li>If you want to completely <b class="grey">disable a heuristic or a separator</b> you have
          to set its frequency to -1 and the sepafreq to -1 for separation by constraint handlers,
          respectively. The commands looks like this:
      <ul>
         <li><code> SCIP&gt; set heuristics &lt;name of a heuristic&gt; freq -1</code></li>
         <li><code> SCIP&gt; set separators &lt;name of a separator&gt; freq -1</code></li>
         <li><code> SCIP&gt; set constraints &lt;name of a constraint handler&gt; sepafreq -1</code></li>
      </ul></li> 
      <li>For <b class="grey">disabling a presolver</b>, you have to set its maxrounds parameter to 0.
      <ul>
         <li><code> SCIP&gt; set presolvers &lt;name of a presolver&gt; maxrounds 0</code></li>
      </ul></li> 
      <li>If you want to <b class="grey">intensify</b> the usage of a <b class="grey">heuristic</b>, 
          you can reduce its frequency to some smaller, positive value, and/or raise the quotient and 
          offset values (maxlpiterquot for diving heuristics, nodes for LNS heuristics).
      <ul>
         <li><code> SCIP&gt; set heuristics &lt;name of a heuristic&gt; freq &lt;some value&gt;</code></li>
         <li><code> SCIP&gt; set heuristic &lt;name of a diving heuristic&gt; maxlpiterquot &lt;some value&gt;</code></li>
         <li><code> SCIP&gt; set heuristic &lt;name of a LNS heuristic&gt; nodesquot &lt;some value&gt;</code></li>
      </ul></li>
      <li>For <b class="grey">intensifying</b> the usage of a <b class="grey">separator</b>, 
          you can raise its maxroundsroot and maxsepacutsroot values.
      <ul>
         <li><code> SCIP&gt; set separators &lt;name of a separator&gt; maxroundsroot &lt;some value&gt;</code></li>
         <li><code> SCIP&gt; set separators &lt;name of a separator&gt; maxrounds &lt;some value&gt;</code></li>
      </ul>
         If you also want to use this separator locally, you have to set its frequency to a positive value 
         and possibly raise maxrounds and maxsepacuts.
      <ul>
         <li><code> SCIP&gt; set separators &lt;name of a separator&gt; freq &lt;some value&gt;</code></li>
         <li><code> SCIP&gt; set separators &lt;name of a separator&gt; maxsepacuts &lt;some value&gt;</code></li>
      </ul>
         Compare the parameters of the heuristic/separator in the appropriate aggressive setting 
         (see <a href="#Q10">previous question</a>).</li>
      <li>For <b class="grey">weakening</b>, you should just do the opposite operation, i.e.,
          reducing the values  you would raise for intensification and vice versa.</li>
    </ul>
    </p>
<h4><a name="Section3">Using SCIP included in another source code</a></h4>
    <p class="q">
      <a name="Q12">How do I construct a problem instance in SCIP?</a>
    </p>
    <p class="a">
      First you have to create a SCIP object via <code>SCIPcreate()</code>, then you start to build the problem via
      <code>SCIPcreateProb()</code>.  Then you create variables via <code>SCIPcreateVar()</code> and add them to the
      problem via <code>SCIPaddVar()</code>.<br>

      The same has to be done for the constraints. For example, if you want to fill in the rows of a general MIP, you
      have to call <code>SCIPcreateConsLinear()</code>, <code>SCIPaddConsLinear()</code> and additionally
      <code>SCIPreleaseCons()</code> after finishing.  If all variables and constraints are present, you can initiate
      the solution process via <code>SCIPsolve()</code>.<br>

      Make sure to also call <code>SCIPreleaseVar()</code> if you do not need the variable pointer anymore.  For an
      explanation of creating and releasing objects, please see the <a
      href="http://scip.zib.de/doc/html/OBJ.html">doxygen documentation.</a>.
    </p>
    <p class="q">
      <a name="Q13">I already know a solution in advance, which I want to pass to SCIP. 
	How do I do this?</a>
    </p>
    <p class="a">
      First you have to build your problem, then you have to transform your problem (SCIP only accepts solutions if it
      is at least in the transformed stage, see <a href="#Q15">here</a>) via calling
      <code>SCIPtransformProb()</code>. Next, you create a new SCIP primal solution by calling
      <code>SCIPcreateSol()</code> and set all nonzero values by calling <code>SCIPsetSolVal()</code>. <br> After that,
      you add this solution by calling <code>SCIPtrySol()</code> (the variable <code>success</code> should be true
      afterwards, if your solution was correct) and then release it by calling <code>SCIPsolFree()</code>.
    </p>
    <p class="q">
      <a name="Q14">What operational stages of SCIP are there and are they important for me?</a>
    </p>
    <p class="a">
      There are ten different stages during a run of SCIP.
      <!--, a flow chart can be seen here: ????????. -->
      There are some methods which cannot be called in all stages, consider for example <code>SCIPtrySol()</code>
      (see <a href="#Q13">previous question</a>).
    </p>
    <p class="q">
      <a name="Q15">What is this thing with the original and the transformed problem about?</a>
    </p>
    <p class="a">
      Before the solving process starts, the original problem is copied.
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
      This can have several reasons. Especially names of binary variables can get different prefixes and suffixes.  Each
      transformed variable and constraint (see <a href="#Q15">here</a>) gets a &quot;t_&quot; as prefix.  Apart from
      that, the meaning of original and transformed variables and constraints is identical.<br>

      General integers with bounds that differ just by 1 will be aggregated to binary variables which get the same name
      with the suffix &quot;_bin&quot; .  E.g. an integer variable <code>t_x</code> with lower bound 4 and upper bound 5
      will be aggregated to a binary variable <code>t_x_bin = t_x - 4</code>.<br>

      Variables can have negated counterparts, e.g. for a binary <code>t_x</code> its (also binary) negated would be
      <code>t_x_neg = 1 - t_x</code>.<br>

      The knapsack constraint handler is able to disaggregate its constraints to cliques, which are set packing
      constraints and create names that consist of the knapsack's name and a suffix
      &quot;<code>_clq_&lt;int&gt;</code>&quot;.  E.g., a knapsack constraint
      <code>knap:&nbsp;x_1&nbsp;+&nbsp;x2&nbsp;+2 x_3&nbsp;&le;&nbsp;2</code> could be disaggregated to the set packing
      constraints <code>knap_clq_1:&nbsp;x_1&nbsp;+&nbsp;x_3&nbsp;&le;&nbsp;1</code> and
      <code>knap_clq_2:&nbsp;x_2&nbsp;+&nbsp;x_3&nbsp;&le;&nbsp;1</code>.
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
    <p class="q">
      <a name="Q17a">I want to stop the solving process after a certain time. How can I do this?</a>
    </p>
    <p class="a">
      Limits are given by parameters in SCIP, for example <code>limits/time</code> for a time limit or
      <code>limits/nodes</code> for a node limit. 
      If you want to set a limit, you have to change these parameters.
      For example, for setting the time limit to one hour, you have to call 
      <code>SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", 3600) )</code>.
      For more examples, please have a look into heur_rens.c.
      
    </p>		   
<h4><a name="Section4">Using SCIP as a Branch-Cut-And-Price-Framework</a></h4>
    <p class="q">
      <a name="Q18">How do I start a project?</a>
    </p>
    <p class="a">
      See <a href="http://scip.zib.de/doc/html/START.html">How to start a new project</a>.
    </p>
    <p class="q">
      <a name="Q19">What types of plugins can I add and how do I do this?</a>
    </p>
    <p class="a">
      See the <a href="http://scip.zib.de/doc/html/index.html">doxygen documentation</a> for a list of 
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
     solutions. If you delete a model constraint, some infeasible integral vectors would
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
      &quot;nosubtour&quot; constraint handler in the TSP
      (see "ConshdlrSubtour.cpp" in the directory "scip/examples/TSP/src/")
      manages &quot;nosubtour&quot; constraints, which demand
      that in a given graph no feasible solution can contain a tour that does not contain all
      cities. In the usual TSP problem, there is only one &quot;nosubtour&quot; constraint, because there
      is only one graph for which subtours have to be ruled out.
      The &quot;nosubtour&quot; constraint handler has various ways of enforcing the &quot;nosubtour&quot; property
      of the solutions. A simple way is to just check each integral solution candidate (in the
      CONSCHECK, CONSENFOLP, and CONSENFOPS callback methods) for subtours. If there is a subtour, the
      solution is rejected. A more elaborate way includes the generation of &quot;subtour elimination
      cuts&quot; in the CONSSEPALP callback method of the constraint handler. Additionally, the
      constraint handler may want to separate other types of cutting planes like comb
      inequalities in its CONSSEPALP callback.
      </p>  
    <p class="q">
      <a name="Q21">Can I remove unnecessary display columns or&mdash;even better&mdash;add my own ones? 
	Can I change the statistics displayed at the end of solving?</a>
    </p>
    <p class="a">
      Setting the status of a display column to 0 turns it off. E.g., type "set display memused status 0" in the
      interactive shell to disable the memory information column, or include the line <code>SCIPsetIntParam(scip,
      &quot;display/memused/status&quot;, 0)</code> into your source code.  Adding an own display column can be realized
      via the <code>SCIPincludeDisp()</code> method, see the <a href="http://scip.zib.de/doc/html/DISP.html">doxygen
      documentation</a>.<br> The statistic display, which is shown by "display statistics" and
      <code>SCIPprintStatistics()</code>, respectively, cannot be changed.
    </p>
    <p class="q">
      <a name="Q22">How do LP-rows look like in SCIP?</a>
    </p>
    <p class="a">
      Each row is of the form <i>lhs</i> &le; &Sigma;(<i>val</i>[<i>j</i>]&middot;<i>col</i>[<i>j</i>]) + <i>const</i>
      &le; <i>rhs</i>.  For now, <i>val</i>[<i>j</i>]&middot;<i>col</i>[<i>j</i>] can be interpreted as
      <i>a</i><sub><i>ij</i></sub>&middot;<i>x</i><sub><i>j</i></sub> (for the difference between columns and variables
      see <a href="#Q24">here</a>).  The constant is essentially needed for collecting the influence of presolving
      reductions like variable fixings and aggregations.<br> The <i>lhs</i> and <i>rhs</i> may take infinite values: a
      less-than inequality would have <i>lhs</i> = -&#8734;, and a greater-than inequality would have <i>rhs</i> =
      +&#8734;.  For equations <i>lhs</i> is equal to <i>rhs</i>.  An infinite left hand side can be recognized by
      <code>SCIPisInfinity(scip, -lhs)</code>, an infinite right hand side can be recognized by
      <code>SCIPisInfinity(scip, rhs)</code>.
    </p>
    <p class="q">
      <a name="Q23">How do I get the data of the current LP-relaxation?</a>
    </p>
    <p class="a">
      You can get all rows in the current LP-relaxation by calling <code>SCIPgetLPRowsData()</code>. The methods
      <code>SCIProwGetConstant()</code>, <code>SCIProwGetLhs()</code>, <code>SCIProwGetRhs()</code>,
      <code>SCIProwGetVals()</code>, <code>SCIProwGetNNonz()</code>, <code>SCIProwGetCols()</code> then give you
      information about each row, see <a href="#Q22">previous question</a>.<br> You get a columnwise representation by
      calling <code>SCIPgetLPColsData()</code>. The methods <code>SCIPcolGetLb()</code> and <code>SCIPcolGetUb()</code>
      give you the locally valid bounds of a column in the LP relaxation of the current branch-and-bound-node.<br> If
      you are interested in global information, you have to call <code>SCIPcolGetVar()</code> to get the variable
      associated to a column (see <a href="#Q24">next question</a>), which you can ask for global bounds via
      <code>SCIPvarGetLbGlobal()</code> and <code>SCIPvarGetUbGlobal()</code> as well as the type of the variable
      (binary, general integer, implicit integer, or continuous) by calling <code>SCIPvarGetType()</code>.  For more
      information, also see <a href="#Q25">this question</a>.
    </p>
    <p class="q">
      <a name="Q24">What is the difference between columns and variables, rows and constraints?</a>
    </p>
    <p class="a">
      The terms columns and rows always refer to the representation in the current LP-relaxation, variables and
      constraints to your global Constraint Integer Program.<br> Each column has an associated variable, which it
      represents, but not every variable must be part of the current LP-relaxation. E.g., it could be already fixed,
      aggregated to another variable, or be priced out if a column generation approach was implemented.<br> Each row has
      either been added to the LP by a constraint handler or by a cutting plane separator.  A constraint handler is able
      to, but does not need to, add one or more rows to the LP as a linear relaxation of each of its constraints.  E.g.,
      in the usual case (i.e. without using dynamic rows) the linear constraint handler adds one row to the LP for each
      linear constraint.
    </p>
    <p class="q">
      <a name="Q25">Are the variables and rows sorted in any particular order?</a>
    </p>
    <p class="a">The variable array which you get by <code>SCIPgetVars()</code> is internally sorted by variable types.
    The ordering is binary, integer, implicit integer and continuous variables, i.e., the binary variables are stored at
    position [0,...,nbinvars-1], the general integers at [nbinvars,...,nbinvars+nintvars-1], and so on. It holds that
    nvars = nbinvars + ninitvars + nimplvars + ncontvars.  There is no further sorting within these sections, as well as
    there is no sorting for the rows.  But each column and each row has a unique index, which can be obtained by
    <code>SCIPcolGetIndex()</code> and <code>SCIProwGetIndex()</code>, respectively.
    </p>
    <p class="q">
      <a name="Q26">When should I use which of the numerical comparison functions?</a>
    </p>
    <p class="a">
      There are various numerical comparison functions available, each of them using a different
      epsilon in its comparisons. Let's take the equality comparison as an example. There are
      the following methods available: <code>SCIPisEQ(), SCIPisSumEQ(), SCIPisFeasEQ(), SCIPisRelEQ(),
      SCIPisSumRelEQ()</code>.
    </p>
    <ul>
    <li>
      <code>SCIPisEQ()</code> should be used to compare two single values that are either results of a simple
      calculation or are input data. The comparison is done w.r.t. the &quot;numerics/epsilon&quot; parameter, which is
      1e-9 in the default settings.
    </li>
    <li>
      <code>SCIPisSumEQ()</code> should be used to compare the results of two scalar products or other &quot;long&quot;
      sums of values. In these sums, numerical inaccuracy can occur due to cancellation of digits in the addition of
      values with opposite sign. Therefore, <code>SCIPisSumEQ()</code> uses a relaxed equality tolerance of
      &quot;numerics/sumepsilon&quot;, which is 1e-6 in the default settings.
    </li>
    <li>
      <code>SCIPisFeasEQ()</code> should be used to check the feasibility of some result, for example after you have
      calculated the activity of a constraint and compare it with the left and right hand sides. The feasibility is
      checked w.r.t. the &quot;numerics/feastol&quot; parameter, and equality is defined in a relative fashion in
      contrast to absolute differences. That means, two values are considered to be equal if their difference divided by
      the larger of their absolute values is smaller than &quot;numerics/feastol&quot;. This parameter is 1e-6 in the
      default settings.
    </li>
    <li>
      <code>SCIPisRelEQ()</code> can be used to check the relative difference between two values, just like what
      <code>SCIPisFeasEQ()</code> is doing. In contrast to <code>SCIPisFeasEQ()</code> it uses
      &quot;numerics/epsilon&quot; as tolerance.
    </li>
    <li>
      <code>SCIPisSumRelEQ()</code> is the same as <code>SCIPisRelEQ()</code> but uses &quot;numerics/sumepsilon&quot;
      as tolerance.  It should be used to compare two results of scalar products or other &quot;long&quot; sums.
    </li>
    </ul>
<h4><a name="Section5">Specific questions about Column Generation and Branch-And-Price with SCIP</a></h4>
    <p class="q">
      <a name="Q29">Why are not all variables in the LP?</a>
    </p>
    <p class="a">
      With <code>SCIPgetLPColsData()</code> you can obtain the columns of the current LP relaxation. It is correct that
      not all variables are necessarily part of the current LP relaxation. In particular, in branch-and-price the
      variables generated at one node in the tree are not necessarily included in the LP relaxation of a different node
      (e.g., if the other node is not a descendant of the first node). But even if you are still at the same node or at
      a descendant node, SCIP can remove columns from the LP, if they are 0 in the LP relaxation. This dynamic column
      deletion can be avoided by setting the "removable" flag to FALSE in the <code>SCIPcreateVar()</code> call.
    </p>
    <p class="q">
      <a name="Q30">I only implemented one pricer, why is there a second one, called variable pricer?</a>
    </p>
    <p class="a">
      As described in the <a href="#Q29">previous question</a>, it may happen, that some variables are not in the
      current LP relaxation. Nevertheless, these variables still exist, and SCIP can calculate their reduced costs and
      add them to the LP again, if necessary. This is the job of the variable pricer. It is called before all other
      pricers.
    </p>
    <p class="q">
      <a name="Q31">How can I delete variables?</a>
    </p>
    <p class="a">
      Currently, it is not possible to delete variables during solving. A variable can be removed from the problem by
      <code>SCIPdelVar()</code>, but it remains in the corresponding constraints. Instead of that, you should either fix
      the variable to zero by <code>SCIPfixVar()</code>, which fixes the variable globally or use
      <code>SCIPchgVarUbNode()</code> and <code>SCIPchgVarLBNode()</code>, which change the bounds only for the current
      subtree.
    </p>
    <p class="a">
      We are working on this problem and expect to provide functionalities for the deletion of variables in the next
      release.
    </p>
    <p class="q">
      <a name="Q32">How can I store branching decisions?</a>
    </p>
    <p class="a">
      This is a very common problem in Branch-And-Price, which you can deal nicely with using SCIP.
      There are basically three different options.
    </p>
    <p class="a">
      The first one is to add binary variables to the problem that encode branching decisions. Then constraints should
      be added that enforce the corresponding branching decisions in the subtrees.
    </p>
    <p class="a">
      If you have complex pricer data like a graph and need to update it after each branching decision,
      you should introduce "marker constraints" that are added to the branching nodes and store
      all the information needed (see the <a href="#Q34">next question</a>).
    </p>
    <p class="a">
      The third way is the usage of an event handler like described in <a href="#Q35">How can an event handler help me
      with my branching?</a>.
    </p>
    <p class="q">
      <a name="Q34">I want to store some information at the nodes and update my pricer's data structures
                    when entering a new node. How can I do that?</a>
    </p>
    <p class="a">
      This can be done by creating a new constraint handler with constraint data that can store the information and
      do/undo changes in the pricer's data structures.
    </p>
    <p class="a">
      In general, all methods of the constraint handler (check, enforcing, separation, ...) should be empty (which means
      to always return the status SCIP_FEASIBLE for the fundamental callbacks), just as if all constraints of this type
      are always feasible. The important callbacks are the CONSACTIVE and CONSDEACTIVE methods.
    </p>
    <p class="a">
      The CONSACTIVE method is always called when a node is entered on which the constraint has been added. Here, you
      need to apply the changes to your pricing data structures. The CONSDEACTIVE method will be called if the node is
      left again. Since the CONSACTIVE and CONSDEACTIVE methods of different constraints are always called in a
      stack-like fashion, this should be exactly what you need.
    </p>
    <p class="a">
      If you have such a constraint handler, just create constraints of this type and add them to the child nodes of
      your branching by <code>SCIPaddConsNode()</code>. Make sure to set the "stickingatnode" flag to TRUE in order to
      prevent SCIP from moving the constraint around in the tree.
    </p>
    <p class="a">
      If you need to fix variables for enforcing your branching decision, this can be done in the propagation callback
      of the constraint handler. Since in general each node is only propagated once, in this case you will have to check
      in your CONSACTIVE method whether new variables were added after your last propagation of this node. If this is
      the case, you will have to mark this node for repropagation by <code>SCIPrepropagateNode()</code>.
    </p>
    <p class="a">
      You can look into the constraint handler of the coloring problem (examples/Coloring/src/cons_storeGraph.c) to get
      an example of a constraint handler that does all these things.
    </p>
    <p class="q">
      <a name="Q35">How can an event handler help me with my branching?</a>
    </p>
    <p class="a">
      An event handler can watch for events like local bound changes on variables. So, if your 
      pricer wants to be informed whenever a local bound of a certain variable changes, add an 
      event handler, catch the corresponding events of the variable, and in the event handler's 
      execution method adjust the data structures of your pricer accordingly.
    </p>

    <p class="q">
      <a name="Q36">How can I add locally valid variables to the problem in my branch-and-price code?</a>
    </p>
    <p class="a">
      Variables in SCIP are always added globally. If you want to add them locally, because they are forbidden
      in another part of the branch-and-bound-tree, you should look to it that they are locally fixed to 0
      in all subtrees where they are not valid. A description of how this can be done is given <a href="#Q34">here</a>.
    </p>


    <p class="q">
      <a name="Q37">My pricer generates the same column twice. How can I solve this problem?</a>
    </p>
    <p class="a">
      First check whether your pricing is correct. Are there upper bounds on variables that you have forgotten to take
      into account? If your pricer cannot cope with other variable bounds than 0 and infinity, you have to mark 
      all constraints containing priced variables as modifiable, and you may have to disable reduced cost
      strengthening by setting propagating/rootredcost/freq to -1.
    </p>
    <p class="a">  
      If your pricer works correctly and makes sure that the same column is added at most once in one pricing round, this
      behavior is probably caused by the PRICER_DELAY property of your pricer.
    </p>
    <p class="a">
      If it is set to FALSE, the following may have happened: The variable pricer (see <a href="#Q30">this question</a>)
      found a variable with negative dual feasibility that was not part of the current LP relaxation and added it to the
      LP. In the same pricing round, your own pricer found the same column and created a new variable for it. This might
      happen, since your pricer uses the same dual values as the variable pricer. To avoid this behavior, set
      PRICER_DELAY to TRUE, so that the LP is reoptimized after the variable pricer added variables to the LP. You can
      find some more information about the PRICER_DELAY property at <a
      href="http://scip.zib.de/doc/html/PRICER.html#PRICER_PROPERTIES">How to add variable pricers </a>.
    </p>

    <p class="q">
      <a name="Q38">Which default plugins should be deactivated in order to get a working branch-and price code?</a>
    </p>
    <p class="a">
      In most cases, you should deactivate separators, since cutting planes that are added to your master problem may
      destroy your pricing problem. Additionally, it may be necessary to deactivate some presolvers, mainly the dual
      fixing presolver. This can be done by not including these plugins into SCIP, namely by not calling
      <code>SCIPincludeSepaXxx()</code> and <code>SCIPincludePresolXxx()</code> in your own plugins-including
      files. Alternatively, you can set the parameters maxrounds and maxroundsroot to zero for all separators and
      maxrounds to zero for the presolvers.
    </p>

    <p class="q">
      <a name="Q39">What are the lazy bounds for variables in SCIP and what do I need them for?</a>
    </p>
    <p class="a">
      In many Branch-and-Price applications, you have binary variables, but you do not want to impose
      upper bounds on these variables in the LP relaxation, because the upper bound is implicitly enforced
      by the problem contraints, anyway. If the upper bounds are explicitly added to the LP, they lead to 
      further dual variables, which may be hard to take into account in the pricing problem.
    </p>
    <p class="a">  
      There are two possibilities how to solve this problem. First, you could change the binary variables 
      to general integer variables, if this does not change the problem. However, if you use special linear 
      constraints like set partitioning/packing/covering, you can only add binary variables to these constraints. 
    </p>
    <p class="a">  
      In order to still allow the usage of these types of constraints in a branch-and-price approach, the 
      concept of lazy bounds was introduced in SCIP 2.0.
      For each variable, you can define lazy upper and lower bounds, i.e. bounds, that are implicitly enforced by
      constraints in the problem. SCIP adds variable bounds to the LP only if the bound is tighter than the 
      corresponding lazy bound. In order to use lazy bounds, you have to make sure, however, that each feasible
      LP solution respects these bounds (even non-optimal LP solutions).
    </p>
    <p class="a">  
      For instance, if you have set partioning constraints in your problem, you can define variables contained in these
      constraints as binary and set the lazy upper bound to 1, which allows you to use the better propagation methods
      of the setppc constraint handler compared to the linear constraint handler without taking care about upper 
      bounds on variables in the master.
    </p>

    <p class="q">
      <a name="Q40">Can I stop the pricing process before the master problem is solved to optimality?</a>
    </p>
    <p class="a">
      In a column generation approach, you usually have to solve the master problem to optimality since otherwise,
      its objective function value is no valid dual bound. However, there is a way in SCIP to
      stop the pricing process earlier, called "early branching".
    </p>
    <p class="a">  
      The reduced cost pricing method of a pricer has a result pointer that should be set each time the method is called.
      In the usual case that the pricer either adds a new variable or ensures that there are no further 
      variables with negative dual feasibility, the result pointer should be set to SCIP_SUCCESS. 
      If the pricer aborts pricing without creating a new variable, but there might exist additional 
      variables with negative dual feasibility, the result pointer should be set to SCIP_DIDNOTRUN. 
      In this case, the lp solution will not be used as a lower bound. 
      The pricer can, however, store a valid lower bound in the lowerbound pointer.
    </p>

<h4><a name="Section6">Specific questions about the copy functionality in SCIP</a></h4>
    <p class="q">
      <a name="Q41">What is <code>SCIPcopy()</code>?</a>
    </p>
    <p class="a">
      The functionality of copying a SCIP model is added to the SCIP version 2.0.0. It gives the possibility to 
      generate a copy of the current SCIP model. This functionality is of interest, for example, for large 
      neighborhood heuristics (such as heur_rens.c). They can now easily copy the complete problem and fix  a certain 
      set of variables and work therefore on a reasonable copy of the original problem.
    </p>
    <p class="q">
      <a name="Q42">How do I get a copy of a variable or a constraint</code>?</a>
    </p>
    <p class="a">
      For the variables and constraints there are the methods <code>SCIPgetVarCopy()</code> and 
      <code>SCIPgetConsCopy()</code> which provide a copy for a variable or a constraint, respectively.
    </p>
    <p class="q">
      <a name="Q39">What does the valid pointer in case of the constraint handler and variable pricer mean?</a>
    </p>
    <p class="a">
      SCIP would like to know if the copied problem is a one to one copy. That means, if all problem defining objects 
      got successfully copied. If this is the case, all reductions made in the copy can be transformed into 
      original instance. The problem defining objects in SCIP are the constraint handlers and the variable pricers.
    </p>
          <!-- FAQ Punkt: Verhaeltnis zu CPLEX callable library, concert, Coin... ????????? -->
          <!-- VBC tool, visualisierung ????????? -->
</div>

