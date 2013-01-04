<div id="faq">
    <ol>
     <li><a href="#Section1"><b>General Questions about SCIP</b></a></li>
      <ul>
	<li><a href="#Q1.1">What is SCIP?</a></li>
        <li><a href="#Q1.2">When should I use SCIP?</a></li>
        <li><a href="#Q1.3">I heard something about licenses. Do I have to pay for using SCIP?</a></li>
	<li><a href="#Q1.4">How do I get started?</a></li>
	<li><a href="#Q1.5">Do I need any extra software?</a></li>
	<li><a href="#Q1.6">I have installation problems. What can I do?</a></li>
	<li><a href="#Q1.7">I changed to a new version of SCIP and now compiling breaks with some error messages which
        I don't understand. Do you have a general hint on that?</a></li>
	<li><a href="#Q1.8">How can I debug in SCIP?</a></li>
	<li><a href="#Q1.9">SCIP prints error messages, aborts, produces segmentation faults, or just behaves strangely.
         What should I do?</a></li>
	<li><a href="#Q1.10"> I would like to check whether some functionality is implemented in SCIP.
        How does the naming of the methods work? Where do I find the most common methods?</a></li>
	<li><a href="#Q1.11">Can I use SCIP as a pure CP-Solver?</a></li>
	<li><a href="#Q1.12">Can I use SCIP as a pure LP-Solver?</a></li>
	<li><a href="#Q1.13">Which kind of MINLPs are supported by SCIP?</a></li>
	<li><a href="#Q1.14">How can I use the SCIP makefiles?</a></li>
	<li><a href="#Q1.15">What is this business with .a and .so libraries in the directory <tt>lib/</tt>?</a></li>
	<li><a href="#Q1.16">Can I compile SCIP as a shared library?</a></li>
	<li><a href="#Q1.17">The methods SCIPgetVarSol() and SCIPvarGetSol() seem to have the same
      functionality. Which one should I use?</a></li>
      </ul>
     <li><a href="#Section2"><b>Using SCIP as a standalone solver</b></a></li>
      <ul>
	<li><a href="#Q2.1">The output is too wide for my terminal window. What can I do?</a></li>
	<li><a href="#Q2.2">What do the cryptic abbreviations for the columns mean which are displayed
        during the solving process of SCIP?</a></li>
	<li><a href="#Q2.3">How do I change the behavior of SCIP?</a></li>
	<li><a href="#Q2.4"> I recognized that one special plugin works very poorly / very well for my problem and
        I want to disable it / weaken its influence / intensify its influence. How do I do this?</a></li>
	<li><a href="#Q2.5">How can I input a MINLP into SCIP?</a></li>
      </ul>
     <li><a href="#Section3"><b>Using SCIP included in another source code</b></a></li>
      <ul>
	<li><a href="#Q3.1">How do I construct a problem instance in SCIP?</a></li>
	<li><a href="#Q3.2">I already know a solution in advance, which I want to pass to SCIP.
        How do I do this?</a></li>
	<li><a href="#Q3.3">What operational stages of SCIP are there and are they important for me?</a></li>
	<li><a href="#Q3.4">What is the difference between the original and the transformed problem?</a></li>
	<li><a href="#Q3.5">Why do the names, e.g., in debug messages often differ from the ones I
        defined?</a></li>
	<li><a href="#Q3.6">What is SCIP_CALL()? Do I need this?</a></li>
	<li><a href="#Q3.7">I want to stop the solving process after a certain time. How can I do this?</a></li>
      </ul>
     <li><a href="#Section4"><b>Using SCIP as a Branch-Cut-And-Price-Framework</b></a></li>
      <ul>
	<li><a href="#Q4.1">How do I start a project?</a></li>
	<li><a href="#Q4.2">What types of plugins can I add and how do I do this?</a></li>
	<li><a href="#Q4.3">When should I implement a constraint handler, when should I implement a separator?</a></li>
	<li><a href="#Q4.4">Can I remove unnecessary display columns or&mdash;even better&mdash;add my own ones?
        Can I change the statistics displayed at the end of solving?</a></li>
	<li><a href="#Q4.5">What do LP-rows look like in SCIP?</a></li>
	<li><a href="#Q4.6">How do I get the data of the current LP-relaxation?</a></li>
	<li><a href="#Q4.7">What is the difference between columns and variables, rows and constraints?</a></li>
	<li><a href="#Q4.8">Are the variables and rows sorted in any particular order?</a></li>
	<li><a href="#Q4.9">When should I use which of the numerical comparison functions?</a></li>
	<li><a href="#Q4.10">How do I solve an LP inside my SCIP plugin?</a></li>
      </ul>
     <li><a href="#Section5"><b>Specific questions about Column Generation and Branch-And-Price with SCIP</b></a></li>
      <ul>
	<li><a href="#Q5.1">Why are not all variables in the LP?</a></li>
	<li><a href="#Q5.2">I only implemented one pricer, why is there a second one, called variable pricer?</a></li>
	<li><a href="#Q5.3">How can I store branching decisions?</a></li>
	<li><a href="#Q5.4">I want to store some information at the nodes and update my pricer's data structures
        when entering a new node. How can I do that?</a></li>
	<li><a href="#Q5.5">How can an event handler help me with my branching?</a></li>
	<li><a href="#Q5.6">How can I add locally valid variables to the problem in my branch-and-price code?</a></li>
	<li><a href="#Q5.7">My pricer generates the same column twice. How can I solve this problem?</a></li>
	<li><a href="#Q5.8">Which default plugins should be deactivated in order to get a working branch-and-price code?</a></li>
	<li><a href="#Q5.9">What are the lazy bounds for variables in SCIP and what do I need them for?</a></li>
	<li><a href="#Q5.10">Can I stop the pricing process before the master problem is solved to optimality?</a></li>
	<li><a href="#Q5.11">How can I delete variables?</a></li>
      </ul>
     <li><a href="#Section6"><b>Specific questions about the copy functionality in SCIP</b></a></li>
      <ul>
	<li><a href="#Q6.1">What is <code>SCIPcopy()</code>?</a></li>
	<li><a href="#Q6.2">How do I get a copy of a variable or a constraint?</a></li>
	<li><a href="#Q6.3">What does the <code>valid</code> pointer in the copy callback of the constraint handler and variable pricer mean?</a></li>
      </ul>
    </ol>


    <h4><a name="Section1">1. General Questions about SCIP</a></h4>
    <p class="q">
      <a name="Q1.1">What is SCIP?</a>
    </p>

    <p class="a">
      SCIP is a constraint integer program solver. It can be used as a framework for <b>branch-cut-and-price</b> 
      and contains all necessary plugins to serve as a <b>standalone solver</b> for MIP, MINLP, and PBO.<br>
    <ul>
      <li>You can use the <a href="http://scip.zib.de/download.shtml">precompiled binaries</a>
          to solve MIPs, MINLPs, or PBO problems. They can read
      <p class="a">
        <ol>
	  <li>CIP - file format (own constraint based file format)</li>
	  <li><a href="http://people.sc.fsu.edu/~jburkardt/data/cnf/cnf.html">CNF - file format</a></li>
	  <li><a href="http://www.gams.com/docs/document.htm">GMS format [write only]</a></li>
	  <li><a href="http://www.g12.cs.mu.oz.au/minizinc/downloads/doc-1.5/flatzinc-spec.pdf">FZN - file format</a></li>
          <li><a href="http://publib.boulder.ibm.com/infocenter/cosinfoc/v12r3/topic/ilog.odms.cplex.help/Content/Optimization/Documentation/Optimization_Studio/_pubskel/ps_reffileformatscplex2159.html">(R)LP - file formats</a></li>
          <li><a href="http://en.wikipedia.org/wiki/MPS_%28format%29">MPS - file format</a></li>
          <li><a href="http://www.cril.univ-artois.fr/PB12/format.pdf">OPB/WBO - file formats</a></li>
	  <li><a href="http://www.coin-or.org/OS/OSiL.html">OSiL format for MINLP</a></li>
          <li><a href="http://polip.zib.de/pipformat.php">PIP - file format</a></li>
	  <li>SOL - file format (solution reader, xml-solutions are also supported)</li>
          <li><a href="http://zimpl.zib.de/download/zimpl.pdf">ZPL - file format</a></li>
        </ol>
      </p>
      <p class="a">
	 See also the page on <a href="group__FILEREADERS.shtml">file readers</a>.
      </p>
      </li>
      <li>You can use SCIP as a subroutine for solving MIPs and more general constraint integer programs from your own source code.</li>
      <li>You can use SCIP as a framework in which you implement your own plugins.</li> 
      <li>You can use SCIP in any combination of the three purposes above.</li>
   </ul>
   <p class="a">
      This FAQ contains separate sections covering each of these usages of SCIP. It further considers specific questions for some features.
   </p>
    <p class="q">
      <a name="Q1.2">When should I use SCIP?</a>
    </p>
    <p class="a">
      If you are either looking for a fast non-commercial MIP-solver or for a 
      branch-cut-and-price-framework in which you can directly implement your own methods 
      &mdash; which can even be for more general purposes than MIP (see, e.g., <a href="#Q1.11">here</a>).
    </p>
    <p class="q">
      <a name="Q1.3">I heard something about licenses. Do I have to pay for using SCIP?</a>
    </p>
    <p class="a"> 
      As long as you use it for academic, non-commercial purposes: No. 
      This will not change. For the other cases, check the explanation of the 
      <a href="http://scip.zib.de/licence.shtml">ZIB academic license</a> and always feel free to ask us.
      If you want to use SCIP commercially, please write an e-mail to koch@zib.de .
    </p>
    <p class="q">
      <a name="Q1.4">How do I get started?</a>
    </p>
    <p class="a">
      An easy way is to use the SCIP-binaries and call SCIP from a shell, see <a href="SHELL.shtml">here</a> for a tutorial.
      For that, you just have to download one of the precompiled binaries from the 
      <a href="http://scip.zib.de/download.shtml">download section</a>, or the zipped source code and compile it 
      with your favorite settings. This is described in detail in the <a href="INSTALL.shtml"><tt>INSTALL</tt></a> file in the SCIP main directory.
    </p>
    <p class="a">
      Another way is to use SCIP as a solver integrated into your own program source code.
      See the directories &quot;examples/MIPsolver/&quot and &quot;examples/Queens/&quot;
      for simple examples and <a href="#Q3.1">this point</a>.
    </p>
    <p class="a">
      A third way is to implement your own plugins into SCIP. 
      This is explained in the HowTos for all plugin types, which you can find in the 
      <a href="index.shtml">doxygen documentation</a>.
      See also <a href="START.shtml">How to start a new project</a>.
    </p>
    <p class="q">
      <a name="Q1.5">Do I need any extra software?</a>
    </p>
    <p class="a">
      Unless you want to use SCIP as a pure CP-Solver (see <a href="#Q1.11">here</a>),
      you need an underlying LP-Solver installed and linked to the libraries
      (see the <a href="INSTALL.shtml"><tt>INSTALL</tt></a> file in the SCIP root directory).
      LP-solvers currently supported by SCIP are: 
    </p>
    <ul>
      <li><a href="https://projects.coin-or.org/Clp/">CLP</a></li>
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
       <a name="Q1.6">I have installation problems. What can I do?</a>
    </p>
    <p class="a">
      Read the <a href="INSTALL.shtml"><tt>INSTALL</tt></a> file in the SCIP root directory. It contains hints of how to get around problems. You can also try
      the binaries available on the <a href="http://scip.zib.de">SCIP page</a>.
    </p>
    <p class="a">
      If you want to use <a href="http://soplex.zib.de/">SoPlex</a> as the underlying LP-solver, you can try the
      following:
      First, download the <a href="http://scip.zib.de">SCIP Optimization Suite</a>.
      Then, extract the file, change into the scipoptsuite directory, and enter 'make'.
      As long as you have all the necessary libraries installed in your system, it should generate a SCIP
      binary linked to ZIMPL and Soplex. The necessary system libraries are:
    </p>
    <ol>
        <li>ZLIB (libz.a)</li>
        <li>GMP (libgmp.a)</li>
        <li>Readline (libreadline.a)</li>
    </ol>
    <p class="a">
      If you do not have all of these libraries, read the <tt>INSTALL</tt> file in the SCIP Optimization
      Suite directory.<br>
      As a summary, the call of <code>make ZIMPL=false ZLIB=false READLINE=false</code> should work on most systems.
    </p>
    <p class="a">
      If you have any problems while using an LP-solver different from SoPlex,
      please read the SCIP <a href="INSTALL.shtml"><tt>INSTALL</tt></a> file first.
    </p>
    <p class="q">
      <a name="Q1.7">I changed to a new version of SCIP and now compiling breaks with some error messages which I don't understand.
      Do you have a general hint on that?</a>
    </p>
    <p class="a">
      Maybe the parameters of a function in SCIP changed. Relevant changes between version are listed below.
    </p>
    <ul>
        <li> <a href="CHG1.shtml">Interface changes between SCIP 0.9 and SCIP 1.0</a></li>
        <li> <a href="CHG2.shtml">Interface changes between SCIP 1.0 and SCIP 1.1</a></li>
        <li> <a href="CHG3.shtml">Interface changes between SCIP 1.1 and SCIP 1.2</a></li>
        <li> <a href="CHG4.shtml">Interface changes between SCIP 1.2 and SCIP 2.0</a></li>
        <li> <a href="CHG5.shtml">Interface changes between SCIP 2.0 and SCIP 2.1</a></li>
        <li> <a href="CHG6.shtml">Interface changes between SCIP 2.1 and SCIP 3.0</a></li>
    </ul>
    <p class="q">
      <a name="Q1.8">How can I debug in SCIP?</a>
    </p>
    <p class="a">
      Compile SCIP in debug mode: <code>make OPT=dbg</code>. Put the binary into a debugger, e.g., <a href="http://www.gnu.org/s/gdb/">gdb</a> 
      and let it run again. If you get an impression which component is causing the trouble, set <code>#define SCIP_DEBUG</code> as the first 
      line of the corresponding *.c file, recompile, let it run again. This will print debug messages from that piece of code. 
      Find a short debugging tutorial <a href="DEBUG.shtml">here</a>.
    </p>
    <p class="q">
      <a name="Q1.9">SCIP prints error messages, aborts, produces segmentation faults, or just behaves strangely. What should I do?</a>
    </p>
    <p class="a">
       <a href="#Q1.8">See above.</a> Often, the asserts that show up in debug mode already help to clarify misunderstandings
       and suggest fixes, if you were calling SCIP functions in an unexpected manner. If you, however, have the impression
       that there is a bug in SCIP, <a href="https://www.edom.mi.uni-erlangen.de/cgi/scipweb/submit.pl">send us a bug report</a>
       and include a log file and ideally a backtrace from a debugger.
    </p>
    <p class="q">
      <a name="Q1.10"> I would like to check whether some functionality is implemented in SCIP.
      How does the naming of the methods work? Where do I find the most common methods?</a>
    </p>
    <p class="a">
      For an explanation of the naming see the 
      <a href="CODE.shtml">coding style guidelines</a>. <br>
      The methods you might want to use should normally be found in scip.h or in pub_*.h,
      see also <a href="DOC.shtml">here</a>.
      In the doxygen documentation of <a href="scip_8h.shtml">scip.h</a>,
      the methods are ordered by topics.
    </p>
    <p class="q">
      <a name="Q1.11">Can I use SCIP as a pure CP-Solver?</a>
    </p>
    <p class="a">
      Yes. SCIP can be used as a pure CP-Solver by typing
      "set emphasis cpsolver" in the shell or by using the function <code>SCIPsetEmphasis()</code>.
      Furthermore, you can compile SCIP without any LP-Solver by <code>make LPS=none</code>.
      See <a href="#Q2.3">here</a> for more on changing the behavior of SCIP.
    </p>    
    <p class="q">
      <a name="Q1.12">Can I use SCIP as a pure LP-Solver?</a>
    </p>
    <p class="a"> 
      Since LPs are only special types of MIPs and CIPs, the principal answer is yes. If you feed a pure LP to SCIP, it
      will first apply presolving and then hand this presolved problem to the underlying LP solver. If solved to
      optimality, you can query the optimal solution values as always.<br/><br/>

      However, there are certain limitations to this: Dual multipliers and reduced costs are not accessible. If the LP
      turns out to be infeasible, you cannot obtain the Farkas proof. Strong LP duality theory does
      not apply to MIP or CIP and there is no special functionality implemented in SCIP for the case
      that the solved problem is an LP.

      Hence, if you need more, "LP specific", information than the primal solution, you are better off using an LP-Solver
      directly. If you are using the SCIP Optimization Suite, you could, e.g., use the included LP solver
      <a href=http://soplex.zib.de>SoPlex</a>.  If you want to solve an LP not from the command line, but within your
      C/C++ program, you could also use SCIP's LP-Interface, see also <a href="#Q4.10">here</a>.
    </p>
    <p class="q">
      <a name="Q1.13">Which kind of MINLPs are supported by SCIP?</a>
    </p>
    <p class="a">
      SCIP supports nonlinear constraints of the form lhs &le; f(x) &le; rhs, where the function f(x) is an algebraic
      expression that can be represented as expression tree. Such an expression tree has constants and variables as terminal
      nodes and operands as non-terminal nodes.
      Expression operands supported by SCIP include addition, subtraction, multiplication, division, exponentiation and
      logarithm. Trigonometric functions are not yet supported by SCIP.<br>
      Nonlinear objective functions are not supported by SCIP and must be modeled as constraint function.<br><br>
      <em>Note</em>, that the support for non-quadratic nonlinear constraints is still in a <b>BETA</b>-stadium and not yet as
      robust as the rest of SCIP.
      Missing bounds on nonlinear variables and tiny or huge coefficients can easily lead to numerical problems,
      which can be avoided by careful modeling.
    </p>
    <p class="q">
       <a name="Q1.14">How can I use the SCIP makefiles?</a>
    </p>
    <p class="a">
      See <a href="MAKE.shtml">Makefiles</a>.
    </p>
    <p class="q">
       <a name="Q1.15">What is this business with .a and .so libraries in the directory <tt>lib/</tt>?</a>
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
       <a name="Q1.16">Can I compile SCIP as a shared library?</a>
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
       needs (see also <a href="MAKE.shtml">Makefiles</a>).
     -->
       If you are using your own building system: The "magic" changes are the <code>-FPIC</code> 
       compiler/linker option and the <code>-Wl,-rpath</code> option.
    </p>
    <p class="q">
      <a name="Q1.17">The methods SCIPgetVarSol() and SCIPvarGetSol() seem to have the same
      functionality. Which one should I use?</a>
    </p>
    <p class="a">
      In fact, there is a slight difference: <code>SCIPvarGetSol()</code> is also able to return pseudo solution values.
      If you do not have an idea, what pseudo solutions are, <code>SCIPgetVarSol()</code> should be just fine.<br> This
      should be the only case of 'duplicate methods'. If you find, however, another one, please contact us.
    </p>


<h4><a name="Section2">2. Using SCIP as a standalone solver</a></h4>
    <p class="q">
      <a name="Q2.1">The output is too wide for my terminal window. What can I do?</a>
    </p>
    <p class="a">
      In the interactive shell you can set the width of the output with the following command
      &quot;set display width&quot; followed by an appropriate number.
      See also the <a href="#Q2.2">next question</a>.
    </p>
    <p class="q">
      <a name="Q2.2">What do the cryptic abbreviations for the columns mean which are displayed
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
      <a name="Q2.3">How do I change the behavior of SCIP?</a>
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
      <a name="Q2.4"> I recognized that one special plugin works very poorly / very well for my problem and
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
         (see <a href="#Q2.3">previous question</a>).</li>
      <li>For <b class="grey">weakening</b>, you should just do the opposite operation, i.e.,
          reducing the values  you would raise for intensification and vice versa.</li>
    </ul>
    </p>
    <p class="q">
      <a name="Q2.5">How can I input a MINLP into SCIP?</a>
    </p>
    <p class="a">
      Mixed-integer nonlinear constrained programs can be input via the reader for OSiL and ZIMPL files and the AMPL interface.
      An interface to the modeling language GAMS is available in the <a href="https://projects.coin-or.org/GAMSlinks">COIN-OR/GAMSlinks</a> project.<br>
      Further, mixed-integer <em>polynomially</em> constrained programs can be input via the readers for PIP files.<br>
      Finally, the readers for *.mps and *.lp files support the CPLEX extensions for <em>quadratic</em> constraints and objective.
    </p>


<h4><a name="Section3">3. Using SCIP included in another source code</a></h4>
    <p class="q">
      <a name="Q3.1">How do I construct a problem instance in SCIP?</a>
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
      href="OBJ.shtml">doxygen documentation.</a>.
    </p>
    <p class="q">
      <a name="Q3.2">I already know a solution in advance, which I want to pass to SCIP.
      How do I do this?</a>
    </p>
    <p class="a">
      First you have to build your problem (at least all variables have to
      exist), then there are two options:
      <ul>
      <li>
        You have the solution in file which fits the solution format of
        SCIP, then you can use <code>SCIPreadSol()</code> to pass that
        solution to SCIP.
      </li>
      <li>
        You create a new SCIP primal solution candidate by
        calling <code>SCIPcreateSol()</code> and set all nonzero values by
        calling <code>SCIPsetSolVal()</code>. After that, you add this
        solution by calling <code>SCIPaddSol()</code> (the variable
        <code>stored</code> should be true afterwards, if your solution was
        added to solution candidate store) and then release it by calling
        <code>SCIPsolFree()</code>. Instead of adding and releasing
        sequentially, you can use <code>SCIPaddSolFree()</code> which
        tries to add the solution to the candidate store and
        free the solution afterwards.
      </li>
      </ul>
    </p>
    <p class="q">
      <a name="Q3.3">What operational stages of SCIP are there and are they important for me?</a>
    </p>
    <p class="a">
      There are ten different stages during a run of SCIP.
      <!--, a flow chart can be seen here: ????????. -->
      There are some methods which cannot be called in all stages, consider for example <code>SCIPtrySol()</code>
      (see <a href="#Q3.2">previous question</a>).
    </p>
    <p class="q">
      <a name="Q3.4">What is the difference between the original and the transformed problem?</a>
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
      <a name="Q3.5">Why do the names, e.g., in debug messages often differ from the ones I
      defined?</a>
    </p>
    <p class="a">
      This can have several reasons. Especially names of binary variables can get different prefixes and suffixes.  Each
      transformed variable and constraint (see <a href="#Q3.4">here</a>) gets a &quot;t_&quot; as prefix.  Apart from
      that, the meaning of original and transformed variables and constraints is identical.<br>

      General integers with bounds that differ just by 1 will be aggregated to binary variables which get the same name
      with the suffix &quot;_bin&quot; .  E.g. an integer variable <code>t_x</code> with lower bound 4 and upper bound 5
      will be aggregated to a binary variable <code>t_x_bin = t_x - 4</code>.<br>

      Variables can have negated counterparts, e.g. for a binary <code>t_x</code> its (also binary) negated would be
      <code>t_x_neg = 1 - t_x</code>.<br>

      The knapsack constraint handler is able to disaggregate its constraints to cliques, which are set packing
      constraints, and create names that consist of the knapsack's name and a suffix
      &quot;<code>_clq_&lt;int&gt;</code>&quot;.  E.g., a knapsack constraint
      <code>knap:&nbsp;x_1&nbsp;+&nbsp;x2&nbsp;+2 x_3&nbsp;&le;&nbsp;2</code> could be disaggregated to the set packing
      constraints <code>knap_clq_1:&nbsp;x_1&nbsp;+&nbsp;x_3&nbsp;&le;&nbsp;1</code> and
      <code>knap_clq_2:&nbsp;x_2&nbsp;+&nbsp;x_3&nbsp;&le;&nbsp;1</code>.
    </p>
    <p class="q">
      <a name="Q3.6">What is SCIP_CALL()? Do I need this?</a>
    </p>
    <p class="a">
      Yes, you do. SCIP_CALL() is a global define, which handles the return codes of all methods
      which return a SCIP_RETCODE and should therefore parenthesize each such method.
      SCIP_OKAY is the code which is returned if everything worked well;
      there are 17 different error codes, see <a href="type__retcode_8h.shtml">type_retcode.h</a>.
      Each method that calls methods which return a SCIP_RETCODE should itself return a SCIP_RETCODE.
      If this is not possible, use SCIP_CALL_ABORT() to catch the return codes of the methods.
      If you do not want to use this either, you have to do the exception handling
      (i.e. the case that the return code is not SCIP_OKAY) on your own.
    </p>
    <p class="q">
      <a name="Q3.7">I want to stop the solving process after a certain time. How can I do this?</a>
    </p>
    <p class="a">
      Limits are given by parameters in SCIP, for example <code>limits/time</code> for a time limit or
      <code>limits/nodes</code> for a node limit. 
      If you want to set a limit, you have to change these parameters.
      For example, for setting the time limit to one hour, you have to call 
      <code>SCIP_CALL( SCIPsetRealParam(scip, "limits/time", 3600) )</code>.
      In the interactive shell, you just enter <code>set limits time 3600</code>.
      For more examples, please have a look into heur_rens.c.
    </p>


<h4><a name="Section4">4. Using SCIP as a Branch-Cut-And-Price-Framework</a></h4>
    <p class="q">
      <a name="Q4.1">How do I start a project?</a>
    </p>
    <p class="a">
      See <a href="START.shtml">How to start a new project</a>.
    </p>
    <p class="q">
      <a name="Q4.2">What types of plugins can I add and how do I do this?</a>
    </p>
    <p class="a">
      See the <a href="index.shtml">doxygen documentation</a> for a list of
      plugin types. There is a HowTo for each of them.
    </p>
    <p class="q">
      <a name="Q4.3">When should I implement a constraint handler, when should I implement a separator?</a>
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
      &quot;nosubtour&quot; constraint handler in the TSP example
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
      <a name="Q4.4">Can I remove unnecessary display columns or&mdash;even better&mdash;add my own ones?
      Can I change the statistics displayed at the end of solving?</a>
    </p>
    <p class="a">
      Setting the status of a display column to 0 turns it off. E.g., type "set display memused status 0" in the
      interactive shell to disable the memory information column, or include the line <code>SCIPsetIntParam(scip,
      &quot;display/memused/status&quot;, 0)</code> into your source code.  Adding your own display column can be done
      by calling the <code>SCIPincludeDisp()</code> method, see the <a href="DISP.shtml">doxygen
      documentation</a>.<br> The statistic display, which is shown by "display statistics" and
      <code>SCIPprintStatistics()</code>, respectively, cannot be changed.
    </p>
    <p class="q">
      <a name="Q4.5">What do LP-rows look like in SCIP?</a>
    </p>
    <p class="a">
      Each row is of the form <i>lhs</i> &le; &Sigma;(<i>val</i>[<i>j</i>]&middot;<i>col</i>[<i>j</i>]) + <i>const</i>
      &le; <i>rhs</i>.  For now, <i>val</i>[<i>j</i>]&middot;<i>col</i>[<i>j</i>] can be interpreted as
      <i>a</i><sub><i>ij</i></sub>&middot;<i>x</i><sub><i>j</i></sub> (for the difference between columns and variables
      see <a href="#Q4.7">here</a>).  The constant is essentially needed for collecting the influence of presolving
      reductions like variable fixings and aggregations.<br> The <i>lhs</i> and <i>rhs</i> may take infinite values: a
      less-than inequality would have <i>lhs</i> = -&#8734;, and a greater-than inequality would have <i>rhs</i> =
      +&#8734;.  For equations <i>lhs</i> is equal to <i>rhs</i>.  An infinite left hand side can be recognized by
      <code>SCIPisInfinity(scip, -lhs)</code>, an infinite right hand side can be recognized by
      <code>SCIPisInfinity(scip, rhs)</code>.
    </p>
    <p class="q">
      <a name="Q4.6">How do I get the data of the current LP-relaxation?</a>
    </p>
    <p class="a">
      You can get all rows in the current LP-relaxation by calling <code>SCIPgetLPRowsData()</code>. The methods
      <code>SCIProwGetConstant()</code>, <code>SCIProwGetLhs()</code>, <code>SCIProwGetRhs()</code>,
      <code>SCIProwGetVals()</code>, <code>SCIProwGetNNonz()</code>, <code>SCIProwGetCols()</code> then give you
      information about each row, see <a href="#Q4.5">previous question</a>.
    </p>
    <p class="a">
      You get a columnwise representation by
      calling <code>SCIPgetLPColsData()</code>. The methods <code>SCIPcolGetLb()</code> and <code>SCIPcolGetUb()</code>
      give you the locally valid bounds of a column in the LP relaxation of the current branch-and-bound-node.
    </p>
    <p class="a">
      If you are interested in global information, you have to call <code>SCIPcolGetVar()</code> to get the variable
      associated to a column (see <a href="#Q4.7">next question</a>), which you can ask for global bounds via
      <code>SCIPvarGetLbGlobal()</code> and <code>SCIPvarGetUbGlobal()</code> as well as the type of the variable
      (binary, general integer, implicit integer, or continuous) by calling <code>SCIPvarGetType()</code>.  For more
      information, also see <a href="#Q4.8">this question</a>.
    </p>
    <p class="q">
      <a name="Q4.7">What is the difference between columns and variables, rows and constraints?</a>
    </p>
    <p class="a">
      The terms columns and rows always refer to the representation in the current LP-relaxation, variables and
      constraints to your global Constraint Integer Program.<br> Each column has an associated variable, which it
      represents, but not every variable must be part of the current LP-relaxation. E.g., it could be already fixed,
      aggregated to another variable, or be priced out if a column generation approach was implemented.
    </p>
    <p class="a">
      Each row has either been added to the LP by a constraint handler or by a cutting plane separator.  A constraint handler is able
      to, but does not need to, add one or more rows to the LP as a linear relaxation of each of its constraints.  E.g.,
      in the usual case (i.e. without using dynamic rows) the linear constraint handler adds one row to the LP for each
      linear constraint.
    </p>
    <p class="q">
      <a name="Q4.8">Are the variables and rows sorted in any particular order?</a>
    </p>
    <p class="a">The variable array which you get by <code>SCIPgetVars()</code> is internally sorted by variable types.
    The ordering is binary, integer, implicit integer and continuous variables, i.e., the binary variables are stored at
    position [0,...,nbinvars-1], the general integers at [nbinvars,...,nbinvars+nintvars-1], and so on. It holds that
    nvars = nbinvars + ninitvars + nimplvars + ncontvars.  There is no further sorting within these sections, as well as
    there is no sorting for the rows.  But each column and each row has a unique index, which can be obtained by
    <code>SCIPcolGetIndex()</code> and <code>SCIProwGetIndex()</code>, respectively.
    </p>
    <p class="q">
      <a name="Q4.9">When should I use which of the numerical comparison functions?</a>
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
    <p class="q">
      <a name="Q4.10">How do I solve an LP inside my SCIP plugin?</a>
    </p>
    <p class="a">
       If the LP is only a slightly modified version of the LP relaxation - changed variable bounds or objective
       coefficients - then you can use SCIP's diving mode: methods <code>SCIPstartDive()</code>,
       <code>SCIPchgVarLbDive()</code>, <code>SCIPsolveDiveLP()</code>, etc.

       Alternatively, SCIP's probing mode allows for a tentative depth first search in the tree and can solve the LP
       relaxations at each node: methods <code>SCIPstartProbing()</code>, <code>SCIPnewProbingNode()</code>,
       <code>SCIPfixVarProbing()</code>, etc. However, you cannot change objective coefficients or enlarge variable bounds
       in probing mode.<br/><br/>

       If you need to solve a separate LP, creating a sub-SCIP is not recommended because of the overhead involved
       and because dual information is not accessible (compare <a href="#Q4.5">here</a>). Instead you can use SCIP's LP
       interface. For this you should include <code>scip/lpi.h</code> and call the methods provided therein.
       Note that the LPi can be used independently from SCIP.
    </p>


<h4><a name="Section5">5. Specific questions about Column Generation and Branch-And-Price with SCIP</a></h4>
    <p class="q">
      <a name="Q5.1">Why are not all variables in the LP?</a>
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
      <a name="Q5.2">I only implemented one pricer, why is there a second one, called variable pricer?</a>
    </p>
    <p class="a">
      As described in the <a href="#Q5.1">previous question</a>, it may happen, that some variables are not in the
      current LP relaxation. Nevertheless, these variables still exist, and SCIP can calculate their reduced costs and
      add them to the LP again, if necessary. This is the job of the variable pricer. It is called before all other
      pricers.
    </p>
    <p class="q">
      <a name="Q5.3">How can I store branching decisions?</a>
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
      all the information needed (see the <a href="#Q5.4">next question</a>).
    </p>
    <p class="a">
      The third way is to use an event handler, which is described <a href="#Q5.5">here</a>.
    </p>
    <p class="q">
      <a name="Q5.4">I want to store some information at the nodes and update my pricer's data structures
                    when entering a new node. How can I do that?</a>
    </p>
    <p class="a">
      This can be done by creating a new constraint handler with constraint data that can store the information and
      do/undo changes in the pricer's data structures.
    </p>
    <p class="a">
      In general, all methods of the constraint handler (check, enforcing, separation, ...) should be empty (which means
      that they always return the status SCIP_FEASIBLE for the fundamental callbacks), just as if all constraints of this type
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
      of the constraint handler. Since, in general, each node is only propagated once, in this case you will have to check
      in your CONSACTIVE method whether new variables were added after your last propagation of this node. If this is
      the case, you will have to mark this node for repropagation by <code>SCIPrepropagateNode()</code>.
    </p>
    <p class="a">
      You can look into the constraint handler of the coloring problem (examples/Coloring/src/cons_storeGraph.c) to get
      an example of a constraint handler that does all these things.
    </p>
    <p class="q">
      <a name="Q5.5">How can an event handler help me with my branching?</a>
    </p>
    <p class="a">
      An event handler can watch for events like local bound changes on variables. So, if your 
      pricer wants to be informed whenever a local bound of a certain variable changes, add an 
      event handler, catch the corresponding events of the variable, and in the event handler's 
      execution method adjust the data structures of your pricer accordingly.
    </p>

    <p class="q">
      <a name="Q5.6">How can I add locally valid variables to the problem in my branch-and-price code?</a>
    </p>
    <p class="a">
      Variables in SCIP are always added globally. If you want to add them locally, because they are forbidden
      in another part of the branch-and-bound-tree, you should ensure that they are locally fixed to 0
      in all subtrees where they are not valid. A description of how this can be done is given <a href="#Q5.4">here</a>.
    </p>

    <p class="q">
      <a name="Q5.7">My pricer generates the same column twice. How can I solve this problem?</a>
    </p>
    <p class="a">
      First check whether your pricing is correct. Are there upper bounds on variables that you have forgotten to take
      into account? If your pricer cannot cope with variable bounds other than 0 and infinity, you have to mark
      all constraints containing priced variables as modifiable, and you may have to disable reduced cost
      strengthening by setting propagating/rootredcost/freq to -1.
    </p>
    <p class="a">  
      If your pricer works correctly and makes sure that the same column is added at most once in one pricing round, this
      behavior is probably caused by the PRICER_DELAY property of your pricer.
    </p>
    <p class="a">
      If it is set to FALSE, the following may have happened: The variable pricer (see <a href="#Q5.2">this question</a>)
      found a variable with negative dual feasibility that was not part of the current LP relaxation and added it to the
      LP. In the same pricing round, your own pricer found the same column and created a new variable for it. This might
      happen, since your pricer uses the same dual values as the variable pricer. To avoid this behavior, set
      PRICER_DELAY to TRUE, so that the LP is reoptimized after the variable pricer added variables to the LP. You can
      find some more information about the PRICER_DELAY property at <a
      href="PRICER.shtml#PRICER_PROPERTIES">How to add variable pricers </a>.
    </p>

    <p class="q">
      <a name="Q5.8">Which default plugins should be deactivated in order to get a working branch-and-price code?</a>
    </p>
    <p class="a">
      In most cases, you should deactivate separators, since cutting planes that are added to your master problem may
      destroy your pricing problem. Additionally, it may be necessary to deactivate some presolvers, mainly the dual
      fixing presolver. This can be done by not including these plugins into SCIP, namely by not calling
      <code>SCIPincludeSepaXyz()</code> and <code>SCIPincludePresolXyz()</code> in your own plugins-including
      files. Alternatively, you can set the parameters maxrounds and maxroundsroot to zero for all separators and
      maxrounds to zero for the presolvers.
    </p>

    <p class="q">
      <a name="Q5.9">What are the lazy bounds for variables in SCIP and what do I need them for?</a>
    </p>
    <p class="a">
      In many Branch-and-Price applications, you have binary variables, but you do not want to impose
      upper bounds on these variables in the LP relaxation, because the upper bound is already implicitly enforced
      by the problem constraints and the objective. If the upper bounds are explicitly added to the LP,
      they lead to further dual variables, which may be hard to take into account in the pricing problem.
    </p>
    <p class="a">  
      There are two possibilities for how to solve this problem. First, you could change the binary variables
      to general integer variables, if this does not change the problem. However, if you use special linear
      constraints like set partitioning/packing/covering, you can only add binary variables to these constraints.
    </p>
    <p class="a">  
      In order to still allow the usage of these types of constraints in a branch-and-price approach, the
      concept of lazy bounds was introduced in SCIP 2.0.
      For each variable, you can define lazy upper and lower bounds, i.e. bounds, that are implicitly enforced by
      constraints and objective. SCIP adds variable bounds to the LP only if the bound is tighter than the
      corresponding lazy bound.
      Note that lazy bounds are explicitly put into and removed from the LP when starting and ending diving mode,
      respectively. This is needed because changing the objective in diving might reverse the implicitly enforced
      bounds.
    </p>
    <p class="a">  
      For instance, if you have set partitioning constraints in your problem, you can define variables contained in these
      constraints as binary and set the lazy upper bound to 1, which allows you to use the better propagation methods
      of the setppc constraint handler compared to the linear constraint handler without taking care about upper
      bounds on variables in the master.
    </p>

    <p class="q">
      <a name="Q5.10">Can I stop the pricing process before the master problem is solved to optimality?</a>
    </p>
    <p class="a">
      In a column generation approach, you usually have to solve the master problem to optimality; otherwise,
      its objective function value is not a valid dual bound. However, there is a way in SCIP to
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

    <p class="q">
      <a name="Q5.11">How can I delete variables?</a>
    </p>
    <p class="a">
      SCIP features the functionality to delete variables from the problem when performing branch-and-price.
      This feature is still in a beta status and can be activated by switching the parameters
      <code>pricing/delvars</code> and <code>pricing/delvarsroot</code> to TRUE in order to allow deletion
      of variables at the root node and at all other nodes, respectively.
      Furthermore, variables have to be marked to be deletable by <code>SCIPvarMarkDeletable()</code>, which has to be done
      before adding the variable to the problem. Then, after a node of the branch-and-bound-tree is processed,
      SCIP automatically deletes variables from the problem that were created at the current node and whose corresponding
      columns were already removed from the LP. Note that due to the way SCIP stores basis information, it is not possible to
      completely delete a variable that was created at another node than the current node.
      You might want to change the parameters <code>lp/colagelimit</code>, <code>lp/cleanupcols</code>, and <code>lp/cleanupcolsroot</code>,
      which have an impact on when and how fast columns are removed from the LP.
    </p>
    <p class="a">
      Constraint handlers support a new callback function that deletes variables from constraints in which they were marked to be deleted.
      Thus, when using automatic variable deletion, you should make sure that all used constraint handlers implement this callback.
      By now, the linear, the set partitioning/packing/covering and the knapsack constraint handler support this callback, which should
      be sufficient for most branch-and-price applications. Note that set covering constraints can be used instead of logicor constraints.
    </p>
    <p class="a">
      Instead of deleting a variable completely, you can also remove it from the problem by either fixing the variable to zero
      using <code>SCIPfixVar()</code>, which fixes the variable globally or using <code>SCIPchgVarUbNode()</code> and
      <code>SCIPchgVarLBNode()</code>, which changes the bounds only for the current subtree.
    </p>

<h4><a name="Section6">6. Specific questions about the copy functionality in SCIP</a></h4>
    <p class="q">
      <a name="Q6.1">What is <code>SCIPcopy()</code>?</a>
    </p>
    <p class="a">
      The functionality of copying a SCIP model was added in SCIP version 2.0.0. It gives the possibility to
      generate a copy of the current SCIP model. This functionality is of interest, for example, in large
      neighborhood heuristics (such as heur_rens.c). They can now easily copy the complete problem and fix a certain
      set of variables to work on a reasonable copy of the original problem.
    </p>
    <p class="q">
      <a name="Q6.2">How do I get a copy of a variable or a constraint?</a>
    </p>
    <p class="a">
      For the variables and constraints there are the methods <code>SCIPgetVarCopy()</code> and 
      <code>SCIPgetConsCopy()</code> which provide a copy for a variable or a constraint, respectively.
    </p>
    <p class="q">
      <a name="Q6.3">What does the <code>valid</code> pointer in the copy callback of the constraint handler and variable pricer mean?</a>
    </p>
    <p class="a">
      SCIP would like to know if the copied problem is a one to one copy. That is, if all problem defining objects
      were successfully copied. If this is the case, all reductions made in the copy can be transferred to the
      original instance. The problem defining objects in SCIP are the constraint handlers and the variable pricers.
    </p>
          <!-- FAQ Punkt: Verhaeltnis zu CPLEX callable library, concert, Coin... ????????? -->
          <!-- VBC tool, visualisierung ????????? -->
</div>

