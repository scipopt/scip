<!-- this file is used to generate the local doxygen documentation -->
<!-- using make doc from within scip or soplex -->

<style>
.reveal:hover {
    text-shadow: 1px 1px 1px #777;
}
.reveal {
    color:#06c;
    padding-top: 1em;
    cursor: pointer;
    margin: 5px 0;
}
.answer {
    background-color: #fff;
    padding-left:   1em;
}
</style>

<h4>General Questions about SCIP</h4><ol>  <li>
    <div id="whatisscip" class="reveal">
      What is SCIP?    </div>
    <div id="whatisscip_ans" class="answer">
       <p>
SCIP is a constraint integer program solver. It can be used as a framework for <b>branch-cut-and-price</b>
and contains all necessary plugins to serve as a <b>standalone solver</b> for MIP, MINLP, and PBO.
</p>
<ul>
<li>You can use the <a href="http://scip.zib.de/#download">precompiled binaries</a>
to solve MIPs, MINLPs, or PBO problems. They can read

<ol>
<li>CIP - file format (own constraint based file format)</li>
<li><a href="http://people.sc.fsu.edu/~jburkardt/data/cnf/cnf.html">CNF - file format</a></li>
<li><a href="http://www.g12.cs.mu.oz.au/minizinc/downloads/doc-1.5/flatzinc-spec.pdf">FZN - file format</a></li>
<li><a href="http://publib.boulder.ibm.com/infocenter/cosinfoc/v12r3/topic/ilog.odms.cplex.help/Content/Optimization/Documentation/Optimization_Studio/_pubskel/ps_reffileformatscplex2159.html">(R)LP - file formats</a></li>
<li><a href="http://en.wikipedia.org/wiki/MPS_%28format%29">MPS - file format</a></li>
<li><a href="http://www.cril.univ-artois.fr/PB12/format.pdf">OPB/WBO - file formats</a></li>
<li><a href="http://www.coin-or.org/OS/OSiL.html">OSiL format for MINLP</a></li>
<li><a href="http://polip.zib.de/pipformat.php">PIP - file format</a></li>
<li>SOL - file format (solution reader, xml-solutions are also supported)</li>
<li><a href="http://zimpl.zib.de/download/zimpl.pdf">ZPL - file format</a></li>
</ol>
See also the page on <a href="./group__FILEREADERS.shtml">file readers</a>.
</li>
<li>You can use SCIP as a subroutine for solving MIPs and more general constraint integer programs from your own source code.</li>
<li>You can use SCIP as a framework in which you implement your own plugins.</li>
<li>You can use SCIP in any combination of the three purposes above.</li>
</ul>
<p>
This FAQ contains separate sections covering each of these usages of SCIP. It further considers specific questions for some features.
</p>
    </div>
  </li>
    <li>
    <div id="whenusescip" class="reveal">
      When should I use SCIP?    </div>
    <div id="whenusescip_ans" class="answer">
       <p>
If you are either looking for a fast non-commercial MIP-solver or for a
branch-cut-and-price-framework in which you can directly implement your own methods
&mdash; which can even be for more general purposes than MIP.
</p>
    </div>
  </li>
    <li>
    <div id="licensefaq" class="reveal">
      I heard something about licenses. Do I have to pay for using SCIP?    </div>
    <div id="licensefaq_ans" class="answer">
       <p>
As long as you use it for academic, non-commercial purposes: No.
This will not change. For the other cases, check the explanation of the
<a href="http://scip.zib.de/#license">ZIB academic license</a> and always feel free to ask us.
If you want to use SCIP commercially, please write an e-mail to koch@zib.de.
</p>
    </div>
  </li>
    <li>
    <div id="howtogetstarted" class="reveal">
      How do I get started?    </div>
    <div id="howtogetstarted_ans" class="answer">
       <p>
An easy way is to use the SCIP-binaries and call SCIP from a shell, see <a href="./SHELL.shtml">here</a> for a tutorial.
For that, you just have to download one of the precompiled binaries from the
<a href="http://scip.zib.de/#download">download section</a>, or the zipped source code and compile it
with your favorite settings. This is described in detail in the <code><a href="./INSTALL.shtml">INSTALL</a></code> file in the SCIP main directory.
</p>
<p>
Another way is to use SCIP as a solver integrated into your own program source code.
See the directories &quot;examples/MIPsolver/&quot; and &quot;examples/Queens/&quot;
for simple examples and <a href="#howtocreateproblem">this point</a>.
</p>
<p>
A third way is to implement your own plugins into SCIP.
This is explained in the HowTos for all plugin types, which you can find in the
<a href="./index.shtml">doxygen documentation</a>.
See also <a href="./START.shtml">How to start a new project</a>.
</p>
    </div>
  </li>
    <li>
    <div id="extrasoftware" class="reveal">
      Do I need any extra software?    </div>
    <div id="extrasoftware_ans" class="answer">
       <p>Unless you want to use SCIP as a pure CP-Solver (see <a href="#scipascpsolver">here</a>),
you need an underlying LP-Solver installed and linked to the libraries
(see the <code><a href="./INSTALL.shtml">INSTALL</a></code> file in the SCIP root directory).
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
<p>
We also provide some <a href="http://scip.zib.de/#download">precompiled binaries</a>. Besides that,
you might need a modeling language like <a href="http://zimpl.zib.de/">ZIMPL</a> to generate *.mps
or *.lp files. ZIMPL files can also directly be read by SCIP.
You can download a package which includes SCIP, SoPlex and ZIMPL <a href="http://scip.zib.de/#download">here</a>.
</p>
<p>
If you want to use SCIP for mixed integer nonlinear programming, you might want to use an underlying NLP solver (e.g.,
<a href="https://projects.coin-or.org/Ipopt">Ipopt</a>) and an expression interpreter (e.g., <a href="http://www.coin-or.org/CppAD/">CppAD</a>).
</p>
    </div>
  </li>
    <li>
    <div id="installationproblems" class="reveal">
      I have installation problems. What can I do?    </div>
    <div id="installationproblems_ans" class="answer">
       <p>
Read the <code><a href="./INSTALL.shtml">INSTALL</a></code> file in the SCIP root directory. It contains hints of how to get around problems. You can also try
the binaries available on the <a href="http://scip.zib.de">SCIP page</a>. <br>
If you want to use <a href="http://soplex.zib.de/">SoPlex</a> as the underlying LP-solver, you can try the
following:
</p>
<p>
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
<p>
If you do not have all of these libraries, read the <code><a href="./INSTALL.shtml">INSTALL</a></code> file
in the SCIP Optimization Suite directory.
</p>
<p>
As a summary, the call of <code>make ZIMPL=false ZLIB=false READLINE=false</code> should work on most systems.
If you have any problems while using an LP-solver different from SoPlex,
please read the SCIP <code><a href="./INSTALL.shtml">INSTALL</a></code> file first.
</p>
<p>
If you encounter compilation problems dealing with the <code>explicit</code> keyword you can either try a newer
compiler or set the flags <code>LEGACY=true</code> for SoPlex and <code>SPX_LEGACY=true</code> for SCIP.
</p>
    </div>
  </li>
    <li>
    <div id="versionchange" class="reveal">
      I changed to a new version of SCIP and now compiling breaks with some error messages which I don't understand.
Do you have a general hint on that?    </div>
    <div id="versionchange_ans" class="answer">
       <p>
Maybe the parameters of a function in SCIP changed. Relevant changes between version are listed below.
</p>
<ul>
<li> <a href="./CHG1.shtml">Interface changes between SCIP 0.9 and SCIP 1.0</a></li>
<li> <a href="./CHG2.shtml">Interface changes between SCIP 1.0 and SCIP 1.1</a></li>
<li> <a href="./CHG3.shtml">Interface changes between SCIP 1.1 and SCIP 1.2</a></li>
<li> <a href="./CHG4.shtml">Interface changes between SCIP 1.2 and SCIP 2.0</a></li>
<li> <a href="./CHG5.shtml">Interface changes between SCIP 2.0 and SCIP 2.1</a></li>
<li> <a href="./CHG6.shtml">Interface changes between SCIP 2.1 and SCIP 3.0</a></li>
<li> <a href="./CHG7.shtml">Interface changes between SCIP 3.0 and SCIP 3.1</a></li>
<li> <a href="./CHG8.shtml">Interface changes between SCIP 3.1 and SCIP 3.2</a></li>
</ul>
    </div>
  </li>
    <li>
    <div id="howtodebug" class="reveal">
      How can I debug in SCIP?    </div>
    <div id="howtodebug_ans" class="answer">
       <p>
Compile SCIP in debug mode: <code>make OPT=dbg</code>. Put the binary into a debugger, e.g., <a href="http://www.gnu.org/s/gdb/">gdb</a>
and let it run again. If you get an impression which component is causing the trouble, set <code>#define SCIP_DEBUG</code> as the first
line of the corresponding *.c file, recompile and let it run again. This will print debug messages from that piece of code.
Find a short debugging tutorial <a href="./DEBUG.shtml">here</a>.
</p>
    </div>
  </li>
    <li>
    <div id="errorhandling" class="reveal">
      SCIP prints error messages, aborts, produces segmentation faults, or just behaves strangely. What should I do?    </div>
    <div id="errorhandling_ans" class="answer">
       <p>
<a href="#howtodebug">See above.</a> Often, the asserts that show up in debug mode already help to clarify misunderstandings
and suggest fixes, if you were calling SCIP functions in an unexpected manner. If you, however, have the impression
that there is a bug in SCIP, <a href="https://www.edom.mi.uni-erlangen.de/cgi/scipweb/submit.pl">send us a bug report</a>
and include a log file and ideally a backtrace from a debugger.
</p>
    </div>
  </li>
    <li>
    <div id="searchdocumentation" class="reveal">
      I would like to check whether some functionality is implemented in SCIP.
How does the naming of the methods work? Where do I find the most common methods?    </div>
    <div id="searchdocumentation_ans" class="answer">
       <p>
For an explanation of the naming see the
<a href="./CODE.shtml">coding style guidelines</a>. <br>
The methods you might want to use should normally be found in scip.h or in pub_*.h,
see also <a href="./DOC.shtml">here</a>.
In the doxygen documentation of <a href="./scip_8h.shtml">scip.h</a>,
the methods are ordered by topics.
</p>
    </div>
  </li>
    <li>
    <div id="scipascpsolver" class="reveal">
      Can I use SCIP as a pure CP-Solver?    </div>
    <div id="scipascpsolver_ans" class="answer">
       <p>
Yes. SCIP can be used as a pure CP-Solver by typing
<code>set emphasis cpsolver</code> in the shell or by using the function <code>SCIPsetEmphasis()</code>.
Furthermore, you can compile SCIP without any LP-Solver by <code>make LPS=none</code>.
See <a href="#howtochangebehaviour">here</a> for more on changing the behavior of SCIP.
</p>
    </div>
  </li>
    <li>
    <div id="scipaslpsolver" class="reveal">
      Can I use SCIP as a pure LP-Solver?    </div>
    <div id="scipaslpsolver_ans" class="answer">
       <p>
Since LPs are only special types of MIPs and CIPs, the principal answer is yes. If you feed a pure LP to SCIP, it
will first apply presolving and then hand this presolved problem to the underlying LP solver. If the LP is solved to
optimality, you can query the optimal solution values as always.
</p>
<p>
However, there are certain limitations to this: Dual multipliers and reduced costs are not accessible. If the LP
turns out to be infeasible, you cannot obtain the Farkas proof. Strong LP duality theory does
not apply to MIP or CIP and there is no special functionality implemented in SCIP for the case
that the solved problem is an LP.
</p>
<p>
Hence, if you need more, "LP specific", information than the primal solution, you are better off using an LP-Solver
directly. If you are using the SCIP Optimization Suite, you could, e.g., use the included LP solver
<a href=http://soplex.zib.de>SoPlex</a>.  If you want to solve an LP not from the command line, but within your
C/C++ program, you could also use SCIP's LP-Interface, see also <a href="#howtosolvelps">here</a>.
</p>
    </div>
  </li>
    <li>
    <div id="minlptypes" class="reveal">
      Which kind of MINLPs are supported by SCIP?    </div>
    <div id="minlptypes_ans" class="answer">
       <p>
SCIP supports nonlinear constraints of the form lhs &le; f(x) &le; rhs, where the function f(x) is an algebraic
expression that can be represented as expression tree. Such an expression tree has constants and variables as terminal
nodes and operands as non-terminal nodes.
Expression operands supported by SCIP include addition, subtraction, multiplication, division, exponentiation and
logarithm. Trigonometric functions are not yet supported by SCIP.
</p>
<p>
Nonlinear objective functions are not supported by SCIP and must be modeled as constraint function.
<em>Note</em>, that the support for non-quadratic nonlinear constraints is still in a <b>BETA</b>-stadium and not yet as
robust as the rest of SCIP.
Missing bounds on nonlinear variables and tiny or huge coefficients can easily lead to numerical problems,
which can be avoided by careful modeling.
</p>
    </div>
  </li>
    <li>
    <div id="howtousemakefiles" class="reveal">
      How can I use the SCIP makefiles?    </div>
    <div id="howtousemakefiles_ans" class="answer">
       <p>
See <a href="./MAKE.shtml">Makefiles</a>.
</p>
    </div>
  </li>
    <li>
    <div id="staticvsshared" class="reveal">
      What is this business with .a and .so libraries in the directory <code>lib/</code>?    </div>
    <div id="staticvsshared_ans" class="answer">
       <p>
When SCIP builds the binary it needs to link with the corresponding libraries of the LP-solver. There are usually
two ways to distribute a library. In the first (.a) the library is linked <em>statically</em> to SCIP; this means
that all of its information is packed into the binary. In the second way (.so) the library is a <em>shared</em>
library. In this case, the code of the library is not inserted into the binary itself, but is loaded at
runtime. This has the advantage that the binaries are smaller, but it comes at the cost that you have to make
sure that the library is found at runtime (for most systems it suffices to put the path of the library into the
<code>LD_LIBRARY_PATH</code> variable).<br/><br/> <em>Note:</em> Depending on the system shared or static libraries
are standard. If both links are present in the <code>lib</code> directory, the linker chooses which version to take
(on newer Linux systems this is usually the shared version). If you do not want this to happen you have to delete
the link that is not intended.
</p>
    </div>
  </li>
    <li>
    <div id="sharedlib" class="reveal">
      Can I compile SCIP as a shared library?    </div>
    <div id="sharedlib_ans" class="answer">
       <p>
You can use the <code>SHARED=true</code> option
when making SCIP. This will generate the libraries of SCIP in shared format. The binary then
also uses this form. Note that the path to the <code>lib</code> directory of SCIP is used to
locate the libraries. If you want to move the libraries, you might to set the
<code>LD_LIBRARY_PATH</code> environment variable to include the new path.

<!--
If you have a different system or compiler (which supports shared libraries) you might try to
copy the file <code>make.linux.x86.gnu.opt-shared</code> and modify it according to your
needs (see also <code><a href="./MAKE.shtml">Makefiles</a></code>).
-->
If you are using your own building system: The "magic" changes are the <code>-FPIC</code>
compiler/linker option and the <code>-Wl,-rpath</code> option.
</p>
    </div>
  </li>
    <li>
    <div id="scipgetvarsol" class="reveal">
      The methods SCIPgetVarSol() and SCIPvarGetSol() seem to have the same
functionality. Which one should I use?    </div>
    <div id="scipgetvarsol_ans" class="answer">
       <p>
In fact, there is a slight difference: <code>SCIPvarGetSol()</code> is also able to return pseudo solution values.
If you do not have an idea, what pseudo solutions are, <code>SCIPgetVarSol()</code> should be just fine. This
should be the only case of 'duplicate methods'. If you find, however, another one, please contact us.
</p>
    </div>
  </li>
    <li>
    <div id="visualizebranchandbound" class="reveal">
      Is there a way to visualize the branch and bound tree?    </div>
    <div id="visualizebranchandbound_ans" class="answer">
       <p>
Yes, currently there are two tree visualizations supported, the <a href="http://www.informatik.uni-koeln.de/ls_juenger/research/vbctool">vbctool</a>,
and the <a href="https://projects.coin-or.org/CoinBazaar/wiki/Projects/BAK">Branch-and-Bound Analysis Kit</a>.
The first comes with a viewer that has an option to uncover the nodes one-by-one (each time you
hit the space key). Additional node information such as its lower bound, depth, and number are accessed through a
context menu. The BAK is a command line tool written in the Python programming language. It offers several
noninteractive visualizations of the tree at a user-defined frequency.
</p>
<p>
For using one of these tools, SCIP lets you define file names <code>set visual vbcfilename somefilename.vbc</code> and
<code>set visual bakfilename somefilename.dat</code>.
</p>
<p>
For those who want to use the step-by-step functionality of vbctool,
it is necessary to use a time-step counter for the visualization instead of the real time. The corresponding parameter
is changed via <code>set visual realtime FALSE</code>.
</p>
For users of the callable library, the corresponding parameters are called "visual/bakfilename", "visual/vbcfilename", and
"visual/realtime".
</p>
    </div>
  </li>
  </ol><br/>
<h4>Using SCIP as a standalone solver</h4><ol>  <li>
    <div id="outputtoolarge" class="reveal">
      The output is too wide for my terminal window. What can I do?    </div>
    <div id="outputtoolarge_ans" class="answer">
       <p>
In the interactive shell you can set the width of the output with the following command
<code>set display width</code> followed by an appropriate number.
See also the <a href="#displaycolumns">next question</a>.
</p>
    </div>
  </li>
    <li>
    <div id="displaycolumns" class="reveal">
      What do the cryptic abbreviations for the columns mean which are displayed
during the solving process of SCIP?    </div>
    <div id="displaycolumns_ans" class="answer">
       <p>
Type <code>display display</code> in the interactive shell to get an explanation of them.<br>
By the way: If a letter appears in front of a display row, it indicates,
which heuristic found the new primal bound, a star representing an integral LP-relaxation.<br>
Typing <code>display statistics</code> after finishing or interrupting the solving process gives
you plenty of extra information about the solving process.<br>
(Typing <code>display heuristics</code> gives you a list of the heuristics including their letters.)
</p>
    </div>
  </li>
    <li>
    <div id="whereareuserparameters" class="reveal">
      Why does SCIP claim that it could not find the user parameters "scip.set"? Where can I get such a file?    </div>
    <div id="whereareuserparameters_ans" class="answer">
       <p>
SCIP comes with default settings that are automatically active when you start the interactive shell.
However, you have the possibility to save customized settings via the <code>set save</code> and <code>set diffsave</code> commands.
Both commands will prompt you to enter a file name and save either all or customized parameters only to the specified file.

A user parameter file that you save as "scip.set" has a special meaning; Whenever you invoke SCIP from a directory containing a file named "scip.set",
the settings therein overwrite the default settings. For more information about customized settings, see the <a href="./SHELL.shtml">
Tutorial on the interactive shell</a>.

Settings files can become incompatible with later releases if we decide to rename/delete a parameter.
Information about this can be found in the CHANGELOG for every release, see also <a href="#versionchange">this related question</a>.

</p>
    </div>
  </li>
    <li>
    <div id="howtochangebehaviour" class="reveal">
      How do I change the behavior of SCIP?    </div>
    <div id="howtochangebehaviour_ans" class="answer">
       <p>
You can switch the settings for all presolving, heuristics, and separation plugins to three different modes
via the <code>set {presolving, heuristics, separation} emphasis</code> parameters in the interactive shell.
<code>off</code> turns off the respective type of plugins, <code>fast</code> chooses settings that lead to less time spent
in this type of plugins, decreasing their impact, and <code>aggressive</code> increases the impact
of this type of plugins.
You can combine these general settings for cuts, presolving, and heuristics arbitrarily.<br>
<code>display parameters</code> shows you which settings currently differ from their default,
<code>set default</code> resets them all.
Furthermore, there are complete settings that can be set by <code>set emphasis</code>, i.e. settings for pure
feasibility problems, solution counting, and CP like search.
</p>
    </div>
  </li>
    <li>
    <div id="displaystatistics" class="reveal">
      How can I learn more about/from the presolve reasoning SCIP applies to my combinatorial optimization problem?    </div>
    <div id="displaystatistics_ans" class="answer">
       <p>
You can look at the statistics (type <code>display statistics</code> in the
interactive shell or call SCIPprintStatistics() when using SCIP as a library).
This way you can see, which of the presolvers, propagators or constraint
handlers performed the reductions.
Then, add a <code>#define SCIP_DEBUG</code> as first line of the corresponding *.c file
in src/scip (e.g. cons_linear.c or presol_probing.c)
Recompile and run again. You will get heaps of information now.
Looking into the code and documentation of the corresponding plugin and ressources from the literature
helps with the investigation.
</p>
    </div>
  </li>
    <li>
    <div id="toggleplugin" class="reveal">
      I recognized that one special plugin works very poorly / very well for my problem and
I want to disable it / weaken its influence / intensify its influence. How do I do this?    </div>
    <div id="toggleplugin_ans" class="answer">
       <ul>
<li>
For using a <b class="grey">non-default branching rule or node selection strategy</b> as
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
These lists give information about the different priorities.
</li>
<li>
If you want to completely <b class="grey">disable a heuristic or a separator</b> you have
to set its frequency to -1 and the sepafreq to -1 for separation by constraint handlers,
respectively. The commands looks like this:
<ul>
<li><code> SCIP&gt; set heuristics &lt;name of a heuristic&gt; freq -1</code></li>
<li><code> SCIP&gt; set separators &lt;name of a separator&gt; freq -1</code></li>
<li><code> SCIP&gt; set constraints &lt;name of a constraint handler&gt; sepafreq -1</code></li>
</ul>
</li>
<li>
For <b class="grey">disabling a presolver</b>, you have to set its maxrounds parameter to 0.
<ul>
<li><code> SCIP&gt; set presolvers &lt;name of a presolver&gt; maxrounds 0</code></li>
</ul>
</li>
<li>
If you want to <b class="grey">intensify</b> the usage of a <b class="grey">heuristic</b>,
you can reduce its frequency to some smaller, positive value, and/or raise the quotient and
offset values (maxlpiterquot for diving heuristics, nodes for LNS heuristics).
<ul>
<li><code> SCIP&gt; set heuristics &lt;name of a heuristic&gt; freq &lt;some value&gt;</code></li>
<li><code> SCIP&gt; set heuristic &lt;name of a diving heuristic&gt; maxlpiterquot &lt;some value&gt;</code></li>
<li><code> SCIP&gt; set heuristic &lt;name of a LNS heuristic&gt; nodesquot &lt;some value&gt;</code></li>
</ul>
</li>
<li>
For <b class="grey">intensifying</b> the usage of a <b class="grey">separator</b>,
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
Compare the parameters of the heuristic/separator in the appropriate <a href="#howtochangebehaviour">aggressive setting</a>
</li>
<li>
For <b class="grey">weakening</b>, you should just do the opposite operation, i.e.,
reducing the values  you would raise for intensification and vice versa.
</li>
</ul>
    </div>
  </li>
    <li>
    <div id="howtoeextendshell" class="reveal">
      How can I use my own functions in the interactive shell/extend the set of available interactive shell commands?    </div>
    <div id="howtoeextendshell_ans" class="answer">
       <p>
If you want to keep the interactive shell functionality, you could add a
dialog handler, that introduces a new SCIP shell command that
</p>
<ol type="a">
<li> solves the problem and calls your function afterwards or </li>
<li>checks whether the stage is SOLVED and only calls your function.</li>
</ol>
<p>
Search <code>SCIPdialogExecOptimize</code> in src/scip/dialog_default.c to see how
the functionality of the "optimize" command is invoked.
Also, in src/scip/cons_countsols.c, you can see an example of a dialog handler
being added to SCIP.
If this is the way you go, please check the <a href="./DIALOG.shtml">How to add dialogs section</a>
of the doxygen documentation.
</p>
    </div>
  </li>
    <li>
    <div id="howtoinputminlp" class="reveal">
      How can I input a MINLP into SCIP?    </div>
    <div id="howtoinputminlp_ans" class="answer">
       <p>
Mixed-integer nonlinear constrained programs can be input via the reader for OSiL and ZIMPL files and the AMPL interface.
An interface to the modeling language GAMS is available in the <a href="https://projects.coin-or.org/GAMSlinks">COIN-OR/GAMSlinks</a> project.
Further, mixed-integer <em>polynomially</em> constrained programs can be input via the readers for PIP files.
Finally, the readers for *.mps and *.lp files support the CPLEX extensions for <em>quadratic</em> constraints and objective.
</p>
    </div>
  </li>
  </ol><br/>
<h4>Using SCIP included in another source code</h4><ol>  <li>
    <div id="howtocreateproblem" class="reveal">
      How do I construct a problem instance in SCIP?    </div>
    <div id="howtocreateproblem_ans" class="answer">
       <p>
First you have to create a SCIP object via <code>SCIPcreate()</code>, then you start to build the problem via
<code>SCIPcreateProb()</code>.  Then you create variables via <code>SCIPcreateVar()</code> and add them to the
problem via <code>SCIPaddVar()</code>.
</p>
<p>
The same has to be done for the constraints. For example, if you want to fill in the rows of a general MIP, you
have to call <code>SCIPcreateConsLinear()</code>, <code>SCIPaddConsLinear()</code> and additionally
<code>SCIPreleaseCons()</code> after finishing.  If all variables and constraints are present, you can initiate
the solution process via <code>SCIPsolve()</code>.
</p>
<p>
Make sure to also call <code>SCIPreleaseVar()</code> if you do not need the variable pointer anymore.  For an
explanation of creating and releasing objects, please see the <a
href="./OBJ.shtml">doxygen documentation.</a>.
</p>
    </div>
  </li>
    <li>
    <div id="howtopasssolution" class="reveal">
      I already know a solution in advance, which I want to pass to SCIP.
How do I do this?    </div>
    <div id="howtopasssolution_ans" class="answer">
       <p>
First you have to build your problem (at least all variables have to
exist), then there are two options:
</p>
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
    </div>
  </li>
    <li>
    <div id="scipstages" class="reveal">
      What operational stages of SCIP are there and are they important for me?    </div>
    <div id="scipstages_ans" class="answer">
       <p>
There are fourteen different stages during a run of SCIP.
<!--,
a flow chart can be seen here: ????????. -->
There are some methods which cannot be called in all stages, consider for example <code>SCIPtrySol()</code>
(see <a href="#howtopasssolution">previous question</a>).
</p>
    </div>
  </li>
    <li>
    <div id="originalvstransformed" class="reveal">
      What is the difference between the original and the transformed problem?    </div>
    <div id="originalvstransformed_ans" class="answer">
       <p>
Before the solving process starts, the original problem is copied.
This copy is called &quot;transformed problem&quot;, and all modifications during the presolving
and solving process are only applied to the transformed problem.<br>
This has two main advantages: first, the user can also modify the problem after partially solving it.
All modifications done by SCIP (presolving, cuts, variable fixings) during the partial solving process will
be deleted together with the transformed problem, the user can modify the original problem and restart solving.
Second, the feasibility of solutions is always tested on the original problem!
</p>
    </div>
  </li>
    <li>
    <div id="namingissues" class="reveal">
      Why do the names, e.g., in debug messages often differ from the ones I
defined?    </div>
    <div id="namingissues_ans" class="answer">
       <p>
This can have several reasons. Especially names of binary variables can get different prefixes and suffixes.  Each
transformed variable and constraint (see <a href="#originalvstransformed">here</a>) gets a &quot;t_&quot; as prefix.  Apart from
that, the meaning of original and transformed variables and constraints is identical.
</p>
<p>
General integers with bounds that differ just by 1 will be aggregated to binary variables which get the same name
with the suffix &quot;_bin&quot; .  E.g. an integer variable <code>t_x</code> with lower bound 4 and upper bound 5
will be aggregated to a binary variable <code>t_x_bin = t_x - 4</code>.
</p>
<p>
Variables can have negated counterparts, e.g. for a binary <code>t_x</code> its (also binary) negated would be
<code>t_x_neg = 1 - t_x</code>.
</p>
<p>
The knapsack constraint handler is able to disaggregate its constraints to cliques, which are set packing
constraints, and create names that consist of the knapsack's name and a suffix
&quot;<code>_clq_&lt;int&gt;</code>&quot;.  E.g., a knapsack constraint
<code>knap:&nbsp;x_1&nbsp;+&nbsp;x2&nbsp;+2 x_3&nbsp;&le;&nbsp;2</code> could be disaggregated to the set packing
constraints <code>knap_clq_1:&nbsp;x_1&nbsp;+&nbsp;x_3&nbsp;&le;&nbsp;1</code> and
<code>knap_clq_2:&nbsp;x_2&nbsp;+&nbsp;x_3&nbsp;&le;&nbsp;1</code>.
</p>
    </div>
  </li>
    <li>
    <div id="scipcall" class="reveal">
      What is SCIP_CALL()? Do I need this?    </div>
    <div id="scipcall_ans" class="answer">
       <p>
Yes, you do. SCIP_CALL() is a global define, which handles the return codes of all methods
which return a SCIP_RETCODE and should therefore parenthesize each such method.
SCIP_OKAY is the code which is returned if everything worked well;
there are 17 different error codes, see <a href="./type__retcode_8h.shtml">type_retcode.h</a>.
Each method that calls methods which return a SCIP_RETCODE should itself return a SCIP_RETCODE.
If this is not possible, use SCIP_CALL_ABORT() to catch the return codes of the methods.
If you do not want to use this either, you have to do the exception handling
(i.e. the case that the return code is not SCIP_OKAY) on your own.
</p>
    </div>
  </li>
    <li>
    <div id="howtointerrupt" class="reveal">
      I want to stop the solving process after a certain time. How can I do this?    </div>
    <div id="howtointerrupt_ans" class="answer">
       <p>
Limits are given by parameters in SCIP, for example <code>limits/time</code> for a time limit or
<code>limits/nodes</code> for a node limit.
If you want to set a limit, you have to change these parameters.
For example, for setting the time limit to one hour, you have to call
<code>SCIP_CALL( SCIPsetRealParam(scip, "limits/time", 3600) )</code>.
In the interactive shell, you just enter <code>set limits time 3600</code>.
For more examples, please have a look into heur_rens.c.
</p>
    </div>
  </li>
  </ol><br/>
<h4>Using SCIP as a Branch-Cut-And-Price-Framework</h4><ol>  <li>
    <div id="howtostartbcpproject" class="reveal">
      How do I start a project?    </div>
    <div id="howtostartbcpproject_ans" class="answer">
       <p>
See <a href="./START.shtml">How to start a new project</a>.
</p>
    </div>
  </li>
    <li>
    <div id="howtoaddplugins" class="reveal">
      What types of plugins can I add and how do I do this?    </div>
    <div id="howtoaddplugins_ans" class="answer">
       <p>
See the <a href="./modules.shtml">doxygen documentation</a> for a list of
plugin types. There is a <a href="./index.shtml#HOWTOADD">HowTo</a> for each of them.
</p>
    </div>
  </li>
    <li>
    <div id="conshdlrvsseparator" class="reveal">
      When should I implement a constraint handler, when should I implement a separator?    </div>
    <div id="conshdlrvsseparator_ans" class="answer">
       <p>
This depends on whether you want to add constraints or only cutting planes. The main
difference is that constraints can be &quot;model constraints&quot;, while cutting planes are only
additional LP rows that strengthen the LP relaxation.
A model constraint is a constraint that is important for the feasibility of the integral
solutions. If you delete a model constraint, some infeasible integral vectors would
suddenly become feasible in the reduced model.
A cutting plane is redundant w.r.t. integral solutions. The set of feasible integral
vectors does not change if a cutting plane is removed. You can, however, relax this
condition slightly and add cutting planes that do cut off feasible solutions, as long as
at least one of the optimal solutions remains feasible.
</p>
<p>
You want to use a constraint handler in the following cases:
</p>
<ol>
<li>
Some of your feasibility conditions can not be expressed by existing constraint types
(e.g., linear constraints), or you would need too many of them. For example, the
&quot;nosubtour&quot; constraint in the TSP is equivalent to exponentially many linear constraints.
Therefore, it is better to implement a &quot;nosubtour&quot; constraint handler that can inspect
solutions for subtours and generate subtour elimination cuts and others (e.g., comb
inequalities) to strengthen the LP relaxation.</li>
<li>Although you can express your feasibility condition by a reasonable number of existing
constraint types, you can represent and process the condition in a more efficient way. For
example, it may be that you can, due to your structural knowledge, implement a stronger or
faster domain propagation or find tighter cutting planes than what one could do with the
sum of the individual &quot;simple&quot; constraints that model the feasibility condition.
</li>
</ol>
<p>
You want to use a cutting plane separator in the following cases:
</p>
<ol>
<li>
You have a general purpose cutting plane procedure that can be applied to any MIP. It
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
<p>
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
    </div>
  </li>
    <li>
    <div id="toggledisplaycolumns" class="reveal">
      Can I remove unnecessary display columns or&mdash;even better&mdash;add my own ones?
Can I change the statistics displayed at the end of solving?    </div>
    <div id="toggledisplaycolumns_ans" class="answer">
       <p>
Setting the status of a display column to 0 turns it off. E.g., type <code>set display memused status 0</code> in the
interactive shell to disable the memory information column, or include the line <code>SCIPsetIntParam(scip,
&quot;display/memused/status&quot;, 0)</code> into your source code.  Adding your own display column can be done
by calling the <code>SCIPincludeDisp()</code> method, see the <a href="./DISP.shtml">doxygen
documentation</a>.<br> The statistic display, which is shown by <code>display statistics</code> and
<code>SCIPprintStatistics()</code>, respectively, cannot be changed.
</p>
    </div>
  </li>
    <li>
    <div id="lprows" class="reveal">
      What do LP-rows look like in SCIP?    </div>
    <div id="lprows_ans" class="answer">
       <p>
Each row is of the form <i>lhs</i> &le; &Sigma;(<i>val</i>[<i>j</i>]&middot;<i>col</i>[<i>j</i>]) + <i>const</i>
&le; <i>rhs</i>.  For now, <i>val</i>[<i>j</i>]&middot;<i>col</i>[<i>j</i>] can be interpreted as
<i>a</i><sub><i>ij</i></sub>&middot;<i>x</i><sub><i>j</i></sub> (for the difference between columns and variables
see <a href="#columnsvsvariables">here</a>).  The constant is essentially needed for collecting the influence of presolving
reductions like variable fixings and aggregations.<br> The <i>lhs</i> and <i>rhs</i> may take infinite values: a
less-than inequality would have <i>lhs</i> = -&#8734;, and a greater-than inequality would have <i>rhs</i> =
+&#8734;.  For equations <i>lhs</i> is equal to <i>rhs</i>.  An infinite left hand side can be recognized by
<code>SCIPisInfinity(scip, -lhs)</code>, an infinite right hand side can be recognized by
<code>SCIPisInfinity(scip, rhs)</code>.
</p>
    </div>
  </li>
    <li>
    <div id="lprelaxationdata" class="reveal">
      How do I get the data of the current LP-relaxation?    </div>
    <div id="lprelaxationdata_ans" class="answer">
       <p>
You can get all rows in the current LP-relaxation by calling <code>SCIPgetLPRowsData()</code>. The methods
<code>SCIProwGetConstant()</code>, <code>SCIProwGetLhs()</code>, <code>SCIProwGetRhs()</code>,
<code>SCIProwGetVals()</code>, <code>SCIProwGetNNonz()</code>, <code>SCIProwGetCols()</code> then give you
information about each row, see <a href="#lprows">previous question</a>.
</p>
<p>
You get a columnwise representation by
calling <code>SCIPgetLPColsData()</code>. The methods <code>SCIPcolGetLb()</code> and <code>SCIPcolGetUb()</code>
give you the locally valid bounds of a column in the LP relaxation of the current branch-and-bound-node.
</p>
<p>
If you are interested in global information, you have to call <code>SCIPcolGetVar()</code> to get the variable
associated to a column (see <a href="#columnsvsvariables">next question</a>), which you can ask for global bounds via
<code>SCIPvarGetLbGlobal()</code> and <code>SCIPvarGetUbGlobal()</code> as well as the type of the variable
(binary, general integer, implicit integer, or continuous) by calling <code>SCIPvarGetType()</code>.  For more
information, also see <a href="#searchdocumentation">this question</a>.
</p>
    </div>
  </li>
    <li>
    <div id="columnsvsvariables" class="reveal">
      What is the difference between columns and variables, rows and constraints?    </div>
    <div id="columnsvsvariables_ans" class="answer">
       <p>
The terms columns and rows always refer to the representation in the current LP-relaxation, variables and
constraints to your global Constraint Integer Program.<br> Each column has an associated variable, which it
represents, but not every variable must be part of the current LP-relaxation. E.g., it could be already fixed,
aggregated to another variable, or be priced out if a column generation approach was implemented.
</p>
<p>
Each row has either been added to the LP by a constraint handler or by a cutting plane separator.  A constraint handler is able
to, but does not need to, add one or more rows to the LP as a linear relaxation of each of its constraints.  E.g.,
in the usual case (i.e. without using dynamic rows) the linear constraint handler adds one row to the LP for each
linear constraint.
</p>
    </div>
  </li>
    <li>
    <div id="variableandrowsorting" class="reveal">
      Are the variables and rows sorted in any particular order?    </div>
    <div id="variableandrowsorting_ans" class="answer">
       <p>
The variable array which you get by <code>SCIPgetVars()</code> is internally sorted by variable types.
The ordering is binary, integer, implicit integer and continuous variables, i.e., the binary variables are stored at
position [0,...,nbinvars-1], the general integers at [nbinvars,...,nbinvars+nintvars-1], and so on. It holds that
nvars = nbinvars + ninitvars + nimplvars + ncontvars.  There is no further sorting within these sections, as well as
there is no sorting for the rows.  But each column and each row has a unique index, which can be obtained by
<code>SCIPcolGetIndex()</code> and <code>SCIProwGetIndex()</code>, respectively.
</p>
    </div>
  </li>
    <li>
    <div id="feasibilitycomparison" class="reveal">
      When should I use which of the numerical comparison functions?    </div>
    <div id="feasibilitycomparison_ans" class="answer">
       <p>
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
    </div>
  </li>
    <li>
    <div id="howtosolvelps" class="reveal">
      How do I solve an LP inside my SCIP plugin?    </div>
    <div id="howtosolvelps_ans" class="answer">
       <p>
If the LP is only a slightly modified version of the LP relaxation - changed variable bounds or objective
coefficients - then you can use SCIP's diving mode: methods <code>SCIPstartDive()</code>,
<code>SCIPchgVarLbDive()</code>, <code>SCIPsolveDiveLP()</code>, etc.
</p>
<p>
Alternatively, SCIP's probing mode allows for a tentative depth first search in the tree and can solve the LP
relaxations at each node: methods <code>SCIPstartProbing()</code>, <code>SCIPnewProbingNode()</code>,
<code>SCIPfixVarProbing()</code>, etc. However, you cannot change objective coefficients or enlarge variable bounds
in probing mode.
</p>
<p>
If you need to solve a separate LP, creating a sub-SCIP is not recommended because of the overhead involved
and because dual information is not accessible (compare <a href="#scipaslpsolver">here</a>). Instead you can use SCIP's LP
interface. For this you should include <code>lpi/lpi.h</code> and call the methods provided therein.
Note that the LPI can be used independently from SCIP.
</p>
    </div>
  </li>
  </ol><br/>
<h4>Specific questions about Column Generation and Branch-And-Price with SCIP</h4><ol>  <li>
    <div id="whattoexpect" class="reveal">
      What can I expect when using SCIP as a Branch-Cut-and-Price framework?    </div>
    <div id="whattoexpect_ans" class="answer">
       <p>
If you want to use SCIP as a branch-and-price framework, you normally need to implement a reader to read in your
problem data and build the problem, a pricer to generate new columns, and a branching rule to do the branching
(see also  <a href="#howtostorebranchingdecisions">this question</a> to see how to store branching decisions, if needed).
SCIP takes care about everything else,
for example the branch-and-bound tree management and LP solving including storage of warmstart bases. Moreover,
many of SCIP's primal heuristics will be used and can help improve your primal bound.
However, this also comes with a few restrictions: You are not allowed to change the objective function coefficients
of variables during the solving process, because that means that previously computed dual bounds might have to be updated.
This prevents the use of dual variable stabilization techniques based on a (more or less strict) bounding box in the dual.
We are working on making this possible and recommend to use a weighted sum stabilization approach until then.
Another point that SCIP does for you is the dynamic removal of columns from the LP due to aging (see also the next two questions).
However, due to the way simplex bases are stored in SCIP, columns can only be removed at the same node where they were created.
</p>
    </div>
  </li>
    <li>
    <div id="wherearethevariables" class="reveal">
      Why are not all variables in the LP?    </div>
    <div id="wherearethevariables_ans" class="answer">
       <p>
With <code>SCIPgetLPColsData()</code> you can obtain the columns of the current LP relaxation. It is correct that
not all variables are necessarily part of the current LP relaxation. In particular, in branch-and-price the
variables generated at one node in the tree are not necessarily included in the LP relaxation of a different node
(e.g., if the other node is not a descendant of the first node). But even if you are still at the same node or at
a descendant node, SCIP can remove columns from the LP, if they are 0 in the LP relaxation. This dynamic column
deletion can be avoided by setting the "removable" flag to FALSE in the <code>SCIPcreateVar()</code> call.
</p>
    </div>
  </li>
    <li>
    <div id="variablepricer" class="reveal">
      I only implemented one pricer, why is there a second one, called variable pricer?    </div>
    <div id="variablepricer_ans" class="answer">
       <p>
As described in the <a href="#wherearethevariables">previous question</a>, it may happen, that some variables are not in the
current LP relaxation. Nevertheless, these variables still exist, and SCIP can calculate their reduced costs and
add them to the LP again, if necessary. This is the job of the variable pricer. It is called before all other
pricers.
</p>
    </div>
  </li>
    <li>
    <div id="howtostorebranchingdecisions" class="reveal">
      How can I store branching decisions?    </div>
    <div id="howtostorebranchingdecisions_ans" class="answer">
       <p>
This is a very common problem in Branch-And-Price, which you can deal nicely with using SCIP.
There are basically three different options.
The first one is to add binary variables to the problem that encode branching decisions. Then constraints should
be added that enforce the corresponding branching decisions in the subtrees.
</p>
<p>
If you have complex pricer data like a graph and need to update it after each branching decision,
you should introduce "marker constraints" that are added to the branching nodes and store
all the information needed (see the <a href="#pricerandnodedata">next question</a>).
</p>
<p>
The third way is to use an event handler, which is described <a href="#eventhdlrduringbranching">here</a>.
</p>
    </div>
  </li>
    <li>
    <div id="pricerandnodedata" class="reveal">
      I want to store some information at the nodes and update my pricer's data structures
when entering a new node. How can I do that?    </div>
    <div id="pricerandnodedata_ans" class="answer">
       <p>
This can be done by creating a new constraint handler with constraint data that can store the information and
do/undo changes in the pricer's data structures.
<p>
Once you have such a constraint handler, just create constraints of this type and add them to the child nodes of
your branching by <code>SCIPaddConsNode()</code>. Make sure to set the "stickingatnode" flag to TRUE in order to
prevent SCIP from moving the constraint around in the tree.
</p>
</p>
In general, all methods of the constraint handler (check, enforcing, separation, ...) should be empty (which means
that they always return the status SCIP_FEASIBLE for the fundamental callbacks), just as if all constraints of this type
are always feasible. The important callbacks are the CONSACTIVE and CONSDEACTIVE methods for communicating the
constraints along the active path to your pricer, and the CONSDELETE callback for deleting data of constraints at nodes
which became obsolete.
<p>
The CONSACTIVE method is always called when a node is entered on which the constraint has been added. Here, you
need to apply the changes to your pricing data structures. The CONSDEACTIVE method will be called if the node is
left again. Since the CONSACTIVE and CONSDEACTIVE methods of different constraints are always called in a
stack-like fashion, this should be exactly what you need.
</p>
<p>
All data of a constraint need to be freed by implementing an appropriate CONSDELETE callback.
</p>
<p>
If you need to fix variables for enforcing your branching decision, this can be done in the propagation callback
of the constraint handler. Since, in general, each node is only propagated once, in this case you will have to check
in your CONSACTIVE method whether new variables were added after your last propagation of this node. If this is
the case, you will have to mark this node for repropagation by <code>SCIPrepropagateNode()</code>.
</p>
<p>
You can look into the constraint handler of the coloring problem (examples/Coloring/src/cons_storeGraph.c) to get
an example of a constraint handler that does all these things.
</p>
    </div>
  </li>
    <li>
    <div id="eventhdlrduringbranching" class="reveal">
      How can an event handler help me with my branching?    </div>
    <div id="eventhdlrduringbranching_ans" class="answer">
       <p>
An event handler can watch for events like local bound changes on variables. So, if your
pricer wants to be informed whenever a local bound of a certain variable changes, add an
event handler, catch the corresponding events of the variable, and in the event handler's
execution method adjust the data structures of your pricer accordingly.
</p>
    </div>
  </li>
    <li>
    <div id="localvariables" class="reveal">
      How can I add locally valid variables to the problem in my branch-and-price code?    </div>
    <div id="localvariables_ans" class="answer">
       <p>
Variables in SCIP are always added globally. If you want to add them locally, because they are forbidden
in another part of the branch-and-bound-tree, you should ensure that they are locally fixed to 0
in all subtrees where they are not valid. A description of how this can be done is given <a href="#pricerandnodedata">here</a>.
</p>
    </div>
  </li>
    <li>
    <div id="doublecolumns" class="reveal">
      My pricer generates the same column twice. How can I solve this problem?    </div>
    <div id="doublecolumns_ans" class="answer">
       <p>
First check whether your pricing is correct. Are there upper bounds on variables that you have forgotten to take
into account? If your pricer cannot cope with variable bounds other than 0 and infinity, you have to mark
all constraints containing priced variables as modifiable, and you may have to disable reduced cost
strengthening by setting propagating/rootredcost/freq to -1.
</p>
<p>
If your pricer works correctly and makes sure that the same column is added at most once in one pricing round, this
behavior is probably caused by the PRICER_DELAY property of your pricer.
</p>
<p>
If it is set to FALSE, the following may have happened: The variable pricer (see <a href="#variablepricer">this question</a>)
found a variable with negative dual feasibility that was not part of the current LP relaxation and added it to the
LP. In the same pricing round, your own pricer found the same column and created a new variable for it. This might
happen, since your pricer uses the same dual values as the variable pricer. To avoid this behavior, set
PRICER_DELAY to TRUE, so that the LP is reoptimized after the variable pricer added variables to the LP. You can
find some more information about the PRICER_DELAY property at <a
href="./PRICER.shtml#PRICER_PROPERTIES">How to add variable pricers </a>.
</p>
    </div>
  </li>
    <li>
    <div id="problematicplugins" class="reveal">
      Which default plugins should be deactivated in order to get a working branch-and-price code?    </div>
    <div id="problematicplugins_ans" class="answer">
       <p>
In most cases you should deactivate separators since cutting planes that are added to your master problem may
destroy your pricing problem. Additionally, it may be necessary to deactivate some presolvers, mainly the dual
fixing presolver. This can be done by not including these plugins into SCIP, namely by not calling
<code>SCIPincludeSepaXyz()</code> and <code>SCIPincludePresolXyz()</code> in your own plugins-including
files. Alternatively, you can set the parameters maxrounds and maxroundsroot to zero for all separators and
maxrounds to zero for the presolvers.
</p>
    </div>
  </li>
    <li>
    <div id="whatarelazybounds" class="reveal">
      What are the lazy bounds for variables in SCIP and what do I need them for?    </div>
    <div id="whatarelazybounds_ans" class="answer">
       <p>
In many Branch-and-Price applications, you have binary variables, but you do not want to impose
upper bounds on these variables in the LP relaxation, because the upper bound is already implicitly enforced
by the problem constraints and the objective. If the upper bounds are explicitly added to the LP,
they lead to further dual variables, which may be hard to take into account in the pricing problem.
</p>
<p>
There are two possibilities for how to solve this problem. First, you could change the binary variables
to general integer variables, if this does not change the problem. However, if you use special linear
constraints like set partitioning/packing/covering, you can only add binary variables to these constraints.
</p>
<p>
In order to still allow the usage of these types of constraints in a branch-and-price approach, the
concept of lazy bounds was introduced in SCIP 2.0.
For each variable, you can define lazy upper and lower bounds, i.e. bounds, that are implicitly enforced by
constraints and objective. SCIP adds variable bounds to the LP only if the bound is tighter than the
corresponding lazy bound.
Note that lazy bounds are explicitly put into and removed from the LP when starting and ending diving mode,
respectively. This is needed because changing the objective in diving might reverse the implicitly enforced
bounds.
</p>
<p>
For instance, if you have set partitioning constraints in your problem, you can define variables contained in these
constraints as binary and set the lazy upper bound to 1, which allows you to use the better propagation methods
of the setppc constraint handler compared to the linear constraint handler without taking care about upper
bounds on variables in the master.
</p>
    </div>
  </li>
    <li>
    <div id="terminatepricing" class="reveal">
      Can I stop the pricing process before the master problem is solved to optimality?    </div>
    <div id="terminatepricing_ans" class="answer">
       <p>
In a column generation approach, you usually have to solve the master problem to optimality; otherwise,
its objective function value is not a valid dual bound. However, there is a way in SCIP to
stop the pricing process earlier, called "early branching".
</p>
<p>
The reduced cost pricing method of a pricer has a result pointer that should be set each time the method is called.
In the usual case that the pricer either adds a new variable or ensures that there are no further
variables with negative dual feasibility, the result pointer should be set to SCIP_SUCCESS.
If the pricer aborts pricing without creating a new variable, but there might exist additional
variables with negative dual feasibility, the result pointer should be set to SCIP_DIDNOTRUN.
In this case, the LP solution will not be used as a lower bound.
Typically, early branching goes along with the computation of a Lagrangian bound in each pricing iteration.
The pricer store store this valid lower bound in the <code>lowerbound pointer</code> in order to update the lower
bound of the current node.
Since SCIP 3.1, it is even possible to state that pricing should be stopped early even though new variables were created
in the last pricing round. For this, the pricer has to set the <code>stopearly</code> pointer to TRUE.
</p>
    </div>
  </li>
    <li>
    <div id="howtodeletevariables" class="reveal">
      How can I delete variables?    </div>
    <div id="howtodeletevariables_ans" class="answer">
       <p>
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
<p>
Constraint handlers support a new callback function that deletes variables from constraints in which they were marked to be deleted.
Thus, when using automatic variable deletion, you should make sure that all used constraint handlers implement this callback.
By now, the linear, the set partitioning/packing/covering and the knapsack constraint handler support this callback, which should
be sufficient for most branch-and-price applications. Note that set covering constraints can be used instead of logicor constraints.
</p>
<p>
Instead of deleting a variable completely, you can also remove it from the problem by either fixing the variable to zero
using <code>SCIPfixVar()</code>, which fixes the variable globally or using <code>SCIPchgVarUbNode()</code> and
<code>SCIPchgVarLbNode()</code>, which changes the bounds only for the current subtree.
</p>
    </div>
  </li>
    <li>
    <div id="howtobranchonconstraints" class="reveal">
      How do I branch on constraints?    </div>
    <div id="howtobranchonconstraints_ans" class="answer">
       <p>
Constraint-based branching is rather straightforward to implement in
SCIP.  You have to add a new branching rule that uses the methods
SCIPcreateChild() and SCIPaddConsNode() in its branching callbacks.
A very good example for this is the Ryan/Foster branching rule that has
been implemented in the binpacking example from the <a href="./EXAMPLES.shtml">examples section</a>.
</p>
<p>
Sometimes it might be more appropriate to implement a constraint
handler instead of a branching rule. This is the case if, e.g., the added constraints alone
do NOT ensure integrality of the integer variables, or if you still want to use the available branching rules.
In the ENFOLP callback of your constraint handler, the
branching really happens. The integrality constraint
handler calls the branching rules within the ENFOLP callback. Give your constraint handler a positive enforcement priority
to trigger your constraint branching before the integrality constraint handler and perform the constraint branching.
</p>
    </div>
  </li>
  </ol><br/>
<h4>Specific questions about the copy functionality in SCIP</h4><ol>  <li>
    <div id="scipcopy" class="reveal">
      What is <code>SCIPcopy()</code>?    </div>
    <div id="scipcopy_ans" class="answer">
       <p>
The functionality of copying a SCIP model was added in SCIP version 2.0.0. It gives the possibility to
generate a copy of the current SCIP model. This functionality is of interest, for example, in large
neighborhood heuristics (such as heur_rens.c). They can now easily copy the complete problem and fix a certain
set of variables to work on a reasonable copy of the original problem.
</p>
    </div>
  </li>
    <li>
    <div id="howtocopyvariablesandconstraints" class="reveal">
      How do I get a copy of a variable or a constraint?    </div>
    <div id="howtocopyvariablesandconstraints_ans" class="answer">
       <p>
For the variables and constraints there are the methods <code>SCIPgetVarCopy()</code> and
<code>SCIPgetConsCopy()</code> which provide a copy for a variable or a constraint, respectively.
</p>
    </div>
  </li>
    <li>
    <div id="useofvalidpointer" class="reveal">
      What does the <code>valid</code> pointer in the copy callback of the constraint handler and variable pricer mean?    </div>
    <div id="useofvalidpointer_ans" class="answer">
       <p>
SCIP would like to know if the copied problem is a one to one copy. That is, if all problem defining objects
were successfully copied. If this is the case, all reductions made in the copy can be transferred to the
original instance. The problem defining objects in SCIP are the constraint handlers and the variable pricers.
</p>    </div>
  </li>
  </ol><br/>
<script>
$(".answer").hide();
$(".reveal").click(function() {
    var questionId = "#"+this.id+"_ans";
    $(questionId).toggle("fast");
});
</script>