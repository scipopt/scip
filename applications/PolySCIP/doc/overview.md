/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*        This file is part of the program PolySCIP                          */
/*                                                                           */
/*    Copyright (C) 2012-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  PolySCIP is distributed under the terms of the ZIB Academic License.     */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with PolySCIP; see the file LICENCE.                               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   overview.md
 * @brief  overview page
 * @author Sebastian Schenker
*/

@mainpage Overview 

About
=====
[PolySCIP] (http://polyscip.zib.de) is a solver for multi-criteria integer programming as well as multi-criteria linear 
programming.
In other words, it solves optimization problems of the form: 
\f{align*}{
\min / \max~ &(c_1^T x, \ldots, c_k^T x) \\
\mbox{s.t. } Ax &\leq b,\\
x &\in \mathbb{Z}^n \lor \mathbb{Q}^n.
\f}
where \f$ k \geq 2,~ A \in \mathbb{Q}^{m \times n},~ b \in
\mathbb{Q}^m \f$.

The name PolySCIP is composed of Poly (from the Greek
&pi;o&lambda;&upsilon;&sigmaf; meaning "many") and SCIP. 
The long-term development goal of PolySCIP is to solve general
multi-criteria mixed-integer programming problems, i.e., problems of
the above mentioned form where \f$x \in \mathbb{Z}^{n_1} \times
\mathbb{Q}^{n_2}$ with $n_1 + n_2 = n \f$.

Please note that there is a more detailed user guide available on the PolySCIP [website] (http://polyscip.zib.de).
  
Installation
============

see INSTALL file 

Usage
=====

* After a succesful built, execute **polyscip -h** for detailed information about its command line arguments

* **polyscip probFile.mop** runs PolySCIP on the given problem file **probFile.mop** without any further parameters



