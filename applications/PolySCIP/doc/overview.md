/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*        This file is part of the program PolySCIP                          */
/*                                                                           */
/*    Copyright (C) 2012-2018 Konrad-Zuse-Zentrum                            */
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
[PolySCIP] (http://polyscip.zib.de) is a solver for multi-criteria integer programming and multi-criteria linear 
programming problems.
In other words, it solves optimization problems of the form: 
\f{align*}{
\min / \max~ (c_1^\top x&, \ldots, c_k^\top x) \\
\mbox{s.t. } Ax &\leq b,\\
x &\in \mathbb{Z}^n \lor \mathbb{Q}^n,
\f}
where \f$ k \geq 2,~ A \in \mathbb{Q}^{m \times n},~ b \in
\mathbb{Q}^m \f$.

The name PolySCIP is composed of Poly (from the Greek
&pi;o&lambda;&upsilon;&sigmaf; meaning "many") and SCIP. 

Please note that there is a [user guide]
(http://polyscip.zib.de/#userguide) available describing command line parameters,
the file format and an easy generation of problem files
via the algebraic modeling language ZIMPL.
  
Installation
============

See the @ref INSTALL "INSTALL file".

Usage
=====

* After a succesful built, execute **polyscip -h** for detailed information about its command line arguments

* **polyscip probFile.mop** runs PolySCIP on the given problem file **probFile.mop** without any further parameters




