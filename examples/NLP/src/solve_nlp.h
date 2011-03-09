/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _NLP_H_
#define _NLP_H_
#include "scip/scip.h"
/**
 * The optimization problem has the following form :
 * minimize f(x) = c^T * x
 *  s.t. lhs_i <= c_i(x) <= rhs_i i=1,...,m
 *       l_i <= x_i <= u_i
 *  where, c_i(x) can be a linear constraint or a ploynomial constraint     
 */

typedef struct MonomialIpoptTag
{
   int nvars;                         /* the number of variables in this monomial */
   int * indicies;	              /* an array of length nvars, indicies[i] means the i-th variable in this monomial is
                                         x[indecies[i]] */
   SCIP_Real * power;                 /* an array of lenght nvars, power[i] means the power of the i-th variable in this monomial,
                                         i.e. the power of x[indicies[x]], NOTICE x is the memeber of struct NLP */
   SCIP_Real coefficient;             /* this is the coefficient of this monomial, you can but usually should not give a 0
                                         this member */
} MonomialIpopt;

typedef struct PolynomialIpoptTag
{
   int nMonomials;               /* this is the number of monomials in this polynomial */
   MonomialIpopt ** monomials;   /* a monomial pointer array of length nMonomials */
} PolynomialIpopt;

typedef struct NLPTag {
   int nvars;                   /* the dimension of variables, this is the total number of the vars */
   int nbinvars;                /* the number of binary variables */
   int nintvars;                /* the number of integer variables */
   int nimplvars;               /* the number of implicite variables */
   int ncontvars;               /* the number of continuous variables */
   int nactivevars;             /* the number of active variables after presolving */
   int nnonactivevars;          /* the number of nonactive variables after presolving */
   int nfixed;                  /* the number of variables which are fixed to constant after presolving */
   int naggr;                   /* the number of variables which are aggregated */
   int nmultaggr;               /* the number of variables which are multiaggregated */
   int nnegation;               /* the number of variables which are the negations of other variables */
   
   SCIP_Real * c;               /* a array of dimension nvars, coefficients of the linear objective function */
   SCIP_Real * x;               /* a array of dimension nvars, variable values used in ipopt, when call ipopt from scip, it store the inintial values for the variables , when ipopt improve the solution it pass the better solution back to scip */
   SCIP_Real * lb;              /* a array of dimension nvars, the lower bound of variables */
   SCIP_Real * ub;              /* a array of dimension nvars,the upper bound of variables */
                                /* NOTICE, for the above member c,x,lb,ub,has the same order, the order is decide by x */  
   
   int mcons;                   /* the number of constriants */
   int m_LP;                    /* the number of linear constriats */
   int m_NLP;                   /* the number of nonlinear constriants , these three member has the relationship mcons = m_LP + m_NLP */
   
   int * nnonz;                 /* it is an array of length m_LP, nnonz[i] tells you how many variables in linear constraint i */
   int ** jCols;                /* it is a pointer array of length m_LP, jCols[j] tells you the indices array of linear constraint j
                                   for j=0, ..., m_LP, here the indices mean, suppose the coefficient of x[i] for i=0,...,nvars-1,
                                   is c[i] != 0 in linear constraint i, then this i should be in the array of jCols[j], and the lengh
                                   array jCols[j] is the nonzero element in the linear constrait j, jCols[j][i] means the j-th linear
                                   constriant's i-th variables is x[jCols[j][i] and we suggest these integer
                                   number should instored in some order, like monotonically increasing */
   SCIP_Real ** values;	        /* it is a pointer array of length m_LP, values[j][i] tells you the i-th nonzero coefficient of the
                                   j-th linear  constriant, i.e. the coefficient of x[jCols[j][i]] in linear constraint j */

   PolynomialIpopt** polynomials; /* it is a array, each element is a polynomial data structure */

   SCIP_Real * lhs;         	  /* it is an array of size mcons, means the left hand side value of each constraint */
   SCIP_Real * rhs;          	  /* it is an array of size mcons, means the right hand side value of each constriant */
} NLP;

extern
int PrintMonomial(MonomialIpopt* monomial);

extern
int PrintPolynomial(PolynomialIpopt* polynomial);

extern
int PrintNLP(SCIP* scip,  NLP* nlp );

extern
int Callipopt(SCIP* scip, NLP* nlp_data );

#endif

