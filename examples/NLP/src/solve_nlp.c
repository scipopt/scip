/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "IpStdCInterface.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "solve_nlp.h"

static NLP* nlp_param = NULL;
/* Function Declarations */
static
Bool eval_f(Index n, Number* x, Bool new_x,
            Number* obj_value, UserDataPtr user_data);
static
Bool eval_grad_f(Index n, Number* x, Bool new_x,
                 Number* grad_f, UserDataPtr user_data);
static
Bool eval_g(Index n, Number* x, Bool new_x,
            Index m, Number* g, UserDataPtr user_data);
static
Bool eval_jac_g(Index n, Number *x, Bool new_x,
                Index m, Index nele_jac,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data);
static
Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data);

/* Interface called by SCIP */
extern
int PrintMonomial(MonomialIpopt* monomial)
{
   int i;
   for( i=0; i<monomial->nvars; ++i)
   {
      if( i==0 )
      {
         printf("%f*x[%d]^%f", monomial->coefficient, monomial->indicies[i], monomial->power[i]);
      }
      else
      {
         printf("x[%d]^%f", monomial->indicies[i], monomial->power[i]);
      }
   }
   return 0;
}

extern
int PrintPolynomial(PolynomialIpopt* polynomial)
{
   int i;

   for( i=0; i<polynomial->nMonomials; ++i)
   {
      if( i ==0 )
      {
         PrintMonomial( polynomial->monomials[i] );
      }
      else
      {
         printf("+");
         PrintMonomial( polynomial->monomials[i] );
      }
   }

   return 0;
}

extern
int PrintNLP( SCIP* scip, NLP* nlp)
{
   int i;
   int j;
   
   /* print the objective function */
   printf("This is the information from the PrintNLP, we print the data in NLP: \n");
   printf("There are %d variables, %d constraints, %d are linear constraits, %d are nonlinear constriants.\n", nlp->nvars, nlp->mcons, nlp->m_LP, nlp->m_NLP);
   
   printf("The objective function is: \n");
   for( i = 0; i<nlp->nvars; ++i)
   {
      if( i == 0 )
      {
         printf("%f*x[%d]",nlp->c[i],i);
      }
      else
      {
         printf("+%f*x[%d]",nlp->c[i],i);
      }
   }
   printf("\n");
   
   /** print the initial values for the variables */

   printf("These are the initial values for the variables: \n");
   for( i=0; i<nlp->nvars; ++i)
   {
      printf("x[%d] = %f \n ", i, nlp->x[i]);
   }
   printf("\n");
   
   /* print the boundary constraints */
   printf("The following are the bound for each variables: \n");
   for( i=0; i<nlp->nvars; ++i)
   {
      if( SCIPisInfinity(scip, -(nlp->lb[i])))
         printf("-Inf");
      else
         printf("%f", nlp->lb[i] );

      printf(" <= x[%d] <= ",i);

      if( SCIPisInfinity(scip, nlp->ub[i] ))
         printf("Inf \n");
      else
         printf("%f \n", nlp->ub[i] );
   }
   printf("\n");


   /* print the linear constraints */
   printf("There are %d linear constraints: \n", nlp->m_LP);
   for( i=0; i< nlp->m_LP; ++i)
   {
      printf("The %d -th linear constriants is: \n", i);

      if( SCIPisInfinity(scip, -(nlp->lhs[i])))
         printf("-Inf <= ");
      else
         printf("%f <= ", nlp->lhs[i] );

      for( j=0; j<nlp->nnonz[i]; ++j)
      {
         if( j==0 )
         {
            printf("%f*x[%d]",nlp->values[i][j],nlp->jCols[i][j]);
         }
         else
         {
            printf("+%f*x[%d]",nlp->values[i][j],nlp->jCols[i][j]);
         
         }
      }

      if( SCIPisInfinity(scip, -(nlp->rhs[i])))
         printf("<= Inf \n");
      else
         printf("<= %f \n", nlp->rhs[i] );
   }
   

   /* print the polynomial constraints */
   printf(" There are %d polynomial constraints: \n ", nlp->m_NLP);

   for( i=0; i<nlp->m_NLP; ++i)
   {
      printf( " The %d -th polynomial constraints( %d -th constraints ) is: \n", i, i+nlp->m_LP);

      if(SCIPisInfinity(scip,nlp->lhs[i+nlp->m_LP]))
         printf("-Inf <=");
      else
         printf("%f <=", nlp->lhs[i+nlp->m_LP]);

      PrintPolynomial(nlp->polynomials[i]);

      if(SCIPisInfinity(scip,nlp->rhs[i+nlp->m_LP]))
         printf("<= Inf \n");
      else
         printf("<= %f \n", nlp->rhs[i+nlp->m_LP]);
   }

   printf("That is the whole information for nlp passed to ipopt .\n");

   return 0;
      
}

extern
int Callipopt(SCIP* scip, NLP* nlp_data )     /**     Index Callipopt(SCIP* scip, NLP* nlp_data )*/
{
  /** get the problem information */
  assert( nlp_data != NULL );
  assert( scip != NULL );
  nlp_param = nlp_data;

  printf("hello, now we are in Callipopt :-)");
  PrintNLP(scip, nlp_param);


  Index n=-1;                          /* number of variables */
  Index m=-1;                          /* number of constraints */
  Number* x_L = NULL;                  /* lower bounds on x */
  Number* x_U = NULL;                  /* upper bounds on x */
  Number* g_L = NULL;                  /* lower bounds on g */
  Number* g_U = NULL;                  /* upper bounds on g */
  IpoptProblem nlp = NULL;             /* IpoptProblem */
  enum ApplicationReturnStatus status; /* Solve return code */
  Number* x = NULL;                    /* starting point and solution vector */
  Number* mult_x_L = NULL;             /* lower bound multipliers
  					  at the solution */
  Number* mult_x_U = NULL;             /* upper bound multipliers
  					  at the solution */
  Number obj;                          /* objective value */
  Index i;                             /* generic counter */

  /* indexing style for matrices */
  Index index_style = 0; /* C-style; start counting of rows and column
  			    indices at 0 */

  /* set the number of variables and allocate space for the bounds */
  n=nlp_param->nvars;
  x_L = (Number*)malloc(sizeof(Number)*n);
  x_U = (Number*)malloc(sizeof(Number)*n);
  /* set the values for the variable bounds */
  for( i=0; i<n; i++)
  {
       x_L[i] = nlp_param->lb[i];
       x_U[i] = nlp_param->ub[i];
  }


  /**  The alternative method
       x_L = nlp_param->lb;
       x_U = nlp_param->ub; */
  
  /** the following code is the example they gave to handle problem  HS_071
   * for (i=0; i<n; i++) {
   *     x_L[i] = 1.0;
   *         x_U[i] = 5.0;
   *           }
   */

  /* set the number of constraints and allocate space for the bounds */
  m=nlp_param->mcons;
  g_L = (Number*)malloc(sizeof(Number)*m);
  g_U = (Number*)malloc(sizeof(Number)*m);
  /* set the values of the constraint bounds */
  /** the following code is the example they gave to handle problem  HS_071
   *   g_L[0] = 25;
   *   g_U[0] = 2e19;
   *   g_L[1] = 40;
   *   g_U[1] = 40;
   */
  for( i=0; i<m; i++)
  {
       g_L[i] = nlp_param->lhs[i];
       g_U[i] = nlp_param->rhs[i];
  }

  /** The following are an alternative message
      g_L = nlp_param->lhs;
      g_U = nlp_param->rhs; */

  /* Number of nonzeros in the Jacobian of the constraints */
  /**  Index nele_jac = n*m; */
  Index  nele_jac = 0;
  Index row;
  Index j;
  Index col;
  printf("Count the number of nonzeros of the jacobian start from %d \n", nele_jac);
  for( row=0; row<nlp_param->m_LP; ++row)
  {
     for(j=0; j<nlp_param->nnonz[row]; ++j)
     {
        nele_jac += 1;
     }
  }
  printf("The number of nozeros from linear constriants is %d \n", nele_jac);
  for( row=nlp_param->m_LP; row<nlp_param->mcons; ++row)
  {
     for(col=0; col<n; ++col)
     {
        if(isVarInPolynomial(nlp_param->polynomials[row-nlp_param->m_LP],col))
        {
           nele_jac += 1;
        }
     }
  }
  
  printf("The number of nozeros from nonlinear constriants is %d \n", nele_jac);
  /* Number of nonzeros in the Hessian of the Lagrangian (lower or
     upper triangual part only) */
  /** Index nele_hess = 0.5*n*(n+1);*/
  Index nele_hess = 0;

  printf("Count the number of nonzeros of hessian, start from %d \n", nele_hess);
  for(row = 0; row < n; ++row)
  {
     for(col = 0; col < row; ++col)
     {
        for( j=0; j< nlp_param->m_NLP; ++j)
        {
           if(is2VarsInPolynomial(nlp_param->polynomials[j], row, col) == TRUE )
           {
              nele_hess++;
           }
           break;
        }
     }
  }
  printf("The nonzeros from the strict lower triangular part is %d \n", nele_hess);
  /* this part is for the diagonal */
  for(row = 0; row < n; ++row)
  {
     for( j=0; j< nlp_param->m_NLP; ++j)
     {
        if(isVarInPolynomialTwice(nlp_param->polynomials[j], row) == TRUE )
        {
           nele_hess++;
        }
        break;
     }
  }
  
  printf("The number of nonzeros for hessian from the diagnal part is %d \n", nele_hess);

  
 
  /* create the IpoptProblem */
  nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
                           index_style, &eval_f, &eval_g, &eval_grad_f,
                           &eval_jac_g, &eval_h);

  /* We can free the memory now - the values for the bounds have been
     copied internally in CreateIpoptProblem */
  free(x_L);
  free(x_U);
  free(g_L);
  free(g_U);

  /* Set some options.  Note the following ones are only examples,
     they might not be suitable for your problem. */
  AddIpoptNumOption(nlp, "tol", 1e-7);
  AddIpoptStrOption(nlp, "mu_strategy", "adaptive");
  AddIpoptStrOption(nlp, "output_file", "ipopt.out");

  /* allocate space for the initial point and set the values */
  x = (Number*)malloc(sizeof(Number)*n);
  /** the following code is the example they gave to handle problem  HS_071
   *   x[0] = 1.0;
   *   x[1] = 5.0;
   *   x[2] = 5.0;
   *   x[3] = 1.0;
   */

  for( i=0; i<n; i++)
  {
     x[i]=nlp_param->x[i];
  }
  
  /**  x = nlp_param->x; */

  /* allocate space to store the bound multipliers at the solution */
  mult_x_L = (Number*)malloc(sizeof(Number)*n);
  mult_x_U = (Number*)malloc(sizeof(Number)*n);

  /* solve the problem */
  status = IpoptSolve(nlp, x, NULL, &obj, NULL, mult_x_L, mult_x_U, NULL);

  if (status == Solve_Succeeded) {
    printf("\n\nSolution of the primal variables, x\n");
    for (i=0; i<n; i++) {
      printf("x[%d] = %e\n", i, x[i]);
    }

    printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
    for (i=0; i<n; i++) {
      printf("z_L[%d] = %e\n", i, mult_x_L[i]);
    }
    for (i=0; i<n; i++) {
      printf("z_U[%d] = %e\n", i, mult_x_U[i]);
    }

    printf("\n\nObjective value\n");
    printf("f(x*) = %e\n", obj);
  }

    for( i=0; i<n; i++)
  {
     nlp_data->x[i]=x[i];
  }
  

  /* free allocated memory */
  FreeIpoptProblem(nlp);
  free(x);
  free(mult_x_L);
  free(mult_x_U);

  return 0;
}

/** this is a local function for given a x and a monomial data to compute this monomial's value */
Number computeMonomialValue( MonomialIpopt * monomial, Number * x)
{
	Number val;
	Index i;
	val = monomial->coefficient;
	for( i=0; i<monomial->nvars; ++i)
	{
		val = val * pow(x[monomial->indicies[i]],monomial->power[i]);
	}
	return val;
}

/** this is a local function for given a x and a polynomial data to compute the value */
static
Number computePolynomialValue(PolynomialIpopt * polynomial, Number * x)
{
	Number val;
        Index i;
	val = 0;
	for( i=0; i<polynomial->nMonomials; ++i)
	{
	    val += computeMonomialValue(polynomial->monomials[i],x);
	}
	return val;
}

/* Function Implementations */
Bool eval_f(Index n, Number* x, Bool new_x,
            Number* obj_value, UserDataPtr user_data)
{
  Index i;
  assert(n == nlp_param->nvars);

  *obj_value = 0;
  for ( i=0; i<n; ++i)
  {
      *obj_value += nlp_param->c[i] * x[i];
  }

  /** the following code is the example they gave to handle problem  HS_071
   *   assert(n == 4);
   *   *obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
   */
  
  printf("The objective function value at\n");
  for( i=0; i<n; ++i )
  {
     printf("x[%d] = %f \n", i, x[i]);
  }
  printf("is %f, please input a char to continue \n",*obj_value);    
     

  return TRUE;
}

Bool eval_grad_f(Index n, Number* x, Bool new_x,
                 Number* grad_f, UserDataPtr user_data)
{
  Index i;
  assert( n == nlp_param->nvars );
  for( i=0; i<n; i++)
  {
	  grad_f[i]=nlp_param->c[i];
  }

  /** the following code is the example they gave to handle problem  HS_071
   *   assert(n == 4);
   *   grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
   *   grad_f[1] = x[0] * x[3];
   *   grad_f[2] = x[0] * x[3] + 1;
   *   grad_f[3] = x[0] * (x[0] + x[1] + x[2]);
   */

  printf("The gradient of the objective function value at\n");
  for( i=0; i<n; ++i )
  {
     printf("x[%d] = %f \n", i, x[i]);
  }

  for (i=0; i<n; ++i)
  {
     printf("grad_f[%d] = %f \n", i, grad_f[i]);
  }
  printf("please input a char to continue \n");    
     

  return TRUE;
}

Bool eval_g(Index n, Number* x, Bool new_x,
            Index m, Number* g, UserDataPtr user_data)
{
  Index i;
  Index j;
  
  assert(n == nlp_param->nvars);
  assert(m == nlp_param->mcons);
  
  /** compute the constraint function value for the linear part */
  for( j=0; j<nlp_param->m_LP; ++j)
  {
      g[j] = 0 ;
      for ( i=0; i<nlp_param->nnonz[j]; ++i)
      {
	      g[j] += x[nlp_param->jCols[j][i]]*nlp_param->values[j][i];
      }
  }

  printf(" There are %d  constraints, from 0 to %d is linear, others are nonlinear", m, nlp_param->m_LP);
  printf("The gradient of the objective function value at\n");
  for( i=0; i<n; ++i )
  {
     printf("x[%d] = %f \n", i, x[i]);
  }

  for( j=0; j<nlp_param->m_LP; ++j)
  {
     printf("g[%d] = %f \n", j, g[j]);
  }

  /** compute the constriant function value for the polynomial part */
  for( j=nlp_param->m_LP; j<m; ++j)
  {
      g[j] = computePolynomialValue(nlp_param->polynomials[j-nlp_param->m_LP],x);
  }
  
  printf("now the polynomial constraints from %d to %d \n", nlp_param->m_LP, m-1);
  for ( j=nlp_param->m_LP; j<m; ++j)
  {
      printf("g[%d] = %f \n", j, g[j]);
  }

     
  /** the following code is the example they gave to handle problem  HS_071
   *   assert(n == 4);
   *   assert(m == 2);
   *   g[0] = x[0] * x[1] * x[2] * x[3];
   *   g[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];
   */

  return TRUE;
}

/** given an index tell whether the variable nlp.x[index] is belong to the given monomial */
Bool isVarInMonomial(MonomialIpopt * monomial, Index varIndex)
{
   Index i;
   for(i=0; i<monomial->nvars; ++i)
   {
      if( varIndex == monomial->indicies[i])
         return TRUE;
   }
   return FALSE;
}

/** given an index tell which variable for 0 to monomial->nvars-1 is nlp.x[index], if none of them is nlp.x[index]
    return -1 */
Index positionVarInMomonial(MonomialIpopt * monomial, Index varIndex)
{
	Index i;
	i = -1;
	for(i=0; i<monomial->nvars; ++i)
	{
	   if( varIndex == monomial->indicies[i])
		   return i;
        }
	return -1;
}
/** given an index that means the variable we want to check is nlp.x[index],
    tell whether the variable corresponding to this index is belong to the given polynomial */
Bool isVarInPolynomial(PolynomialIpopt * polynomial, Index varIndex)
{
   Index i;
   for(i=0; i<polynomial->nMonomials; ++i)
   {
      if (isVarInMonomial(polynomial->monomials[i],varIndex) == TRUE )
         return TRUE;
   }
   return FALSE;
}

/** given a monomial, a index and a point compute the gradient value of this index at this point */
Number computeMonomialGradientElement(MonomialIpopt * monomial, Index varIndex, Number *x)
{
	Number val;
	Index i;
	Index j;
	i =  positionVarInMomonial( monomial,varIndex );
	if( i >= 0)
	{
           if( monomial->power[i] == 1 )
           {
              val = monomial->coefficient;
           }
           else
           {
              val = monomial->coefficient*(monomial->power[i])*pow(x[monomial->indicies[i]],monomial->power[i]-1);
           }
           for( j=0; j<monomial->nvars && j != i ; ++j)
           {
              val = val * pow(x[monomial->indicies[j]],monomial->power[j]);
           }

	}
	else
	{ 
		val = 0;
	}
	return val;

}

/** given a monomial, a index and a point compute the gradient value of this index at this point */
Number computePolynomialGradientElement(PolynomialIpopt * polynomial, Index varIndex, Number *x)
{
	Number val;
	Index i;
	val = 0;
	for(i=0;i<polynomial->nMonomials;++i)
	{
		val += computeMonomialGradientElement(polynomial->monomials[i],varIndex, x);
	}

	return val;
}

Bool eval_jac_g(Index n, Number *x, Bool new_x,
                Index m, Index nele_jac,
                Index *iRow, Index *jCol, Number *values,
                UserDataPtr user_data)
{
  Index row; /* row counter for loop , i.e. tell you which constriant we are consider now */
  Index col; /* col counter for loop , i.e. tell you which variables we are consider now */
  Index j;  /* j is the coounter for the LP part */
  Index k; /* k is the counter for the nonzero elemnet in the jacobi of constriants */
  Index i;  /* for printing information */

  assert( n == nlp_param->nvars);
  assert( m == nlp_param->mcons);
  
  if (values == NULL) {
    /* return the structure of the jacobian */
    k = 0;
    printf("Print the structure of the jacobian, no values \n");
    for( row=0; row<nlp_param->m_LP; ++row)
    {
	 for(j=0; j<nlp_param->nnonz[row]; ++j)
	 {
	     iRow[k] = row;
	     jCol[k] = nlp_param->jCols[row][j];
             printf("the %d-th nonzeros in jacobian is in cons %d for  vars %d \n", k, row, jCol[k]);
	     k += 1;
	 }
    }
    for( row=nlp_param->m_LP; row<nlp_param->mcons; ++row)
    {
	    for(col=0; col<n; ++col)
	    {
		    if(isVarInPolynomial(nlp_param->polynomials[row-nlp_param->m_LP],col))
		    {
			    iRow[k] = row;
			    jCol[k] = col;
                            printf("the %d-th nonzeros in jacobian is in cons %d for  vars %d \n", k, row, col);
			    k += 1;
		    }
	    }
	}
    /** QUESTION is it necessay for me, the following line is the one I write follow the style of eval_h of HS071 */
    assert(k == nele_jac);
   
    /** the following is a particular example for HS071 
     *     this particular jacobian is dense 
     *     iRow[0] = 0;
     *     jCol[0] = 0;
     *     iRow[1] = 0;
     *     jCol[1] = 1;
     *     iRow[2] = 0;
     *     jCol[2] = 2;
     *     iRow[3] = 0;
     *     jCol[3] = 3;
     *     iRow[4] = 1;
     *     jCol[4] = 0;
     *     iRow[5] = 1;
     *     jCol[5] = 1;
     *     iRow[6] = 1;
     *     jCol[6] = 2;
     *     iRow[7] = 1;
     *     jCol[7] = 3;
     */
  }
  else {
    /* return the values of the jacobian of the constraints */
    k = 0;
    printf("The values of the jacobian of the constraints at point \n");
    for( i=0; i<n; ++i )
    {
       printf("x[%d] = %f \n", i, x[i]);
    }

    
    printf("first in linear constraints \n");
    for( row=0; row<nlp_param->m_LP; ++row)
    {
	    for(j=0;j<nlp_param->nnonz[row];++j)
	    {
		    values[k] = nlp_param->values[row][j];
                    printf("the %d-th nonzeros in jacobian is in cons %d for  vars %d, value is %f \n", k, row, nlp_param->jCols[row][j], values[k]);
		    k += 1;
	    }
    }

    
    printf("Then for the nonlinear constriants \n");
    for( row=nlp_param->m_LP; row<nlp_param->mcons; ++row)
    {
	    for(col=0;col<n;++col)
	    {
		if(isVarInPolynomial(nlp_param->polynomials[row-nlp_param->m_LP],col))
		{
                   values[k] = computePolynomialGradientElement(nlp_param->polynomials[row-nlp_param->m_LP], col, x);
                   printf("the %d-th nonzeros in jacobian is in cons %d for  vars %d, value is %f \n", k, row, col , values[k]);
			k +=1;
		}
	    }
    }
    assert(k == nele_jac);

    /** the following is a particular example for HS071 
    /* values[0] = x[1]*x[2]*x[3]; /* 0,0 */
    /* values[1] = x[0]*x[2]*x[3]; /* 0,1 */
    /* values[2] = x[0]*x[1]*x[3]; /* 0,2 */
    /* values[3] = x[0]*x[1]*x[2]; /* 0,3 */

    /* values[4] = 2*x[0];         /* 1,0 */
    /* values[5] = 2*x[1];         /* 1,1 */
    /* values[6] = 2*x[2];         /* 1,2 */
    /* values[7] = 2*x[3];         /* 1,3 */
  }
  
  printf(" This is the end of eval_jac");
  /**  getchar(); */

  return TRUE;
}

/** given an index tell whether the variable corresponding to this index is belong to the given monomial at least in power two */
Bool isVarInMonomialTwice(MonomialIpopt * monomial, Index varIndex)
{
   Index i;
   for(i=0; i<monomial->nvars; ++i)
   {
      if( varIndex == monomial->indicies[i])
      {
         if(monomial->power[i] >= 2)
         {
            return TRUE;
         }
         else
         {
            return FALSE;
         }
      }
   }
   return FALSE;
}
/** given an index tell whether the variable corresponding to this index is belong to the given monomial at least in power two */
Bool is2VarsInMonomial(MonomialIpopt * monomial, Index varIndex1, Index varIndex2)
{
   Index i;
   Bool is1;
   Bool is2;
   is1 = FALSE;
   is2 = FALSE;

   for(i=0; i<monomial->nvars; ++i)
   {
      if( varIndex1 == monomial->indicies[i] )
      {
         is1 = TRUE;
      }
      if( varIndex2 == monomial->indicies[i] )
      {
         is2 = TRUE;
      }
   }
   if( is1 == TRUE && is2 == TRUE )
      return TRUE;
   else
      return FALSE;
}
/** given an index tell whether the variable corresponding to this index is belong to the given monomial at least in power two */
Bool isVarInPolynomialTwice(PolynomialIpopt * polynomial, Index varIndex)
{
   Index i;
   for(i=0; i<polynomial->nMonomials; ++i)
   {
      if(isVarInMonomialTwice(polynomial->monomials[i], varIndex) == TRUE )
      {
         return TRUE;
      }
   }
   return FALSE;
}

/** given an index tell whether the variable corresponding to this index is belong to the given monomial at least in power two */
Bool is2VarsInPolynomial(PolynomialIpopt * polynomial, Index varIndex1, Index varIndex2)
{
   Index i;
   for(i=0; i<polynomial->nMonomials; ++i)
   {
      if(is2VarsInMonomial(polynomial->monomials[i], varIndex1, varIndex2) == TRUE )
      {
         return TRUE;
      }
   }
   return FALSE;
}

/** this is a local method to compute the value of a given element in the hession of a monomial at a given point, but this element should not be in diagnonal*/
Number computeHessMonomial(MonomialIpopt * monomial, Index varIndex1, Index varIndex2, Number* x)
{
   Number val;
   Index i;
   Index j;
   Index k;
   i =  positionVarInMomonial(monomial, varIndex1);
   j =  positionVarInMomonial(monomial, varIndex2);
   if( i >= 0 && j >= 0)
   {
      if(monomial->power[i] == 1 && monomial->power[j] == 1 )
      {
         val = monomial->coefficient;
      }
      else if(monomial->power[i] == 1 )
      {
         val = monomial->coefficient*(monomial->power[j])*pow(x[monomial->indicies[j]],monomial->power[j]-1);
      }
      else if(monomial->power[j] == 1 )
      {
         val = monomial->coefficient*(monomial->power[i])*pow(x[monomial->indicies[i]],monomial->power[i]-1);
      }
      else
      {
         val = monomial->coefficient*(monomial->power[i])*pow(x[monomial->indicies[i]],monomial->power[i]-1)
                                 *(monomial->power[j])*pow(x[monomial->indicies[j]],monomial->power[j]-1);
      }
      for( k=0; k<monomial->nvars && k != i && k != j ; ++k)
      {
         val = val * pow(x[monomial->indicies[k]],monomial->power[k]);
      }

   }
   else
   {
      val = 0;
   }

   return val;
}

/** this is a local method to compute the value of a given element in the hession of a monomial at a given point, but this element should be in diagnonal*/
Number computeHessMonomialDiag(MonomialIpopt * monomial, Index varIndex, Number *x)
{
   Number val;
   Index i;
   Index k;
   i =  positionVarInMomonial(monomial, varIndex);
   if( i >= 0 && monomial->power[i] >= 2)
   {
      if(monomial->power[i] == 2 )
      {
         val = 2*monomial->coefficient;
      }
      else
      {
         val = monomial->coefficient*(monomial->power[i])*(monomial->power[i]-1)*pow(x[monomial->indicies[i]],monomial->power[i]-2);
      }
      for( k=0; k<monomial->nvars && k != i ; ++k)
      {
         val = val * pow(x[monomial->indicies[k]],monomial->power[k]);
      }

   }
   else
   {
      val = 0;
   }
   return val;
}

/** this is a local method to compute the value of a given element in the hession of a polynomial at a given point, but this element should not be in diagnonal*/
Number computeHessPolynomial(PolynomialIpopt * polynomial, Index varIndex1, Index varIndex2, Number * x)
{
   Number val;
   Index i;
   val = 0;
   for(i=0; i<polynomial->nMonomials; ++i)
   {
      val += computeHessMonomial(polynomial->monomials[i], varIndex1, varIndex2, x);
   }

   return val;
}

/** this is a local method to compute the value of a given element in the hession of a monomial at a given point, but this element should be in diagnonal*/
Number computeHessPolynomialDiag(PolynomialIpopt * polynomial, Index varIndex, Number * x)
{
   Number val;
   Index i;
   val = 0;
   for(i=0; i<polynomial->nMonomials; ++i)
   {
      val += computeHessMonomialDiag(polynomial->monomials[i], varIndex, x);
   }

   return val;
}

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
            Index m, Number *lambda, Bool new_lambda,
            Index nele_hess, Index *iRow, Index *jCol,
            Number *values, UserDataPtr user_data)
{
  Index idx = 0; /* nonzero element counter */
  Index row = 0; /* row counter for loop */
  Index col = 0; /* col counter for loop */
  Index i;       /* i counter for polynomial */
  Index j;       /* j index for printing */
  Bool isAdd = FALSE;    /* to see whether record this element */
  
  assert( n == nlp_param->nvars);
  assert( m == nlp_param->mcons);
  
  if (values == NULL)
  {
     /* return the structure. This is a symmetric matrix, fill the lower left
      * triangle only. */
     idx = 0;
     /* this part is for the strict lower left triangle */
     printf("give the structure of hessian, give the strict lower left triangle firstly \n");
     for(row = 0; row < n; ++row)
     {
        for(col = 0; col < row; ++col)
        {
           for( i=0; i< nlp_param->m_NLP; ++i)
           {
              if(is2VarsInPolynomial(nlp_param->polynomials[i], row, col) == TRUE )
              {
                 iRow[idx] = row;
                 jCol[idx] = col;
                 printf("The %d-th nonzeros in hessian is x[%d][%d] \n ", idx, row, col);
                 idx++;
              }
              break;
           }
        }
     }
     /* this part is for the diagonal */
     for(row = 0; row < n; ++row)
     {
        for( i=0; i< nlp_param->m_NLP; ++i)
        {
           if(isVarInPolynomialTwice(nlp_param->polynomials[i], row) == TRUE )
           {
              iRow[idx] = row;
              jCol[idx] = row;
              printf("The %d-th nonzeros in hessian is x[%d][%d] \n ", idx, row, row);
              idx++;
           }
           break;
        }
     }
     assert(idx == nele_hess);
     
    /* the folling is the example code for hs071
    /* the hessian for this problem is actually dense 
    idx=0;
    for (row = 0; row < 4; row++) {
      for (col = 0; col <= row; col++) {
        iRow[idx] = row;
        jCol[idx] = col;
        idx++;
      }
    }

    assert(idx == nele_hess); */
  }
  else
  {
     /* return the values. This is a symmetric matrix, fill the lower left
      * triangle only */

     /* for our current problem the hessian of linear obj is zero matrix, the hessian of linear constraint is zero matrix
      * the only part we need to care about is the polynomial constriants part */

     printf("This is the values of hessian at point \n ");
     for( j=0; j<n; ++j )
     {
        printf("x[%d] = %f \n", j, x[j]);

     }


       
     idx = 0;

     /* this part is for the strict lower left triangle */
     printf("give the values of hessian, give the strict lower left triangle firstly \n");
     
     for(row = 0; row < n; ++row)
     {
        for(col = 0; col < row; ++col)
        {
           if( idx < nele_hess)
           {
              isAdd = FALSE;
              values[idx] = 0;
              for( i=0; i< nlp_param->m_NLP; ++i)
              {
                 if(is2VarsInPolynomial(nlp_param->polynomials[i], row, col) == TRUE )
                 {
                    isAdd = TRUE;
                    values[idx] += lambda[i+nlp_param->m_LP]*computeHessPolynomial(nlp_param->polynomials[i], row, col, x);
                 }
              }           if(isAdd == TRUE)
              {
                 printf("The %d-th nonzeros in hessian is x[%d]x[%d], its value is %f \n ", idx, row, col, values[idx]);
                 idx++;
              }
           }
        }
     }
     
     
     /* this part is for the diagonal */
     for(row = 0; row < n; ++row)
     {
        if (idx < nele_hess)
        {
        isAdd = FALSE;
        values[idx] = 0;
        for( i=0; i< nlp_param->m_NLP; ++i)
        {
           if(isVarInPolynomialTwice(nlp_param->polynomials[i], row) == TRUE )
           {
              isAdd = TRUE;
              values[idx] += lambda[i+nlp_param->m_LP]*computeHessPolynomialDiag(nlp_param->polynomials[i], row, x);
           }
        }
        if(isAdd == TRUE)
        {
           printf("The %d-th nonzeros in hessian is x[%d]x[%d], its value is %f \n ", idx, row, row, values[idx]);
           idx++;
        }
        }
     }

     
     assert(idx == nele_hess);


    /** the following are the codes for hs071 example  
    /* fill the objective portion */
    /* values[0] = obj_factor * (2*x[3]);               /* 0,0 */

    /* values[1] = obj_factor * (x[3]);                 /* 1,0 */
    /* values[2] = 0;                                   /* 1,1 */

    /* values[3] = obj_factor * (x[3]);                 /* 2,0 */
    /* values[4] = 0;                                   /* 2,1 */
    /* values[5] = 0;                                   /* 2,2 */

    /* values[6] = obj_factor * (2*x[0] + x[1] + x[2]); /* 3,0 */
    /* values[7] = obj_factor * (x[0]);                 /* 3,1 */
    /* values[8] = obj_factor * (x[0]);                 /* 3,2 */
    /* values[9] = 0;                                   /* 3,3 */


    /* add the portion for the first constraint */
    /* values[1] += lambda[0] * (x[2] * x[3]);          /* 1,0 */

    /* values[3] += lambda[0] * (x[1] * x[3]);          /* 2,0 */
    /* values[4] += lambda[0] * (x[0] * x[3]);          /* 2,1 */

    /* values[6] += lambda[0] * (x[1] * x[2]);          /* 3,0 */
    /* values[7] += lambda[0] * (x[0] * x[2]);          /* 3,1 */
    /* values[8] += lambda[0] * (x[0] * x[1]);          /* 3,2 */

    /* add the portion for the second constraint */
    /* values[0] += lambda[1] * 2;                      /* 0,0 */

    /* values[2] += lambda[1] * 2;                      /* 1,1 */

    /* values[5] += lambda[1] * 2;                      /* 2,2 */

    /* values[9] += lambda[1] * 2;                      /* 3,3 */
  }
  printf("This is the end of eval_h \n ");
  /**  getchar();*/
  return TRUE;
}
