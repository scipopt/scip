/**
 * @file model.h
 * @brief Basic classes to describe a MIP model
 *
 * @author Domenico Salvagnin
 */

#ifndef MODEL_H
#define MODEL_H

#include "gmputils.h"
#include <string>
#include <vector>
#include <iosfwd>
#include <map>

/**
 * @brief Class representing a problem variable.
 * Holds global variable information (such as bounds and objective coefficent)
 * and a solution value (default to zero)
 */

class Var
{
public:
   /**
    * @enum VarType Var domain type
    */
   enum VarType
   {
      BINARY,
      INTEGER,
      CONTINUOUS
   };

   /**
    * Constructor
    * @param _name name of the variable
    * @param _type type of domain
    * @param _lb global lower bound
    * @param _ub global upper bound
    * @param _obj objective coefficent
    */
   Var(const char* _name, VarType _type, const Rational& _lb, const Rational& _ub, const Rational& _obj);

   /**
    * Check if the current value of the variable is within its bounds.
    * @param boundTolerance absolute tolerance for bounds check
    * @return true if the variable satisfies its bounds, false otherwise
    */
   bool checkBounds(const Rational& boundTolerance) const;

   /**
    * Check if the current value of the variable satisfies the integrality requirement.
    * This check is always true if the variable is continuous.
    * @param intTolerance absolute tolerance for integrality check
    * @return true if the value is integral or the variable is continuous, false otherwise
    */
   bool checkIntegrality(const Rational& intTolerance) const;

   /**
    * Calculate the bounds violation of the variable (0 if not violated).
    * Return value is in @param boundViol
    */
   void boundsViolation(Rational& boundViol) const;

   /**
    * Calculate the integrality violation of the variable (0 if not violated).
    * Return value is in @param intViol
    */
   void integralityViolation(Rational& intViol) const;

   /**
    * Print a description of the variable to an output stream
    */
   void print(std::ostream& out) const;

   /** @name public data */
   // \{
   std::string name; //< name of the variable
   VarType type; //< type of domain
   Rational lb; //< global lower bound
   Rational ub; //< global upper bound
   Rational objCoef; //< objective coefficent
   Rational value; //< solution value as reported by the solver (default is zero)
   // \}
};


/**
 * @brief Base class representing a generic constraint.
 */

class Constraint
{
public:
	/**
    * Constructor
    * @param _name name of the constraint
    * @param _redundant this flag is true if the constraint is redundant for the problems (e.g. is a cutting plane)
    */
   Constraint(const char* _name, bool _redundant = false);

	/**
	 * Virtual Destructor
	 */

	virtual ~Constraint() {}

	/**
	 * Check if the constraint is satisfied by the current values of its variables.
	 * @param tolerance tolerance for checking feasibility
	 * @return true if the constraint is satisfied, false otherwise.
	 */
	virtual bool check(const Rational& tolerance) const =0;

	/**
    * Calculate the constraint violation
    * Return value is in @param viol
    */
   virtual void violation(Rational& viol) const =0;

	/**
	 * Print a description of the constraint to an output stream
	 */
	virtual void print(std::ostream& out) const =0;

	/** @name public data */
	// \{
	std::string name; //< name of the constraint
	std::string type; //< constraint type
	bool redundant; //< is the constraint redundant in the model?
	// \}
};


/**
 * @brief Class representing a linear constraint.
 * Constraints are always stored as ranged constraints, i.e.
 * as lhs <= a^T x <= rhs, and in sparse form.
 */

class LinearConstraint : public Constraint
{
public:
   /**
    * @enum LinearType Linear Constraint type (i.e. sense)
    */
   enum LinearType
   {
      LESS_THAN,
      GREATER_THAN,
      EQUAL,
      RANGED
   };

   /**
    * Constructor
    * @param _name name of the constraint
    * @param _lintype type of constraint
    * @param _lhs left hand side of the constraint
    * @param _rhs right hand side of the constraint
    * @param _redundant this flag is true if the constraint is redundant for the problems (e.g. is a cutting plane)
    */
   LinearConstraint(const char* _name, LinearType _lintype, const Rational& _lhs, const Rational& _rhs, bool _redundant = false);

   /**
    * Add a variable to the constraint. Does NOT check for duplicates.
    * @param v pointer to variable object
    * @param c coefficient of the variable in the constraint
    */
   void push(Var* v, const Rational& c);

   /**
    * Check if the constraint is satisfied by the current values of its variables.
    * The check is mostly a relative tolerance one.
    *
    * Let's say we have a constraint of the form ax <= rhs, a solution x* and a tolerance eps.
    * Then the constraint is considered as violated if:
    * ax^* - rhs > eps * max {(ax^*)_+, (ax^*)_-, |rhs|, 1}
    * Note that we split the activity into positive and negative part to avoid cancellation effects,
    * and we still put 1 into the arguments of max in order to get back to an absolute tolerance
    * in case of small values. The same reasoning applies for lhs, of course.
    *
    * @param tolerance tolerance for checking feasibility
    * @return true if the constraint is satisfied, false otherwise.
    */
   bool check(const Rational& tolerance) const;

   /**
    * Calculate the constraint violation
    * Return value is in @param viol
    */
   void violation(Rational& viol) const;

   /**
    * Print a description of the constraint to an output stream
    */
   void print(std::ostream& out) const;

   /** @name public data */
   // \{
   LinearType lintype; //< constraint type
   std::vector<Var*> vars; //< list of variables
   std::vector<Rational> coefs; //< list of corresponding coefficients
   Rational lhs; //< left hand side
   Rational rhs; //< right hand side
   // \}
};


/**
 * @brief Class representing a SOS constraint.
 * The order of the variables is given by the order in which they are added to the constraints
 * No weight information is used for now!
 */

class SOSConstraint : public Constraint
{
public:
   /**
    * @enum SOSType SOS type (i.e. 1 or 2)
    */
   enum SOSType
   {
      TYPE_1,
      TYPE_2
   };

   /**
    * Constructor
    * @param _name name of the constraint
    * @param _sostype type of constraint
    * @param _redundant this flag is true if the constraint is redundant for the problems (e.g. is a cutting plane)
    */
   SOSConstraint(const char* _name, SOSType _sostype, bool _redundant = false);

   /**
    * Add a variable to the SOS constraint. Does NOT check for duplicates.
    * @param v pointer to variable object
    */
   void push(Var* v);

   /**
    * Check if the constraint is satisfied by the current values of its variables.
    * @param tolerance tolerance for checking feasibility
    * @return true if the constraint is satisfied, false otherwise.
    */
   bool check(const Rational& tolerance) const;

   /**
    * Calculate the constraint violation
    * Since this is a "combinatorial" constraint, we always return 0 for now...FIXME
    * Return value is in @param viol
    */
   void violation(Rational& viol) const;

   /**
    * Print a description of the constraint to an output stream
    */
   void print(std::ostream& out) const;

   /** @name public data */
   // \{
   SOSType sostype; //< constraint type
   std::vector<Var*> vars; //< list of variables
   // \}
protected:
   // helpers
   bool checkType1(const Rational& tolerance) const;
   bool checkType2(const Rational& tolerance) const;
};


/**
 * @brief Class representing a indicator constraint.
 * An indicator constraint is a constraint of the form y = 0/1 -> ax <=/=/>= b
 * where y is a binary variable and the consequence is an arbitrary linear constraint.
 */

class IndicatorConstraint : public Constraint
{
public:
   /**
    * Constructor
    * @param _name name of the constraint
    * @param _ifvar binary variable premise of the implication
    * @param _ifvalue value that ifvar must take to triggers the implication
    * @param _thencons consequence of the implication (an arbitrary constraint)
    * @param _redundant this flag is true if the constraint is redundant for the problems (e.g. is a cutting plane)
    */
   IndicatorConstraint(const char* _name, Var* _ifvar, bool _ifvalue, Constraint* _thencons, bool _redundant = false);

   ~IndicatorConstraint();

   /**
    * Check if the constraint is satisfied by the current values of its variables.
    * @param tolerance tolerance for checking feasibility
    * @return true if the constraint is satisfied, false otherwise.
    */
   bool check(const Rational& tolerance) const;

   /**
    * Calculate the constraint violation
    * Return value is in @param viol
    */
   void violation(Rational& viol) const;

   /**
    * Print a description of the constraint to an output stream
    */
   void print(std::ostream& out) const;

   /** @name public data */
   // \{
   Var* ifvar;
   bool ifvalue;
   Constraint* thencons;
   // \}
};


/**
 * @brief Class representing a MIP problem.
 * Holds the list of variables and constraints of the model.
 * The model owns all variable and constraint objects, and frees
 * them on destruction.
 */

class Model
{
public:
   /**
    * @enum ObjSense Optimization direction (i.e. minimization or maximization)
    */
   enum ObjSense
   {
      MINIMIZE,
      MAXIMIZE
   };

   /** Constructor */
   Model();
   /** Destructor */
   ~Model();

   /**
    * Get a variable by name
    * @return a pointer to the variable with name @param name if found, NULL otherwise
    */
   Var* getVar(const char* name) const;

   /**
    * Get a constraint by name
    * @return a pointer to the constraint with name @param name if found, NULL otherwise
    */
   Constraint* getCons(const char* name) const;

   /**
    * Add a variable to the model.
    * If a variable with the same name exists it is replaced by the new one.
    * @param var variable to add
    */
   void pushVar(Var* var);

   /**
    * Add a constraint to the model
    * If a constraint with the same name exists it is replaced by the new one.
    * @param cons constraint to add
    */
   void pushCons(Constraint* cons);

  /**
   * Remove a constraint (if it exists) from the model
   * @param name name of the constraint to remove
   */
  void removeCons(const char* name);

   /**
    * Get number of variables
    */
   unsigned int numVars() const;

   /**
    * Get number of constraints
    */
   unsigned int numConss() const;

   /**
    * Read solution values for the variables from file.
    * The format of the solution is very simple. In each line
    * we expect a pair of <variable name, value> separated by spaces. \n
    *
    * Example: \n
    * x1     3.23232 \n
    * x2    34.00000 \n
    * ... \n
    *
    * @param filename path to the file with the solution values
    * @return true if successful, false otherwise
    */
   bool readSol(const char* filename);

   /**
    * Check if the model is satisfied by the current values of the variables.
    * Checks both domains and linear constraints.
    * @param intTolerance tolerance for integrality check
    * @param linearTolerance tolerance for constraint (and bound) checks and
    * for comparing the real objective value with the one given by the solver (if any)
    * @param intFeasible stores if the solution satisfies the integrality requirements
    * @param linearFeasible stores if the solution satisfies the linear constraints
    * @param correctObj stores if the objective value computed by the solver is correct
    */
   void check(
      const Rational& intTolerance,
      const Rational& linearTolerance,
      bool& intFeasible,
      bool& linearFeasible,
      bool& correctObj) const;

   /**
    * Calculate the maximum integrality, linear and objective violations with the
    * current values of the variables
    */
   void maxViolations(
      Rational& intViol,
      Rational& linearViol,
      Rational& objViol) const;

   /**
    * Check w.r.t. exact solver information
    */
   void checkWrtExact(
      const std::string& solverStatus,
      const std::string& exactStatus,
      const Rational& exactObjVal,
      const Rational& linearTolerance,
      bool& feasibility,
      bool& objective) const;

   /**
    * Print the model to an output stream (for debugging)
    */
   void print(std::ostream& out) const;

   /**
    * Print the current solution to an output stream (for debugging)
    */
   void printSol(std::ostream& out) const;

   /** @name public data */
   // \{
   std::string modelName; //< name of the model
   std::string objName; //< name of the objective function
   ObjSense objSense; //< optimization direction
   bool hasObjectiveValue; //< do we have the objective value of the solution as reported by the solver?
   Rational objectiveValue; //< objective value of the solution as reported by the solver
   // \}
protected:
   std::map<std::string, Var*> vars; //< set of variables
   std::map<std::string, Constraint*> conss; //< set of constraints
};

#endif /* MODEL_H */
