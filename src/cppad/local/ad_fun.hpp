/* $Id: ad_fun.hpp 2178 2011-10-30 06:52:58Z bradbell $ */
# ifndef CPPAD_AD_FUN_INCLUDED
# define CPPAD_AD_FUN_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
/*
$begin ADFun$$
$spell 
	xk
	Ind
	bool
	taylor_
	sizeof
	const
	std
	ind_taddr_
	dep_taddr_
$$

$spell
$$

$section ADFun Objects$$

$index ADFun, object$$
$index object, ADFun$$

$head Purpose$$
An AD of $italic Base$$
$xref/glossary/Operation/Sequence/operation sequence/1/$$
is stored in an $code ADFun$$ object by its $xref/FunConstruct/$$.
The $code ADFun$$ object can then be used to calculate function values,
derivative values, and other values related to the corresponding function.

$childtable%
	cppad/local/independent.hpp%
	cppad/local/fun_construct.hpp%
	cppad/local/dependent.hpp%
	cppad/local/abort_recording.hpp%
	omh/seq_property.omh%
	cppad/local/fun_eval.hpp%
	cppad/local/drivers.hpp%
	cppad/local/fun_check.hpp%
	cppad/local/optimize.hpp
%$$

$end
*/

CPPAD_BEGIN_NAMESPACE
/*!
\file ad_fun.hpp
File used to define the ADFun<Base> class.
*/

/*!
Class used to hold function objects

\tparam Base
A function object has a recording of <tt>AD<Base></tt> operations.
It does it calculations using \c Base operations.
*/

template <class Base>
class ADFun {
// ------------------------------------------------------------
// Private member variables
private:
	/// debug checking number of comparision operations that changed
	size_t compare_change_;

	/// number of taylor_ coefficieint per variable (currently stored)
	size_t taylor_per_var_;

	/// number of columns currently allocated for taylor_ array
	size_t taylor_col_dim_;

	/// number of rows (variables) in the recording (play_)
	size_t total_num_var_;

	/// tape address for the independent variables
	CppAD::vector<size_t> ind_taddr_;

	/// tape address and parameter flag for the dependent variables
	CppAD::vector<size_t> dep_taddr_;

	/// which dependent variables are actually parameters
	CppAD::vector<bool>   dep_parameter_;

	/// the operation sequence corresponding to this object
	player<Base> play_;

	/// results of the forward mode calculations
	pod_vector<Base> taylor_;

	/// packed results of the forward mode Jacobian sparsity calculations
	/// (\c for_jac_sparse_pack_.n_set() != 0  implies
	/// for_jac_sparse_set_.n_set() == 0)
	sparse_pack      for_jac_sparse_pack_;

	/// set results of the forward mode Jacobian sparsity calculations
	/// (\c for_jac_sparse_set_.n_set() != 0  implies
	/// for_jac_sparse_pack_.n_set() == 0)
	sparse_set       for_jac_sparse_set_;

// ------------------------------------------------------------
// Private member functions

	/// change the operation sequence corresponding to this object
	template <typename ADvector>
	void Dependent(ADTape<Base> *tape, const ADvector &y);

	// ------------------------------------------------------------
	// vector of bool version of ForSparseJac
	// (see doxygen in for_sparse_jac.hpp)
	template <class VectorSet>
	void ForSparseJacCase(
		bool               set_type  ,
		size_t             q         ,
		const VectorSet&   r         ,  
		VectorSet&         s
	);
	// vector of std::set<size_t> version of ForSparseJac
	// (see doxygen in for_sparse_jac.hpp)
	template <class VectorSet>
	void ForSparseJacCase(
		const std::set<size_t>&  set_type  ,
		size_t                   q         ,
		const VectorSet&         r         ,  
		VectorSet&               s
	);
	// ------------------------------------------------------------
	// vector of bool version of RevSparseJac
	// (see doxygen in rev_sparse_jac.hpp)
	template <class VectorSet>
	void RevSparseJacCase(
		bool               set_type  ,
		size_t             p         ,
		const VectorSet&   s         ,  
		VectorSet&         r
	);
	// vector of std::set<size_t> version of RevSparseJac
	// (see doxygen in rev_sparse_jac.hpp)
	template <class VectorSet>
	void RevSparseJacCase(
		const std::set<size_t>&  set_type  ,
		size_t                   p         ,
		const VectorSet&         s         ,  
		VectorSet&               r
	);
	// ------------------------------------------------------------
	// vector of bool version of RevSparseHes
	// (see doxygen in rev_sparse_hes.hpp)
	template <class VectorSet>
	void RevSparseHesCase(
		bool               set_type  ,
		size_t             q         ,
		const VectorSet&   s         ,  
		VectorSet&         h
	);
	// vector of std::set<size_t> version of RevSparseHes
	// (see doxygen in rev_sparse_hes.hpp)
	template <class VectorSet>
	void RevSparseHesCase(
		const std::set<size_t>&  set_type  ,
		size_t                   q         ,
		const VectorSet&         s         ,  
		VectorSet&               h
	);
	// ------------------------------------------------------------
	// vector of bool version of SparseJacobian
	// (see doxygen in sparse_jacobian.hpp)
	template <class VectorBase, class VectorSet>
	void SparseJacobianCase(
		bool                     set_type    ,
		const VectorBase&        x           ,
		const VectorSet&         p           ,
		VectorBase&              jac
	);
	// vector of std::set<size_t> version of SparseJacobian
	// (see doxygen in sparse_jacobian.hpp)
	template <class VectorBase, class VectorSet>
	void SparseJacobianCase(
		const std::set<size_t>&  set_type    ,
		const VectorBase&        x           ,
		const VectorSet&         p           ,
		VectorBase&              jac
	);
	// ------------------------------------------------------------
	// vector of bool version of SparseHessian
	// (see doxygen in sparse_hessian.hpp)
	template <class VectorBase, class VectorSet>
	void SparseHessianCase(
		bool                     set_type    ,
		const VectorBase&        x           ,
		const VectorBase&        w           ,
		const VectorSet&         p           ,
		VectorBase&              hes
	);
	// vector of std::set<size_t> version of SparseHessian
	// (see doxygen in sparse_hessian.hpp)
	template <class VectorBase, class VectorSet>
	void SparseHessianCase(
		const std::set<size_t>&  set_type    ,
		const VectorBase&        x           ,
		const VectorBase&        w           ,
		const VectorSet&         p           ,
		VectorBase&              hes
	);
// ------------------------------------------------------------
public:
	/// copy constructor
	ADFun(const ADFun& g) 
	: total_num_var_(0)
	{	CppAD::ErrorHandler::Call(
		true,
		__LINE__,
		__FILE__,
		"ADFun(const ADFun& g)",
		"Attempting to use the ADFun<Base> copy constructor.\n"
		"Perhaps you are passing an ADFun<Base> object "
		"by value instead of by reference."
		);
	 }

	/// default constructor
	ADFun(void); 

	// assignment operator
	// (see doxygen in fun_construct.hpp)
	void operator=(const ADFun& f);

	/// sequence constructor
	template <typename ADvector>
	ADFun(const ADvector &x, const ADvector &y);

	/// destructor
	~ADFun(void)
	{ }

	/// deprecated: assign a new operation sequence
	template <typename ADvector>
	void Dependent(const ADvector &y);

	/// assign a new operation sequence
	template <typename ADvector>
	void Dependent(const ADvector &x, const ADvector &y);

	/// forward mode sweep
	template <typename VectorBase>
	VectorBase Forward(
		size_t p, const VectorBase& x, std::ostream& s = std::cout);

	/// reverse mode sweep
	template <typename VectorBase>
	VectorBase Reverse(size_t p, const VectorBase &v);

	// forward mode Jacobian sparsity 
	// (see doxygen documentation in for_sparse_jac.hpp)
	template <typename VectorSet>
	VectorSet ForSparseJac(
		size_t q, const VectorSet &r
	);
	// reverse mode Jacobian sparsity 
	// (see doxygen documentation in rev_sparse_jac.hpp)
	template <typename VectorSet>
	VectorSet RevSparseJac(
		size_t q, const VectorSet &s
	);
	// reverse mode Hessian sparsity 
	// (see doxygen documentation in rev_sparse_hes.hpp)
	template <typename VectorSet>
	VectorSet RevSparseHes(
		size_t q, const VectorSet &s
	);

	/// amount of memeory used for Jacobain sparsity pattern
	size_t size_forward_bool(void) const
	{	return for_jac_sparse_pack_.memory(); }

	/// free memeory used for Jacobain sparsity pattern
	void size_forward_bool(size_t zero) 
	{	CPPAD_ASSERT_KNOWN(
			zero == 0,
			"size_forward_bool: argument not equal to zero"
		);
		for_jac_sparse_pack_.resize(0, 0); 
	}

	/// total number of elements used for Jacobian sparsity pattern
	size_t size_forward_set(void) const
	{	return for_jac_sparse_set_.number_elements(); }

	/// free memeory used for Jacobain sparsity pattern
	void size_forward_set(size_t zero)
	{	CPPAD_ASSERT_KNOWN(
			zero == 0,
			"size_forward_bool: argument not equal to zero"
		);
		for_jac_sparse_set_.resize(0, 0); 
	}

	/// number of operators in the operation sequence
	size_t size_op(void) const
	{	return play_.num_rec_op(); }

	/// number of operator arguments in the operation sequence
	size_t size_op_arg(void) const
	{	return play_.num_rec_op_arg(); }

	/// amount of memory required for the operation sequence
	size_t size_op_seq(void) const
	{	return play_.Memory(); }

	/// number of parameters in the operation sequence
	size_t size_par(void) const
	{	return play_.num_rec_par(); }

	/// number of taylor_ coefficients currently calculated (per variable)
	size_t size_taylor(void) const
	{	return taylor_per_var_; } 

	/// number of characters in the operation sequence
	size_t size_text(void) const
	{	return play_.num_rec_text(); }

	/// number of variables in opertion sequence
	size_t size_var(void) const
	{	return total_num_var_; }

	/// number of VecAD indices in the operation sequence
	size_t size_VecAD(void) const
	{	return play_.num_rec_vecad_ind(); }

	/// set number of coefficients currently allocated (per variable)
	void capacity_taylor(size_t per_var);   

	/// number of independent variables
	size_t Domain(void) const
	{	return ind_taddr_.size(); }

	/// number of dependent variables
	size_t Range(void) const
	{	return dep_taddr_.size(); }

	/// is variable a parameter
	bool Parameter(size_t i)
	{	CPPAD_ASSERT_KNOWN(
			i < dep_taddr_.size(),
			"Argument to Parameter is >= dimension of range space"
		);
		return dep_parameter_[i]; 
	}

# ifndef NDEBUG
	/// in not NDEBUG case, number of comparison operations that change
	size_t CompareChange(void) const
	{	return compare_change_; }
# endif

	/// calculate entire Jacobian
	template <typename VectorBase>
	VectorBase Jacobian(const VectorBase &x); 

	/// calculate Hessian for one component of f
	template <typename VectorBase>
	VectorBase Hessian(const VectorBase &x, const VectorBase &w); 
	template <typename VectorBase>
	VectorBase Hessian(const VectorBase &x, size_t i); 

	/// forward mode calculation of partial w.r.t one domain component
	template <typename VectorBase>
	VectorBase ForOne(
		const VectorBase   &x ,
		size_t              j );

	/// reverse mode calculation of derivative of one range component
	template <typename VectorBase>
	VectorBase RevOne(
		const VectorBase   &x ,
		size_t              i );

	/// forward mode calculation of a subset of second order partials
	template <typename VectorBase, typename VectorSize_t>
	VectorBase ForTwo(
		const VectorBase   &x ,
		const VectorSize_t &J ,
		const VectorSize_t &K );

	/// reverse mode calculation of a subset of second order partials
	template <typename VectorBase, typename VectorSize_t>
	VectorBase RevTwo(
		const VectorBase   &x ,
		const VectorSize_t &I ,
		const VectorSize_t &J );

	/// calculate sparse Jacobians 
	template <typename VectorBase>
	VectorBase SparseJacobian(const VectorBase &x); 

	/// calculate sparse Jacobians 
	template <typename VectorBase, typename VectorSet>
	VectorBase SparseJacobian(const VectorBase &x, const VectorSet &p); 

	/// calculate sparse Hessians 
	template <typename VectorBase>
	VectorBase SparseHessian(const VectorBase &x, const VectorBase &w); 

	/// calculate sparse Hessians 
	template <typename VectorBase, typename VectorBool>
	VectorBase SparseHessian(
		const VectorBase &x, const VectorBase &w, const VectorBool &p
	); 

	// Optimize the tape
	// (see doxygen documentation in optimize.hpp)
	void optimize(void);
	// ------------------- Deprecated -----------------------------

	/// deprecated: number of variables in opertion sequence
	size_t Size(void) const
	{	return total_num_var_; }

	/// deprecated: # taylor_ coefficients currently stored (per variable)
	size_t Order(void) const
	{	return taylor_per_var_ - 1; }

	/// Deprecated: amount of memory for this object 
	/// Note that an approximation is used for the std::set<size_t> memory
	size_t Memory(void) const
	{	size_t pervar  = taylor_col_dim_ * sizeof(Base)
		+ for_jac_sparse_pack_.memory()
		+ 3 * sizeof(size_t) * for_jac_sparse_set_.number_elements();
		size_t total   = total_num_var_ * pervar + play_.Memory();
		return total;
	}

	/// deprecated: # taylor_ coefficients stored (per variable)
	size_t taylor_size(void) const
	{	return taylor_per_var_; } 

	/// deprecated: Does this AD operation sequence use 
	//VecAD<Base>::reference operands
	bool use_VecAD(void) const
	{	return play_.num_rec_vecad_ind() > 0; }
};
// ---------------------------------------------------------------------------

CPPAD_END_NAMESPACE

// non-user interfaces
# include <cppad/local/forward0sweep.hpp>
# include <cppad/local/forward_sweep.hpp>
# include <cppad/local/reverse_sweep.hpp>
# include <cppad/local/for_jac_sweep.hpp>
# include <cppad/local/rev_jac_sweep.hpp>
# include <cppad/local/rev_hes_sweep.hpp>

// user interfaces
# include <cppad/local/parallel_ad.hpp>
# include <cppad/local/independent.hpp>
# include <cppad/local/dependent.hpp>
# include <cppad/local/fun_construct.hpp>
# include <cppad/local/abort_recording.hpp>
# include <cppad/local/fun_eval.hpp>
# include <cppad/local/drivers.hpp>
# include <cppad/local/fun_check.hpp>
# include <cppad/local/omp_max_thread.hpp>
# include <cppad/local/optimize.hpp> 

# endif
