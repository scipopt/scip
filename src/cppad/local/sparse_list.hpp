/* $Id: sparse_list.hpp 2910 2013-10-07 13:27:58Z bradbell $ */
# ifndef CPPAD_SPARSE_LIST_INCLUDED
# define CPPAD_SPARSE_LIST_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-12 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Eclipse Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
# include <cppad/local/define.hpp>
# include <list>

namespace CppAD { // BEGIN_CPPAD_NAMESPACE
/*!
\defgroup sparse_list_hpp sparse_list.hpp
\{
\file sparse_list.hpp
Vector of sets of positive integers stored as singly linked lists.
*/

/*!
Vector of sets of positive integers, each set stored as a singly
linked list.
*/
class sparse_list {
private:
	/// type used for each set in the vector of sets
	/// (note that next is index in data_ for next element in this set).
	struct Element { size_t value; size_t next; };
	/// Number of sets that we are representing
	/// (set by constructor and resize).
	size_t n_set_;
	/// Possible elements in each set are 0, 1, ..., end_ - 1;
	size_t end_;
	/// The data for all the elements in all the sets
	CppAD::pod_vector<Element> data_;
	/// next element
	Element next_element_;
public:
	// -----------------------------------------------------------------
	/*! Default constructor (no sets)
	*/
	sparse_list(void) : 
	n_set_(0)                     , 
	end_(0) 
	{	next_element_.value = end_; }
	// -----------------------------------------------------------------
	/*! Make use of copy constructor an error
 
	\param v
	vector that we are attempting to make a copy of.
 	*/
	sparse_list(const sparse_list& v)
	{	// Error: 
		// Probably a sparse_list argument has been passed by value
		CPPAD_ASSERT_UNKNOWN(false); 
	}
	// -----------------------------------------------------------------
	/*! Destructor 
	*/
	~sparse_list(void)
	{ } 
	// -----------------------------------------------------------------
	/*! Change number of sets, end marker, and initialize all sets as empty

	If \c n_set_in is zero, any memory currently allocated for this object 
	is freed. Otherwise, new memory may be allocated for the sets (if needed).

	\param n_set_in
	is the number of sets in this vector of sets.

	\param end_in
	is the maximum element plus one (the minimum element is 0).
	*/
	void resize(size_t n_set_in, size_t end_in) 
	{	n_set_                 = n_set_in;
		end_                   = end_in;
		next_element_.value    = end_in;
		if( n_set_in == 0 )
		{	// free all memory connected with this object
			data_.free();
			return;
		}
		// now start a new vector with empty sets
		data_.erase();
		data_.extend(n_set_);
		for(size_t i = 0; i < n_set_; i++)
			data_[i].value = end_;
	}
	// -----------------------------------------------------------------
	/*! Add one element to a set.

	\param index
	is the index for this set in the vector of sets.

	\param element
	is the element we are adding to the set.

	\par Checked Assertions
	\li index    < n_set_
	\li element  < end_
	*/
	void add_element(size_t index, size_t element)
	{	CPPAD_ASSERT_UNKNOWN( index   < n_set_ );
		CPPAD_ASSERT_UNKNOWN( element < end_ );

		size_t value = data_[index].value;

		// case of inserting at beginning
		if( element < value )
		{	size_t insert       = data_.extend(1);
			data_[insert]       = data_[index];
			data_[index].value  = element;
			data_[index].next   = insert;
			return;
		}

		// search list for place to insert 
		size_t  previous = index;
		size_t  current  = data_[previous].next;
		value            = data_[current].value;
		while( value < element )
		{	previous = current;
			current = data_[previous].next;
			value   = data_[current].value; 
		}
		if( element != value )
		{	CPPAD_ASSERT_UNKNOWN( element < value ); 
			size_t insert         = data_.extend(1);
			//
			data_[insert].next    = data_[previous].next;
			data_[previous].next  = insert;
			data_[insert].value   = element;
		}
	}
	// -----------------------------------------------------------------
	/*! Is an element of a set.

	\param index
	is the index for this set in the vector of sets.

	\param element
	is the element we are adding to the set.

	\par Checked Assertions
	\li index    < n_set_
	\li element  < end_
	*/
	bool is_element(size_t index, size_t element)
	{	CPPAD_ASSERT_UNKNOWN( index   < n_set_ );
		CPPAD_ASSERT_UNKNOWN( element < end_ );

		size_t i     = index;
		size_t value = data_[i].value;
		while( value < element )
		{	i     = data_[i].next;
			value = data_[i].value;
		}

		return element == value;
	}
	// -----------------------------------------------------------------
	/*! Begin retrieving elements from one of the sets.
	
	\param index
	is the index for the set that is going to be retrieved.
	The elements of the set are retrieved in increasing order.

	\par Checked Assertions
	\li index  < n_set_

	*/
	void begin(size_t index)
	{	// initialize element to search for in this set
		CPPAD_ASSERT_UNKNOWN( index < n_set_ );
		next_element_  = data_[index]; 

		return;
	}
	// -----------------------------------------------------------------
	/*! Get the next element from the current retrieval set.
	
	\return
	is the next element in the set with index
	specified by the previous call to \c begin.
	If no such element exists, \c this->end() is returned.

	\par Assumption
	There is no call to \c add_element since the previvious \c begin
	*/
	size_t next_element(void)
	{	size_t element = next_element_.value;
		if( element != end_ )
			next_element_ = data_[next_element_.next];

		return element;
	}
	// -----------------------------------------------------------------
	/*! Assign the empty set to one of the sets.

	\param target
	is the index of the set we are setting to the empty set.

	\par Checked Assertions
	\li target < n_set_
	*/
	void clear(size_t target)
	{	CPPAD_ASSERT_UNKNOWN( target < n_set_ );
		data_[target].value = end_;
	}
	// -----------------------------------------------------------------
	/*! Assign one set equal to another set.

	\param this_target
	is the index (in this \c sparse_list object) of the set being assinged.

	\param other_source
	is the index (in the other \c sparse_list object) of the 
	that we are using as the value to assign to the target set.

	\param other
	is the other \c sparse_list object (which may be the same as this
	\c sparse_list object).

	\par Checked Assertions
	\li this_target  < n_set_
	\li other_index  < other.n_set_
	*/
	void assignment(
		size_t               this_target  , 
		size_t               other_source , 
		const sparse_list&   other        )
	{	CPPAD_ASSERT_UNKNOWN( this_target  <   n_set_        );
		CPPAD_ASSERT_UNKNOWN( other_source <   other.n_set_  );
		CPPAD_ASSERT_UNKNOWN( end_        == other.end()   );

		size_t this_index        = this_target;
		size_t other_index       = other_source;
		size_t value             = other.data_[other_index].value;
		while( value != end_ )
		{	size_t next             = data_.extend(1);
			data_[this_index].value = value;
			data_[this_index].next  = next;
			this_index              = next;
			other_index             = other.data_[other_index].next;
			value                   = other.data_[other_index].value;
		}
		data_[this_index].value = end_;
	}
	// -----------------------------------------------------------------
	/*! Assign a set equal to the union of two other sets.

	\param this_target
	is the index (in this \c sparse_list object) of the set being assinged.

	\param this_left
	is the index (in this \c sparse_list object) of the 
	left operand for the union operation.
	It is OK for \a this_target and \a this_left to be the same value.

	\param other_right
	is the index (in the other \c sparse_list object) of the 
	right operand for the union operation.
	It is OK for \a this_target and \a other_right to be the same value.

	\param other
	is the other \c sparse_list object (which may be the same as this
	\c sparse_list object).

	\par Checked Assertions
	\li this_target <  n_set_
	\li this_left   <  n_set_
	\li other_right <  other.n_set_
	*/
	void binary_union(
		size_t                  this_target  , 
		size_t                  this_left    , 
		size_t                  other_right  , 
		const sparse_list&      other        )
	{	
		CPPAD_ASSERT_UNKNOWN( this_target < n_set_         );
		CPPAD_ASSERT_UNKNOWN( this_left   < n_set_         );
		CPPAD_ASSERT_UNKNOWN( other_right < other.n_set_   );
		CPPAD_ASSERT_UNKNOWN( end_        == other.end()   );

		// merge left and right into target set
		size_t left_value  = data_[this_left].value;
		size_t left_next   = data_[this_left].next;
		size_t right_value = other.data_[other_right].value;
		size_t right_next  = other.data_[other_right].next;
		size_t current       = this_target;
		while( (left_value < end_) | (right_value < end_) )
		{	if( left_value == right_value )
			{	right_value = other.data_[right_next].value;
				right_next  = other.data_[right_next].next;
			}
			size_t next         = data_.extend(1);
			data_[current].next = next;
			if( left_value < right_value )
			{	data_[current].value = left_value;
				left_value           = data_[left_next].value;
				left_next            = data_[left_next].next;
			}
			else
			{	data_[current].value = right_value;
				right_value = other.data_[right_next].value;
				right_next  = other.data_[right_next].next;
			}
			current = next;
		}
		data_[current].value = end_;
	}
	// -----------------------------------------------------------------
	/*! Sum over all sets of the number of elements
 
	/return
	The the total number of elements
 	*/
	size_t number_elements(void) const
	{	size_t i, count;
		count = 0;
		for(i = 0; i < n_set_; i++)
		{	size_t index = i;
			while( data_[index].value != end_ )
			{	count++;
				index = data_[index].next;
			}
		}
		return count;
	}
	// -----------------------------------------------------------------
	/*! Fetch n_set for vector of sets object.
	
	\return
	Number of from sets for this vector of sets object
	*/
	size_t n_set(void) const
	{	return n_set_; }
	// -----------------------------------------------------------------
	/*! Fetch end for this vector of sets object. 
	
	\return
	is the maximum element value plus one (the minimum element value is 0).
	*/
	size_t end(void) const
	{	return end_; }
};

/*! 
Copy a user vector of sets sparsity pattern to an internal sparse_list object.

\tparam VectorSet
is a simple vector with elements of type \c std::list<set>.

\param internal
The input value of sparisty does not matter.
Upon return it contains the same sparsity pattern as \c user
(or the transposed sparsity pattern).

\param user
sparsity pattern that we are placing \c internal.

\param n_row
number of rows in the sparsity pattern in \c user
(range dimension).

\param n_col
number of columns in the sparsity pattern in \c user
(domain dimension).

\param transpose
if true, the sparsity pattern in \c internal is the transpose
of the one in \c user. 
Otherwise it is the same sparsity pattern.
*/
template<class VectorSet>
void sparsity_user2internal(
	sparse_list&            internal  , 
	const VectorSet&        user      ,
	size_t                  n_row     ,
	size_t                  n_col     ,
	bool                    transpose )
{	CPPAD_ASSERT_UNKNOWN( n_row == size_t(user.size()) );

	CPPAD_ASSERT_KNOWN(
		size_t( user.size() ) == n_row,
		"Size of this vector of lists sparsity pattern is not equal "
		"the range dimension for the corresponding function."
	);

	size_t i, j;
	std::set<size_t>::const_iterator itr;

	// transposed pattern case
	if( transpose )
	{	internal.resize(n_col, n_row);
		for(i = 0; i < n_row; i++)
		{	itr = user[i].begin();
			while(itr != user[i].end())
			{	j = *itr++;
				CPPAD_ASSERT_UNKNOWN( j < n_col );
				internal.add_element(j, i);
			}
		}
		return;
	}

	// same pattern case
	internal.resize(n_row, n_col);
	for(i = 0; i < n_row; i++)
	{	itr = user[i].begin();
		while(itr != user[i].end())
		{	j = *itr++;
			CPPAD_ASSERT_UNKNOWN( j < n_col );
			internal.add_element(i, j);
		}
	}
	return;
}


/*! \} */
} // END_CPPAD_NAMESPACE
# endif
