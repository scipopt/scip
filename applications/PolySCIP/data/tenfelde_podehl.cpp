#include <algorithm>
#include <assert.h>
#include <functional>
#include <iostream>
#include <list>
#include <numeric>
#include <vector>

using std::vector;
using std::list;
using std::begin;
using std::end;

vector<int> obj1 = {3, 6, 4, 5, 2, 3, 5, 4, 3, 5, 4, 2, 4, 5, 3, 6};

vector<int> obj2 = {2, 3, 5, 4, 5, 3, 4, 3, 5, 2, 6, 4, 4, 5, 2, 5};

vector<int> obj3 = {4, 2, 4, 2, 4, 2, 4, 6, 4, 2, 6, 3, 2, 4, 5, 3};

void print_vec(const vector<int>& vec) {
  std::cout << "[ ";
  for (auto i : vec)
    std::cout << i << " ";
  std::cout << "]\n";
}

vector<int> getFeasibleOutcome(const vector<int>& sol) {
  auto outcome = vector<int>{};
  outcome.push_back(std::inner_product(begin(sol),
				       end(sol),
				       begin(obj1),
				       0.));
  outcome.push_back(std::inner_product(begin(sol),
				       end(sol),
				       begin(obj2),
				       0.));
  outcome.push_back(std::inner_product(begin(sol),
				       end(sol),
				       begin(obj3),
				       0.));
  return outcome;
}

void add_vecs_with_ones(vector<int> base,
			int no_of_ones,
			int start,
			vector<vector<int>>& all) {
  if (no_of_ones == 0) {
    auto obj3_sum_const = std::inner_product(begin(base),
					     end(base),
					     begin(obj3),
					     0.);
    if (obj3_sum_const <= 15) {
      all.push_back(base);
    }
  }
  else {
    for (auto i=start; i<base.size(); ++i) {
      auto changed = base;
      changed[i] = 1;
      add_vecs_with_ones(std::move(changed), no_of_ones-1, i+1, all);
    }
  }
}

int main() {
  assert (obj1.size() == 16);
  assert (obj2.size() == 16);
  assert (obj3.size() == 16);

  auto feasible_outcomes = list<vector<int>>{};

  for (int first = 0; first<4; ++first) {
    for (int sec = 0; sec<4; ++sec) {
      if (sec == first)
	continue;
      for (int third = 0; third<4; ++third) {
	if (third == first || third == sec)
	  continue;
	for (int forth = 0; forth<4; ++forth) {
	  if (forth == third || forth == sec || forth == first)
	    continue;
	  auto sol = vector<int>(16,0.);
	  sol[first] = 1;
	  sol[4+sec] = 1;
	  sol[8+third] = 1;
	  sol[12+forth] = 1;
	  if (std::inner_product(begin(sol),
				 end(sol),
				 begin(obj2),
				 0.) <= 15) {
	    feasible_outcomes.push_back( getFeasibleOutcome(sol) );
	  }
	}
      }
    }
  }

  std::cout << "NUMBER OF FEASIBLE OUTCOMES: " << feasible_outcomes.size() << "\n";

  auto it = begin(feasible_outcomes);
  while (it != end(feasible_outcomes)) {
    auto dominated_eq = false;
    for (auto other=begin(feasible_outcomes);
	 other!=end(feasible_outcomes); ++other) {
      if (other == it) {
	continue;
      }
      else {
	dominated_eq = std::equal(begin(*it),
				  end(*it),
				  begin(*other),
				  std::greater_equal<int>());
      }
      if (dominated_eq)
	break;
    }
    if (dominated_eq) {
      it = feasible_outcomes.erase(it);
    }
    else {
      ++it;
    }
  }
  // output nondominated outcomes
  std::cout << "NONDOMINATED OUTCOMES: \n";
  for (const auto& outcome : feasible_outcomes) {
    print_vec(outcome);
  }
  return 0;
}
