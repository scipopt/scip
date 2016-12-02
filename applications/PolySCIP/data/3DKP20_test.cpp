#include <algorithm>
#include <assert.h>
#include <functional>
#include <iostream>
#include <list>
#include <numeric> // std::inner_product
#include <vector>

using std::list;
using std::vector;

  vector<int> obj1 = {12, 11, 60, 100, 34, 63, 10, 19, 92, 52,
		   61, 60, 80, 96, 30, 77, 92, 26, 69, 21};

  vector<int> obj2 = {38, 56, 63, 29,  41, 83, 85, 74, 21, 17,
		   59, 49, 69, 49, 69, 60, 17, 74, 100, 37};

  vector<int> obj3 = {69, 64, 34, 83, 31, 10, 25, 53, 92, 24,
		   92, 96, 94, 46, 21, 89, 84, 25, 92,  66};

  vector<int> cap1 = {20, 75, 55, 66, 12, 73, 95, 17, 41, 19,
		   95, 13, 30, 60, 74, 97, 20, 68, 49, 84};

  vector<int> cap2 = {95, 51, 53, 23, 97, 45, 85, 21, 59, 31,
		   74, 62, 96, 95, 70, 25, 66, 59, 80, 56};

  vector<int> cap3 = {76, 44, 96, 32, 91, 72, 64, 93, 96, 19,
		   72, 52, 15, 47, 82, 83, 45, 62, 68, 58};

  double weight1 = 531.5;

  double weight2 = 621.5;

  double weight3 = 633.5;

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
    auto cap1_sum = std::inner_product(begin(base),
					end(base),
					begin(cap1),
					0.);
    auto cap2_sum = std::inner_product(begin(base),
				       end(base),
				       begin(cap2),
				       0.);
    auto cap3_sum = std::inner_product(begin(base),
				       end(base),
				       begin(cap3),
				       0.);
    if (cap1_sum <= weight1 &&
	cap2_sum <= weight2 &&
	cap3_sum <= weight3) {
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
  assert (obj1.size() == 20);
  assert (obj2.size() == 20);
  assert (obj3.size() == 20);
  assert (cap1.size() == 20);
  assert (cap2.size() == 20);
  assert (cap3.size() == 20);

  auto feasible_solutions = vector<vector<int>>{};
  auto sol_size = 20;
  for (auto i=1; i<=sol_size; ++i) {
    add_vecs_with_ones(vector<int>(sol_size,0),
		       i, 0, feasible_solutions);
  }
  std::cout << "NUMBER OF FEASIBLE SOLS: " << feasible_solutions.size() << "\n";
  auto feasible_outcomes = list<vector<int>>{};
  for (const auto& sol : feasible_solutions) {
    feasible_outcomes.push_back(getFeasibleOutcome(sol));
  }

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
				  std::less_equal<int>());
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
