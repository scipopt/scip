#pragma once

/*
  Copyright (c) 2003-2021 Tommi Junttila
  Released under the GNU Lesser General Public License version 3.

  This file is part of bliss.

  bliss is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, version 3 of the License.

  bliss is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with bliss.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include <algorithm>

namespace bliss {

/**
 * \brief A min-heap of unsigned integers.
 */

class Heap
{
  std::vector<unsigned int> contents;
  /** \nternal
   * Less-than operator for building min-heaps.
   */
  struct {
    /** \internal Use greater-than as less-than -> min-heaps. */
    bool operator()(const unsigned int e1, const unsigned int e2) {return e1 > e2; }
  } gt;
public:

  /**
   * Is the heap empty?
   * Time complexity is O(1).
   */
  bool is_empty() const {return contents.empty(); }

  /**
   * Remove all the elements in the heap.
   * Time complexity is O(1).
   */
  void clear() {contents.clear(); }

  /**
   * Insert the element \a e in the heap.
   * Time complexity is O(log(N)), where N is the number of elements
   * currently in the heap.
   */
  void insert(const unsigned int e) {
    contents.push_back(e);
    std::push_heap(contents.begin(), contents.end(), gt);
  }

  /**
   * Return the smallest element in the heap.
   * Time complexity is O(1).
   */
  unsigned int smallest() const {return contents.front(); }

  /**
   * Remove and return the smallest element in the heap.
   * Time complexity is O(log(N)), where N is the number of elements
   * currently in the heap.
   */
  unsigned int remove() {
    const unsigned int result = smallest();
    std::pop_heap(contents.begin(),contents.end(), gt);
    contents.pop_back();
    return result;
  }

  /**
   * Get the number of elements in the heap.
   */
  size_t size() const {return contents.size(); }
};

} // namespace bliss
