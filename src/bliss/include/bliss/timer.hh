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


#include <chrono>

namespace bliss {

/**
 * \brief A simple helper class for measuring elapsed time.
 */
class Timer
{
  std::chrono::high_resolution_clock::time_point start_time_;
public:
  /**
   * \brief Create and start a new timer.
   */
  Timer() {
    reset();
  }

  /**
   * \brief Reset the timer.
   */
  void reset() {
    start_time_ = std::chrono::high_resolution_clock::now();
  }

  /**
   * \brief Get the time (in seconds) elapsed since
   *        the creation or the last reset() call of the timer.
   */
  double get_duration() const {
    return std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::high_resolution_clock::now() - start_time_).count();
  }
};

} // namespace bliss
