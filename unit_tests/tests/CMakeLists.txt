# Quantum++
#
# Copyright (c) 2013 - 2017 Vlad Gheorghiu (vgheorgh@gmail.com)
#
# This file is part of Quantum++.
#
# Quantum++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Quantum++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.

INCLUDE_DIRECTORIES(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
ADD_EXECUTABLE(qpp_testing
        classes/gates.cpp
        classes/timer.cpp
        entanglement.cpp
        entropies.cpp
        functions.cpp
        input_output.cpp
        instruments.cpp
        number_theory.cpp
        operations.cpp
        random.cpp
        statistics.cpp
        testing_main.cpp
        traits.cpp)
TARGET_LINK_LIBRARIES(qpp_testing gmock)
