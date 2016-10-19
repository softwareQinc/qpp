/*
 * Quantum++
 *
 * Copyright (c) 2013 - 2017 Vlad Gheorghiu (vgheorgh@gmail.com)
 *
 * This file is part of Quantum++.
 *
 * Quantum++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Quantum++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "classes/states.h"

/******************************************************************************/
/// BEGIN qpp::States::mes(idx d = 2) const
TEST(qpp_States_mes, AllTests)
{
    // d = 1 (number)
    ket ket_one(1);
    ket_one << 1;
    EXPECT_NEAR(0, norm(qpp::st.mes(1) - ket_one), 1e-7);

    // qubits
    EXPECT_NEAR(0, norm(qpp::st.mes() - st.b00), 1e-7);

    // qutrits
    idx d = 3;
    ket psi = mket({0, 0}, {d, d}) / std::sqrt(d);
    for (idx i = 1; i < d; ++i)
    {
        psi += mket({i, i}, {d, d});
    }
    ket mes_qutrit(3, 3);
    mes_qutrit << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    mes_qutrit /= std::sqrt(d);
    EXPECT_NEAR(0, norm(qpp::st.mes(d) - mes_qutrit), 1e-7);

    // ququads
    d = 4;
    psi = mket({0, 0}, {d, d}) / std::sqrt(d);
    for (idx i = 1; i < d; ++i)
    {
        psi += mket({i, i}, {d, d});
    }
    ket mes_ququad(4, 4);
    mes_ququad << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
    mes_ququad /= std::sqrt(d);
    EXPECT_NEAR(0, norm(qpp::st.mes(d) - mes_ququad), 1e-7);
}
/******************************************************************************/
