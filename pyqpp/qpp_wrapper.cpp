/*
 * This file is part of pyqpp.
 *
 * Copyright (c) 2019 - 2023 softwareQ Inc. All rights reserved.
 *
 * MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "pyqpp_common.h"

#include "constants_bind.hpp"
#include "functions_bind.hpp"
#include "instruments_bind.hpp"
#include "random_bind.hpp"

#include "classes/gates_bind.hpp"
#include "classes/reversible_bind.hpp"
#include "classes/states_bind.hpp"

#include "classes/circuits/circuits_bind.hpp"
#include "classes/circuits/engines_bind.hpp"

#include "qasm/qasm_bind.hpp"

#include "pyqpp_specific_bind.hpp"

PYBIND11_MODULE(pyqpp, m) {
    m.doc() =
        "Python 3 wrapper for Quantum++ (https://github.com/softwareQinc/qpp)";

    init_constants(m);
    init_functions(m);
    init_instruments(m);
    init_random(m);

    init_classes_gates(m);
    init_classes_reversible(m);
    init_classes_states(m);

    init_classes_circuits_circuits(m);
    init_classes_circuits_engines(m);

    init_qasm_qasm(m);

    init_pyqpp_specific(m);
}
