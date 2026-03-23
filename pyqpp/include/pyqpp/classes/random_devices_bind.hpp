/*
 * This file is part of pyqpp.
 *
 * Copyright (c) 2017 - 2026 softwareQ Inc. All rights reserved.
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

/**
 * @file <pyqpp/classes/random_devices_bind.hpp>
 * @brief Bindings for <qpp/classes/random_devices.hpp>
 */

#ifndef PYQPP_CLASSES_RANDOM_DEVICES_BIND_HPP_
#define PYQPP_CLASSES_RANDOM_DEVICES_BIND_HPP_

#include "pyqpp/pyqpp_common.hpp"

inline void init_classes_random_devices(py::module_& m) {
    using namespace qpp;

    auto random_devices = m.def_submodule("random_devices");
    // We expose the functionality of the RandomDevices singleton
    // directly as functions in the 'random' submodule

    // Wrap the mt19937 type so Python recognizes the return type of
    // get_prng()
    py::class_<std::mt19937>(random_devices, "mt19937")
        .def(
            "seed", [](std::mt19937& self, unsigned int s) { self.seed(s); },
            py::arg("s"))
        .def(
            "__call__", [](std::mt19937& self) { return self(); },
            "Generate a random number");

    random_devices.def(
        "get_prng",
        []() -> std::mt19937& {
            return RandomDevices::get_instance().get_prng();
        },
        py::return_value_policy::reference, // Crucial: don't let Python delete
                                            // the global PRNG
        "Returns a reference to the internal Mersenne Twister (mt19937) "
        "object");

    random_devices.def(
        "save",
        []() {
            std::ostringstream oss;
            RandomDevices::get_instance().save(oss);
            return oss.str();
        },
        "Saves the current state of the PRNG to a string");

    random_devices.def(
        "load",
        [](const std::string& state) {
            std::istringstream iss(state);
            RandomDevices::get_instance().load(iss);
        },
        "Loads the PRNG state from a string", py::arg("state"));

    // Seed the PRNG manually by accessing the underlying mt19937
    random_devices.def(
        "seed",
        [](unsigned int seed_val) {
            RandomDevices::get_instance().get_prng().seed(seed_val);
        },
        "Seeds the internal Mersenne Twister PRNG", py::arg("seed_val"));
}

#endif /* PYQPP_CLASSES_RANDOM_DEVICES_BIND_HPP_ */
