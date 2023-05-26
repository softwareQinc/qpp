#include <vector>

#include "gtest/gtest.h"

#include "qpp.h"

using namespace qpp;

// Unit testing "classes/circuits/circuits.hpp"

/******************************************************************************/
/// BEGIN QCircuit& QCircuit::add_circuit(QCircuit other, bigint pos_qudit,
///       std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_QCircuit_add_circuit, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::add_dit(idx n = 1, idx i = -1)
TEST(qpp_QCircuit_add_dit, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::add_qudit(idx n = 1, idx i = -1)
TEST(qpp_QCircuit_add_qudit, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::adjoint()
TEST(qpp_QCircuit_adjoint, AllTests) {}
/******************************************************************************/
/// BEGIN iterator QCircuit::begin()
TEST(qpp_QCircuit_begin, Iterator) {}
/******************************************************************************/
/// BEGIN const_iterator QCircuit::begin() const noexcept
TEST(qpp_QCircuit_begin, ConstIterator) {}
/******************************************************************************/
/// BEGIN const_iterator QCircuit::cbegin() const noexcept
TEST(qpp_QCircuit_cbegin, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::cCTRL(const cmat& U,
///       const std::vector<idx>& ctrl_dits, const std::vector<idx>& target,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_cCTRL, MultipleCtrlsMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::cCTRL(const cmat& U,
///       const std::vector<idx>& ctrl_dits, idx target,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_cCTRL, MultipleCtrlsSingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::cCTRL(const cmat& U, idx ctrl_dit,
///       const std::vector<idx>& target,
///       std::optional<idx> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_cCTRL, SingleCtrlMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::cCTRL(const cmat& U, idx ctrl_dit, idx target,
///       std::optional<idx> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_cCTRL, SingleCtrlSingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::cCTRL_fan(const cmat& U,
///       const std::vector<idx>& ctrl_dits, const std::vector<idx>& target,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_cCTRL_fan, MultipleCtrlsMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::cCTRL_fan(const cmat& U, idx ctrl_dit,
///       const std::vector<idx>& target,
///       std::optional<idx> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_cCTRL_fan, SingleCtrlMultipleTargets) {}
/******************************************************************************/
/// BEGIN const_iterator QCircuit::cend() const noexcept
TEST(qpp_QCircuit_cend, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::compress()
TEST(qpp_QCircuit_compress, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::CTRL(const cmat& U,
///       const std::vector<idx>& ctrl, const std::vector<idx>& target,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_CTRL, MultipleCtrlsMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::CTRL(const cmat& U,
///       const std::vector<idx>& ctrl, idx target,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_CTRL, MultipleCtrlsSingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::CTRL(const cmat& U, idx ctrl,
///       const std::vector<idx>& target,
///       std::optional<idx> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_CTRL, SingleCtrlMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::CTRL(const cmat& U, idx ctrl, idx target,
///       std::optional<idx> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_CTRL, SingleCtrlSingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::CTRL_fan(const cmat& U,
///       const std::vector<idx>& ctrl, const std::vector<idx>& target,
///       std::optional<std::vector<idx>> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_CTRL_fan, MultipleCtrlsMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::CTRL_fan(const cmat& U, idx ctrl,
///       const std::vector<idx>& target,
///       std::optional<idx> shift = std::nullopt,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_CTRL_fan, SingleCtrlMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::discard(const std::vector<idx>& target,
///       std::optional<std::string> name = "discard")
TEST(qpp_QCircuit_discard, MultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::discard(idx target,
///       std::optional<std::string> name = "discard"")
TEST(qpp_QCircuit_discard, SingleTarget) {}
/******************************************************************************/
/// BEGIN iterator QCircuit::end()
TEST(qpp_QCircuit_end, Iterator) {}
/******************************************************************************/
/// BEGIN const_iterator QCircuit::end() const noexcept
TEST(qpp_QCircuit_end, ConstIterator) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::gate(const cmat& U, idx i, idx j, idx k,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_gate, ThreeQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::gate(const cmat& U, idx i, idx j,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_gate, TwoQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::gate(const cmat& U, idx i,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_gate, SingleQudit) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::gate(const cmat& U,
///       const std::vector<idx>& target,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_gate, JointQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::gate_fan(const cmat& U,
///       const std::vector<idx>& target,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_gate_fan, SpecificQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::gate_fan(const cmat& U,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_gate_fan, AllQudits) {}
/******************************************************************************/
/// BEGIN std::vector<idx> QCircuit::get_clean_dits() const
TEST(qpp_QCircuit_get_clean_dits, AllTests) {}
/******************************************************************************/
/// BEGIN std::vector<idx> QCircuit::get_clean_qudits() const
TEST(qpp_QCircuit_get_clean_qudits, AllTests) {}
/******************************************************************************/
/// BEGIN const std::unordered_map<std::size_t, cmat>&
///       QCircuit::get_cmat_hash_tbl() const noexcept
TEST(qpp_QCircuit_get_cmat_hash_tbl, AllTests) {}
/******************************************************************************/
/// BEGIN idx QCircuit::get_d() const noexcept
TEST(qpp_QCircuit_get_d, AllTests) {}
/******************************************************************************/
/// BEGIN idx QCircuit::get_depth() const
TEST(qpp_QCircuit_get_depth, AllTests) {}
/******************************************************************************/
/// BEGIN std::vector<idx> QCircuit::get_dirty_dits() const
TEST(qpp_QCircuit_get_dirty_dits, AllTests) {}
/******************************************************************************/
/// BEGIN std::vector<idx> QCircuit::get_dirty_qudits() const
TEST(qpp_QCircuit_get_dirty_qudits, AllTests) {}
/******************************************************************************/
/// BEGIN idx QCircuit::get_gate_count(std::optional<cmat> U = std::nullopt)
///       const
TEST(qpp_QCircuit_get_gate_count, TotalGateCount) {}
TEST(qpp_QCircuit_get_gate_count, SpecificGateCount) {}
/******************************************************************************/
/// BEGIN idx QCircuit::get_gate_depth(std::optional<cmat> U = std::nullopt)
///       const
TEST(qpp_QCircuit_get_gate_depth, TotalGateDepth) {}
TEST(qpp_QCircuit_get_gate_depth, SpecificGateDepth) {}
/******************************************************************************/
/// BEGIN std::vector<idx> QCircuit::get_measured() const
TEST(qpp_QCircuit_get_measured, AllTests) {}
/******************************************************************************/
/// BEGIN std::vector<idx> QCircuit::get_measured_nd() const
TEST(qpp_QCircuit_get_measured_nd, AllTests) {}
/******************************************************************************/
/// BEGIN idx QCircuit::get_measurement_count(std::optional<cmat> V =
///       std::nullopt) const
TEST(qpp_QCircuit_get_measurement_count, TotalMeasurementCount) {}
TEST(qpp_QCircuit_get_measurement_count, SpecificMeasurementCount) {}
/******************************************************************************/
/// BEGIN idx QCircuit::get_measurement_depth(std::optional<cmat> V =
///       std::nullopt) const
TEST(qpp_QCircuit_get_measurement_depth, TotalMeasurementDepth) {}
TEST(qpp_QCircuit_get_measurement_depth, SpecificMeasurementDepth) {}
/******************************************************************************/
/// BEGIN std::vector<idx> QCircuit::get_measurement_dits() const
TEST(qpp_QCircuit_get_measurement_dits, AllTests) {}
/******************************************************************************/
/// BEGIN std::optional<std::string> QCircuit::get_name() const
TEST(qpp_QCircuit_get_name, AllTests) {}
/******************************************************************************/
/// BEGIN idx QCircuit::get_nc() const noexcept
TEST(qpp_QCircuit_get_nc, AllTests) {}
/******************************************************************************/
/// BEGIN std::vector<idx> QCircuit::get_non_measured() const
TEST(qpp_QCircuit_get_non_measured, AllTests) {}
/******************************************************************************/
/// BEGIN idx QCircuit::get_nop_count() const
TEST(qpp_QCircuit_get_nop_count, AllTests) {}
/******************************************************************************/
/// BEGIN idx QCircuit::get_nq() const noexcept
TEST(qpp_QCircuit_get_nq, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit::Resources QCircuit::get_resources() const
TEST(qpp_QCircuit_get_resources, AllTests) {}
/******************************************************************************/
/// BEGIN idx QCircuit::get_step_count() const noexcept
TEST(qpp_QCircuit_get_step_count, AllTests) {}
/******************************************************************************/
/// BEGIN bool QCircuit::has_measurements() const noexcept
TEST(qpp_QCircuit_has_measurements, AllTests) {}
/******************************************************************************/
/// BEGIN inline static bool QCircuit::is_cCTRL(const GateStep& gate_step)
TEST(qpp_QCircuit_is_cCTRL, AllTests) {}
/******************************************************************************/
/// BEGIN bool QCircuit::is_clean_dit(idx i) const
TEST(qpp_QCircuit_is_clean_dit, AllTests) {}
/******************************************************************************/
/// BEGIN bool QCircuit::is_clean_qudit(idx i) const
TEST(qpp_QCircuit_is_clean_qudit, AllTests) {}
/******************************************************************************/
/// BEGIN inline static bool QCircuit::is_CTRL(const GateStep& gate_step)
TEST(qpp_QCircuit_is_CTRL, AllTests) {}
/******************************************************************************/
/// BEGIN bool QCircuit::is_measurement_dit(idx i) const
TEST(qpp_QCircuit_is_measurement_dit, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::kron(QCircuit qc)
TEST(qpp_QCircuit_kron, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::match_circuit_left(QCircuit other,
///       const std::vector<idx>& target,
///       std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_QCircuit_match_circuit_left, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::match_circuit_right(QCircuit other,
///       const std::vector<idx>& target,
///       std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_QCircuit_match_circuit_right, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::measureV(const cmat& V,
///       const std::vector<idx>& target, idx c_reg, bool destructive = true,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_measureV, JointMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::measureV(const cmat& V, idx target, idx c_reg,
///       bool destructive = true,
///       std::optional<std::string> name = std::nullopt)
TEST(qpp_QCircuit_measureV, SingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::measure(const std::vector<idx>& target,
///       idx c_reg = 0, bool destructive = true,
///       std::optional<std::string> name = "mZ")
TEST(qpp_QCircuit_measure, MultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::measure(idx target, idx c_reg,
///       bool destructive = true, std::optional<std::string> name = "mZ")
TEST(qpp_QCircuit_measure, SingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::measure_all(idx c_reg = 0,
///       bool destructive = true, std::optional<std::string> name = "mZ")
TEST(qpp_QCircuit_measure_all, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::nop()
TEST(qpp_QCircuit_nop, AllTests) {}
/******************************************************************************/
/// BEGIN bool QCircuit::operator!=(const QCircuit& rhs) const noexcept
TEST(qpp_QCircuit_operator_noneq, AllTests) {}
/******************************************************************************/
/// BEGIN bool QCircuit::operator==(const QCircuit& rhs) const noexcept
TEST(qpp_QCircuit_operator_eq, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::QFT(bool swap = true)
TEST(qpp_QCircuit_QFT, AllQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::QFT(const std::vector<idx>& target,
///       bool swap = true)
TEST(qpp_QCircuit_QFT, SpecificQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::remove_clean_dit(idx target)
TEST(qpp_QCircuit_remove_clean_dit, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::remove_clean_dits(std::vector<idx> target)
TEST(qpp_QCircuit_remove_clean_dits, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::remove_clean_qudit(idx target)
TEST(qpp_QCircuit_remove_clean_qudit, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::remove_clean_qudits(std::vector<idx> target)
TEST(qpp_QCircuit_remove_clean_qudits, AllTests) {}
/******************************************************************************/
/// BEGIN bool QCircuit::removes_qudits() const noexcept
TEST(qpp_QCircuit_removes_qudits, AllTests) {}
/// BEGIN QCircuit& QCircuit::replicate(idx n)
TEST(qpp_QCircuit_replicate, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::reset(const std::vector<idx>& target,
///       std::optional<std::string> name = "reset")
TEST(qpp_QCircuit_reset, MultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::reset(idx target,
///       std::optional<std::string> name = "reset")
TEST(qpp_QCircuit_reset, SingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::set_name(const std::string& name)
TEST(qpp_QCircuit_set_name, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::TFQ(bool swap = true)
TEST(qpp_QCircuit_TFQ, AllQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& QCircuit::TFQ(const std::vector<idx>& target,
///       bool swap = true)
TEST(qpp_QCircuit_TFQ, SpecificQudits) {}
/******************************************************************************/
/// BEGIN std::string QCircuit::to_JSON(
///       bool enclosed_in_curly_brackets = true)
///       const override
TEST(qpp_QCircuit_to_JSON, AllTests) {}
/******************************************************************************/
/// BEGIN bool QCircuit::was_measured(idx i) const
TEST(qpp_QCircuit_was_measured, AllTests) {}
/******************************************************************************/
/// BEGIN bool QCircuit::was_measured_nd(idx i) const
TEST(qpp_QCircuit_was_measured_nd, AllTests) {}
/******************************************************************************/

// free functions
/******************************************************************************/
/// BEGIN inline QCircuit add_circuit(QCircuit qc1, const QCircuit& qc2,
///       bigint pos_qudit, std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_add_circuit, AllTests) {}
/******************************************************************************/
/// BEGIN inline QCircuit adjoint(QCircuit qc)
TEST(qpp_adjoint, QCircuitAllTests) {}
/******************************************************************************/
/// BEGIN inline QCircuit kron(QCircuit qc1, const QCircuit& qc2)
TEST(qpp_kron, QCircuitAllTests) {}
/******************************************************************************/
/// BEGIN inline QCircuit qpe_circuit(cmat U, qpp::idx n,
///       bool omit_measurements = true, idx d = 2,
///       std::optional<std::string> name = "qpe")
TEST(qpp_qpe_circuit, AllTests) {}
/******************************************************************************/
/// BEGIN inline QCircuit match_circuit_left(QCircuit qc1, const QCircuit& qc2,
///       const std::vector<idx>& target,
///       std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_match_circuit_left, AllTests) {}
/******************************************************************************/
/// BEGIN inline QCircuit match_circuit_right(QCircuit qc1, const QCircuit& qc2,
///       const std::vector<idx>& target,
///       std::optional<idx> pos_dit = std::nullopt)
TEST(qpp_match_circuit_right, AllTests) {}
/******************************************************************************/
/// BEGIN inline QCircuit replicate(QCircuit qc, idx n)
TEST(qpp_replicate, AllTests) {}
/******************************************************************************/
/// BEGIN inline QCircuit random_circuit_count(idx nq, idx d, idx gate_count,
///       std::optional<realT> p_two,
///       std::optional<cmat> with_respect_to_gate = std::nullopt,
///       std::optional<std::vector<cmat>> one_qudit_gate_set = std::nullopt,
///       std::optional<std::vector<cmat>> two_qudit_gate_set = std::nullopt,
///       std::optional<std::vector<std::string>> one_qudit_gate_names =
///           std::nullopt,
///       std::optional<std::vector<std::string>> two_qudit_gate_names =
///           std::nullopt)
TEST(qpp_random_circuit_count, AllTests) {}
/******************************************************************************/
/// BEGIN inline QCircuit random_circuit_depth(idx nq, idx d, idx gate_depth,
///       std::optional<realT> p_two,
///       std::optional<cmat> with_respect_to_gate = std::nullopt,
///       std::optional<std::vector<cmat>> one_qudit_gate_set = std::nullopt,
///       std::optional<std::vector<cmat>> two_qudit_gate_set = std::nullopt,
///       std::optional<std::vector<std::string>> one_qudit_gate_names =
///           std::nullopt,
///       std::optional<std::vector<std::string>> two_qudit_gate_names =
///           std::nullopt)
TEST(qpp_random_circuit_depth, AllTests) {}
/******************************************************************************/
