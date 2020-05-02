#include <vector>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "classes/circuits/circuits.h"

/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::add_circuit(QCircuit other, bigint pos_qudit,
///       idx pos_dit = -1)
TEST(qpp_QCircuit_add_circuit, AllTests) {}
/******************************************************************************/
/// BEGIN friend QCircuit qpp::add_circuit(QCircuit qc1, const QCircuit& qc2,
///       bigint pos_qudit, idx pos_dit = -1)
TEST(qpp_friend_QCircuit_add_circuit, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::add_dit(idx n = 1, idx i = -1)
TEST(qpp_QCircuit_add_dit, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::add_qudit(idx n = 1, idx i = -1)
TEST(qpp_QCircuit_add_qudit, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::adjoint()
TEST(qpp_QCircuit_adjoint, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit qpp::adjoint(QCircuit qc)
TEST(qpp_friend_QCircuit_adjoint, AllTests) {}
/******************************************************************************/
/// BEGIN iterator qpp::QCircuit::begin()
TEST(qpp_QCircuit_begin, Iterator) {}
/******************************************************************************/
/// BEGIN const_iterator qpp::QCircuit::begin() const noexcept
TEST(qpp_QCircuit_begin, ConstIterator) {}
/******************************************************************************/
/// BEGIN const_iterator qpp::QCircuit::cbegin() const noexcept
TEST(qpp_QCircuit_cbegin, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::cCTRL(const cmat& U,
///       const std::vector<idx>& ctrl_dits, const std::vector<idx>& target,
///       const std::vector<idx>& shift = {}, std::string name = {})
TEST(qpp_QCircuit_cCTRL, MultipleCtrlsMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::cCTRL(const cmat& U,
///       const std::vector<idx>& ctrl_dits, idx target,
///       const std::vector<idx>& shift = {}, std::string name = {})
TEST(qpp_QCircuit_cCTRL, MultipleCtrlsSingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::cCTRL(const cmat& U, idx ctrl_dit,
///       const std::vector<idx>& target, idx shift = 0, std::string name = {})
TEST(qpp_QCircuit_cCTRL, SingleCtrlMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::cCTRL(const cmat& U, idx ctrl_dit,
///       idx target, idx shift = 0, std::string name = {})
TEST(qpp_QCircuit_cCTRL, SingleCtrlSingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::cCTRL_custom(const cmat& U,
///       const std::vector<idx>& ctrl_dits, const std::vector<idx>& target,
///       const std::vector<idx>& shift = {}, std::string name = {})
TEST(qpp_QCircuit_cCTRL_custom, AllTests) {}
/******************************************************************************/
/// BEGIN const_iterator qpp::QCircuit::cend() const noexcept
TEST(qpp_QCircuit_cend, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::compress()
TEST(qpp_QCircuit_compress, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::CTRL(const cmat& U,
///       const std::vector<idx>& ctrl, const std::vector<idx>& target,
///       const std::vector<idx>& shift = {}, std::string name = {})
TEST(qpp_QCircuit_CTRL, MultipleCtrlsMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::CTRL(const cmat& U,
///       const std::vector<idx>& ctrl, idx target,
///       const std::vector<idx>& shift = {}, std::string name = {})
TEST(qpp_QCircuit_CTRL, MultipleCtrlsSingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::CTRL(const cmat& U, idx ctrl,
///       const std::vector<idx>& target, idx shift = 0, std::string name = {})
TEST(qpp_QCircuit_CTRL, SingleCtrlMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::CTRL(const cmat& U, idx ctrl, idx target,
///       idx shift = 0, std::string name = {})
TEST(qpp_QCircuit_CTRL, SingleCtrlSingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::CTRL_custom(const cmat& U,
///       const std::vector<idx>& ctrl,const std::vector<idx>& target,
///       const std::vector<idx>& shift = {}, std::string name = {})
TEST(qpp_QCircuit_CTRL_custom, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::discard(const std::vector<idx>& target,
///       std::string name = {})
TEST(qpp_QCircuit_discard, MultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::discard(idx target, std::string name = {})
TEST(qpp_QCircuit_discard, SingleTarget) {}
/******************************************************************************/
/// BEGIN iterator qpp::QCircuit::end()
TEST(qpp_QCircuit_end, Iterator) {}
/******************************************************************************/
/// BEGIN const_iterator qpp::QCircuit::end() const noexcept
TEST(qpp_QCircuit_end, ConstIterator) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::gate(const cmat& U, idx i, idx j, idx k,
///       std::string name = {})
TEST(qpp_QCircuit_gate, ThreeQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::gate(const cmat& U, idx i, idx j,
///       std::string name = {})
TEST(qpp_QCircuit_gate, TwoQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::gate(const cmat& U, idx i,
///       std::string name = {})
TEST(qpp_QCircuit_gate, SingleQudit) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::gate_custom(const cmat& U,
///       const std::vector<idx>& target, std::string name = {})
TEST(qpp_QCircuit_gate_custom, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::gate_fan(const cmat& U,
///       const std::vector<idx>& target, std::string name = {})
TEST(qpp_QCircuit_gate_fan, SpecificQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::gate_fan(const cmat& U,
///       std::string name = {})
TEST(qpp_QCircuit_gate_fan, AllQudits) {}
/******************************************************************************/
/// BEGIN std::vector<idx> qpp::QCircuit::get_clean() const
TEST(qpp_QCircuit_get_clean, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_d() const noexcept
TEST(qpp_QCircuit_get_d, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_depth() const
TEST(qpp_QCircuit_get_depth, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_gate_count() const
TEST(qpp_QCircuit_get_gate_count, TotalGateCount) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_gate_count(const std::string& name) const
TEST(qpp_QCircuit_get_gate_count, SpecificGateCount) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_gate_depth() const
TEST(qpp_QCircuit_get_gate_depth, TotalGateDepth) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_gate_depth(const std::string& name) const
TEST(qpp_QCircuit_get_gate_depth, SpecificGateDepth) {}
/******************************************************************************/
/// BEGIN std::vector<idx> qpp::QCircuit::get_measured() const
TEST(qpp_QCircuit_get_measured, AllQudits) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_measured(idx i) const
TEST(qpp_QCircuit_get_measured, SpecificQudit) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_measurement_count() const
TEST(qpp_QCircuit_get_measurement_count, TotalMeasurementCount) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_measurement_count(
///       const std::string& name) const
TEST(qpp_QCircuit_get_measurement_count, SpecificMeasurementCount) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_measurement_depth() const
TEST(qpp_QCircuit_get_measurement_depth, TotalMeasurementDepth) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_measurement_depth(
///       const std::string& name) const
TEST(qpp_QCircuit_get_measurement_depth, SpecificMeasurementDepth) {}
/******************************************************************************/
/// BEGIN std::string qpp::QCircuit::get_name() const
TEST(qpp_QCircuit_get_name, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_nc() const noexcept
TEST(qpp_QCircuit_get_nc, AllTests) {}
/******************************************************************************/
/// BEGIN std::vector<idx> qpp::QCircuit::get_non_measured() const
TEST(qpp_QCircuit_get_non_measured, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_nop_count() const
TEST(qpp_QCircuit_get_nop_count, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_nq() const noexcept
TEST(qpp_QCircuit_get_nq, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::QCircuit::get_step_count() const noexcept
TEST(qpp_QCircuit_get_step_count, AllTests) {}
/******************************************************************************/
/// BEGIN bool qpp::QCircuit::is_clean(idx i) const
TEST(qpp_QCircuit_is_clean, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::kron(const QCircuit& qc)
TEST(qpp_QCircuit_kron, AllTests) {}
/******************************************************************************/
/// BEGIN friend QCircuit kron(const QCircuit& qc1, const QCircuit& qc2)
TEST(qpp_friend_QCircuit_kron, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::measureV(const cmat& V,
///       const std::vector<idx>& target, idx c_reg, bool destructive = true,
///       std::string name = {})
TEST(qpp_QCircuit_measureV, JointMultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::measureV(const cmat& V, idx target,
///       idx c_reg, bool destructive = true, std::string name = {})
TEST(qpp_QCircuit_measureV, SingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::measureZ(const std::vector<idx>& target,
///       idx c_reg, bool destructive = true, std::string name = {})
TEST(qpp_QCircuit_measureZ, MultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::measureZ(idx target, idx c_reg,
///       bool destructive = true, std::string name = {})
TEST(qpp_QCircuit_measureZ, SingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::nop()
TEST(qpp_QCircuit_nop, AllTests) {}
/******************************************************************************/
/// BEGIN bool qpp::QCircuit::operator!=(const QCircuit& rhs) const noexcept
TEST(qpp_QCircuit_operator_noneq, AllTests) {}
/******************************************************************************/
/// BEGIN bool qpp::QCircuit::operator==(const QCircuit& rhs) const noexcept
TEST(qpp_QCircuit_operator_eq, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::QFT(bool swap = true)
TEST(qpp_QCircuit_QFT, AllQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::QFT(const std::vector<idx>& target,
///       bool swap = true)
TEST(qpp_QCircuit_QFT, SpecificQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::remove_clean(idx target)
TEST(qpp_QCircuit_remove_clean, SingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::remove_clean(std::vector<idx> target)
TEST(qpp_QCircuit_remove_clean, MultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::replicate(idx n)
TEST(qpp_QCircuit_replicate, AllTests) {}
/******************************************************************************/
/// BEGIN friend QCircuit replicate(QCircuit qc, idx n)
TEST(qpp_friend_QCircuit_replicate, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::reset(const std::vector<idx>& target,
///       std::string name = {})
TEST(qpp_QCircuit_reset, MultipleTargets) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::reset(idx target, std::string name = {}
TEST(qpp_QCircuit_reset, SingleTarget) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::set_name(const std::string& name)
TEST(qpp_QCircuit_set_name, AllTests) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::TFQ(bool swap = true)
TEST(qpp_QCircuit_TFQ, AllQudits) {}
/******************************************************************************/
/// BEGIN QCircuit& qpp::QCircuit::TFQ(const std::vector<idx>& target,
///       bool swap = true)
TEST(qpp_QCircuit_TFQ, SpecificQudits) {}
/******************************************************************************/
/// BEGIN std::string qpp::QCircuit::to_JSON(
///       bool enclosed_in_curly_brackets = true) const override
TEST(qpp_QCircuit_to_JSON, AllTests) {}
/******************************************************************************/
