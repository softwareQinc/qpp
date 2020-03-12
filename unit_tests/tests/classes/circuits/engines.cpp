#include <vector>
#include "gtest/gtest.h"
#include "qpp.h"

using namespace qpp;

// Unit testing "classes/circuits/engines.h"

/******************************************************************************/
/// BEGIN QEngine& qpp::QEngine::execute(const QCircuit::iterator& it)
TEST(qpp_QEngine_execute, IteratorStep) {}
/******************************************************************************/
/// BEGIN virtual QEngine& qpp::QEngine::execute(
///       const QCircuit::iterator::value_type& elem)
TEST(qpp_QEngine_execute, CircuitStep) {}
/******************************************************************************/
/// BEGIN QEngine& qpp::QEngine::execute(idx reps = 1, bool clear_stats = true)
TEST(qpp_QEngine_execute, AllCircuitWithRepetitions) {}
/******************************************************************************/
/// BEGIN const QEngine& qpp::QEngine::get_circuit() const& noexcept AND
///       QCircuit get_circuit() const&& noexcept
TEST(qpp_QEngine_get_circuit, AllTests) {}
/******************************************************************************/
/// BEGIN idx qpp::QEngine::get_dit(idx i) const
TEST(qpp_QEngine_get_dit, AllTests) {}
/******************************************************************************/
/// BEGIN std::vector<idx> qpp::QEngine::get_dits() const
TEST(qpp_QEngine_get_dits, AllTests) {}
/******************************************************************************/
/// BEGIN std::vector<idx> qpp::QEngine::get_measured() const
TEST(qpp_QEngine_get_measured, AllQudits) {}
/******************************************************************************/
/// BEGIN bool qpp::QEngine::get_measured(idx i) const
TEST(qpp_QEngine_get_measured, SpecificQudit) {}
/******************************************************************************/
/// BEGIN std::vector<idx> qpp::QEngine::get_non_measured() const
TEST(qpp_QEngine_get_non_measured, AllTests) {}
/******************************************************************************/
/// BEGIN std::vector<double> qpp::QEngine::get_probs() const
TEST(qpp_QEngine_get_probs, AllTests) {}
/******************************************************************************/
/// BEGIN ket qpp::QEngine::get_psi() const
TEST(qpp_QEngine_get_psi, AllTests) {}
/******************************************************************************/
/// BEGIN std::map<std::string, idx, internal::EqualSameSizeStringDits>
///       qpp::QEngine::get_stats() const
TEST(qpp_QEngine_get_stats, AllTests) {}
/******************************************************************************/
/// BEGIN QEngine& qpp::QEngine::reset()
TEST(qpp_QEngine_reset, AllTests) {}
/******************************************************************************/
/// BEGIN QEngine& qpp::QEngine::reset_stats()
TEST(qpp_QEngine_reset_stats, AllTests) {}
/******************************************************************************/
/// BEGIN QEngine& qpp::QEngine::set_dit(idx i, idx value)
TEST(qpp_QEngine_set_dit, AllTests) {}
/******************************************************************************/
/// BEGIN QEngine& qpp::QEngine::set_psi(const ket& psi)
TEST(qpp_QEngine_set_psi, AllTests) {}
/******************************************************************************/
/// BEGIN std::string qpp::QEngine::to_JSON(
///       bool enclosed_in_curly_brackets = true) const override
TEST(qpp_QEngine_to_JSON, AllTests) {}
/******************************************************************************/

/******************************************************************************/
/// BEGIN QNoisyEngine& qpp::QNoisyEngine::execute(
///       const QCircuit::iterator::value_type& elem) override
TEST(qpp_QNoisyEngine_execute, CircuitStep) {}
/******************************************************************************/
/// BEGIN std::vector<std::vector<idx>>
///       qpp::QNoisyEngine::get_noise_results() const
TEST(qpp_QNoisyEngine_get_noise_results, AllTests) {}
/******************************************************************************/
