// Standalone example, assumes quantum++ is installed in a system-wide visible
// directory

#include <iostream>

#include <qpp/qpp.h>

int main() {
    using namespace qpp;

    QCircuit qc{1, 1, 2, "coin flip"};
    qc.gate(gt.H, 0);
    qc.measure_all();

    std::cout << qc << "\n\n" << qc.get_resources() << "\n\n";
    std::cout << QEngine{qc}.execute(100) << "\n";

    cmat A = st.z0;

    std::cout << A.size() << '\n';
    std::cout << A(0, 0) << '\n';
    std::cout << A(1, 0) << '\n';

    auto res = n2multiidx(0, std::vector<idx>{1});
    std::cout << disp(res, " ") << '\n';

    idx n_rows = internal::get_num_subsys(static_cast<idx>(A.rows()), 2);
    std::cout << n_rows << '\n';

    idx n_cols = internal::get_num_subsys(static_cast<idx>(A.cols()), 2);
    std::cout << n_cols << '\n';

    std::vector<idx> col_dims(n_cols, 2);
    std::cout << disp(col_dims, " ") << '\n';
}
