#include <vector>
#include <tuple>

std::tuple<long, long, int> runCusparseCublasCG(std::vector<int> &I, std::vector<int> &J, std::vector<double> &val, std::vector<double> &rhs, std::vector<double> &x, int M, int N, int nz, double tol, int max_iter);