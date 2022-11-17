#include <vector>
#include <tuple>

std::tuple<long, long, int> runCusparseCublasCG(std::vector<int> &I, std::vector<int> &J, std::vector<float> &val, std::vector<float> &rhs, std::vector<float> &x, int M, int N, int nz, float tol, int max_iter);