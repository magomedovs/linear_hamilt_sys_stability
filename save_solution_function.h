#ifndef SAVE_SOLUTION
#define SAVE_SOLUTION

#include <vector>
#include <string>
#include <fstream>

//template <typename Container_T, typename Var_T>
//void SaveSolutionIntoFile(const std::vector<Container_T>& x_vec, const std::vector<Var_T>& t_vec, std::string filename);

template <typename Container_T, typename Var_T>
void SaveSolutionIntoFile(const std::vector<Container_T>& x_vec, const std::vector<Var_T>& t_vec, std::string filename="ode_output.txt") {
    if (!x_vec.empty() && !t_vec.empty()) {
        std::ofstream output(filename);
        size_t dim = x_vec[0].size();
        for (size_t i = 0; i < t_vec.size(); ++i) {
            output << t_vec[i];
            for (size_t j = 0; j < dim; ++j) {
                output << " " << x_vec[i][j];
            }
            output << "\n";
        }
        output.close();
    } else {
        std::cout << "Solution is not saved into file.\n";
    }
}

#endif
