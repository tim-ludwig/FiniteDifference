#include <iostream>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>

#include "FiniteDifference.h"

void solveProblem(FiniteDifferenceProblem& problem, double xStep, double yStep, std::string fileName = "result.csv") {
    std::cout << "Step size (" << xStep << ", " << yStep << ")" << std::endl;
    std::cout << "Preparing system of linear equations." << std::endl;

    FiniteDifferenceSystem system = problem.prepareSystem(xStep, yStep);

    std::cout << "Grid is (" << system.rows() << " * " << system.cols() << ") with " << system.unknowns() << " unknowns." << std::endl;
    std::cout << "solving..." << std::endl;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    Eigen::MatrixXd solution = system.solve();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Took " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " seconds." << std::endl;

    std::cout << "Writing result to 'simulation/" << fileName << "'" << std::endl;

    Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n");
    std::ofstream file("simulation/" + fileName);
    if (!file.is_open()) {
        std::cerr << "Opening file for writing failed." << std::endl;
    } else {
        file << solution.format(CSVFormat);
        file.close();
        std::cout << "Done." << std::endl;
    }
}

FiniteDifferenceProblem generateHeatFlowProblem() {
    auto stencilGen = [](double xStep, double yStep) {
        double alpha = 0.25; // parameter of heat flow equation
        double r = alpha * yStep / (xStep * xStep);

        return Stencil{1, 1,
            Eigen::MatrixXd{
                {0,          1, 0},
                {r, -1 - 2 * r, r}
            }, [](double x, double y) { return 0; }
        };
    };
    auto initial = std::optional<Stencil::Generator>{
        [](double xStep, double yStep) {
            return Stencil{0, 0,
                Eigen::MatrixXd{
                    { 1 }
                }, [](double x, double y) { return sin(2 * M_PI * x); }
            };
        }
    };
    auto left = std::optional<Stencil::Generator>{
        [](double xStep, double yStep) {
            return Stencil{0, 0,
                Eigen::MatrixXd{
                    { -1, 1 }
                }, [](double x, double y) { return 0; }
            };
        }
    };
    auto right = std::optional<Stencil::Generator>{
        [](double xStep, double yStep) {
            return Stencil{0, 1,
                Eigen::MatrixXd{
                    { -1, 1 }
                }, [](double x, double y) { return 0; }
            };
        }
    };

    BoundaryConditions b{left, right, initial, std::optional<Stencil::Generator>{}};
    return FiniteDifferenceProblem{Domain{0, 1, 0, 1}, b, stencilGen};
}

int main() {
    std::cout << "Heatflow" << std::endl;
    FiniteDifferenceProblem heatFlowProblem = generateHeatFlowProblem();
    solveProblem(heatFlowProblem, 0.001, 0.01, "heatflow.csv");
    std::cout << std::endl;
}