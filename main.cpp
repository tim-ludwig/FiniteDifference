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
        std::exit(1);
    }

    file << solution.format(CSVFormat);
    file.close();
    std::cout << "Done." << std::endl;
}

FiniteDifferenceProblem generateHeatFlowProblem() {
    auto g = [](double x, double y) {
        return sin(2 * M_PI * x);
    };
    auto stencilGen = [](double xStep, double yStep) {
        double alpha = 0.06; // parameter of heat flow equation
        double r = alpha * yStep / (xStep * xStep);

        return Stencil{1, 1,
            Eigen::MatrixXd{
                {0,          1, 0},
                {r, -1 - 2 * r, r}
            }
        };
    };
    Boundary b{0, 1, 0, 1, g};

    return FiniteDifferenceProblem{b, stencilGen};
}

FiniteDifferenceProblem generateHeatFlowCrankNicolsonProblem() {
    auto g = [](double x, double y) {
        return sin(2 * M_PI * x);
    };
    auto stencilGen = [](double xStep, double yStep) {
        double alpha = 0.06; // parameter of heat flow equation
        double r = alpha * yStep / (xStep * xStep);

        return Stencil{1, 1,
            Eigen::MatrixXd{
                {r,  2 - 2 * r, r},
                {r, -2 - 2 * r, r}
            }
        };
    };
    Boundary b{0, 1, 0, 1, g};

    return FiniteDifferenceProblem{b, stencilGen};
}

FiniteDifferenceProblem generateTransportProblem() {
    auto g = [](double x, double y) {
        if (y == 0) return sin(M_PI * x);
        return 0.0;
    };
    auto stencilGen = [](double xStep, double yStep) {
        double c = 0.5;
        return Stencil{1, 1,
            Eigen::MatrixXd{
                {         0, -1 / yStep            },
                {-c / xStep,  1 / yStep + c / xStep}
            }
        };
    };
    Boundary b{0, 1, 0, 1, g};

    return FiniteDifferenceProblem{b, stencilGen};
}

FiniteDifferenceProblem generatePoissonProblemFivePointStencil() {
    auto f = [](double x, double y) {
        return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
    };
    auto g = [](double x, double y) {
        return 0;
    };
    auto stencilGen = [](double xStep, double yStep) {
        return Stencil{1, 1,
            Eigen::MatrixXd{
                { 0, -1,  0},
                {-1,  4, -1},
                { 0, -1,  0}
            }, xStep * xStep
        };
    };
    Boundary b{0, 1, 0, 1, g};

    return FiniteDifferenceProblem{b, stencilGen, f};
}

FiniteDifferenceProblem generatePoissonProblemNinePointStencil() {
    auto f = [](double x, double y) {
        return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
    };
    auto g = [](double x, double y) {
        return 0;
    };
    auto stencilGen = [](double xStep, double yStep) {
        return Stencil{1, 1,
            Eigen::MatrixXd{
                {-1, -4, -1},
                {-4, 20, -4},
                {-1, -4, -1}
            }, 6 * xStep * xStep
        };
    };
    Boundary b{0, 1, 0, 1, g};

    return FiniteDifferenceProblem{b, stencilGen, f};
}

int main() {
    std::cout << "Heatflow" << std::endl;
    FiniteDifferenceProblem heatFlowProblem = generateHeatFlowProblem();
    solveProblem(heatFlowProblem, 0.002, 0.01, "heatflow.csv");
    std::cout << std::endl;

    std::cout << "Heatflow Crank Nicolson" << std::endl;
    FiniteDifferenceProblem heatFlowCrankNicolsonProblem = generateHeatFlowCrankNicolsonProblem();
    solveProblem(heatFlowCrankNicolsonProblem, 0.002, 0.01, "heatflowCN.csv");
    std::cout << std::endl;

    std::cout << "Poisson five point stencil" << std::endl;
    FiniteDifferenceProblem poisson5 = generatePoissonProblemFivePointStencil();
    solveProblem(poisson5, 0.05, 0.05, "poisson5.csv");
    std::cout << std::endl;

    std::cout << "Poisson nine point stencil" << std::endl;
    FiniteDifferenceProblem poisson9 = generatePoissonProblemNinePointStencil();
    solveProblem(poisson9, 0.05, 0.05, "poisson9.csv");
    std::cout << std::endl;
    
    std::cout << "Transport" << std::endl;
    FiniteDifferenceProblem transportProblem = generateTransportProblem();
    solveProblem(transportProblem, 0.002, 0.01, "transport.csv");
    std::cout << std::endl;
}