#ifndef FINITE_DIFFERENCE_SYSTEM_H
#define FINITE_DIFFERENCE_SYSTEM_H

#include <Eigen/Sparse>
#include "Domain.h"

class FiniteDifferenceSystem {
private:
    int _rows;
    int _cols;

public:
    Eigen::SparseMatrix<double, Eigen::RowMajor> _a;
    Eigen::VectorXd _b;

public:
    FiniteDifferenceSystem(int rows, int cols, Eigen::SparseMatrix<double, Eigen::RowMajor> a, Eigen::VectorXd b)
    : _rows{rows},
      _cols{cols},
      _a{a},
      _b{b}
    {

    }
    
    int rows() const { return _rows; }
    int cols() const { return _cols; }
    int unknowns() const { return _a.rows(); }

    Eigen::MatrixXd solve() const {
        #define NODE_IDX(row, col) ((row) * _cols + (col))

        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::NaturalOrdering<int>> solver;
        solver.compute(_a);
        if (solver.info() != Eigen::Success) { std::cerr << "solver::compute failed" << std::endl; }
        Eigen::VectorXd x = solver.solve(_b);
        if (solver.info() != Eigen::Success) { std::cerr << "solver::solve failed" << std::endl; }
        Eigen::MatrixXd solution{_rows, _cols};

        for (int row = 0, nodeIdx = 0; row < _rows; ++row) {
            for (int col = 0; col < _cols; ++col, ++nodeIdx) {
                solution(row, col) = x(NODE_IDX(row, col));
            }
        }

        return solution;

        #undef NODE_IDX
    }
};

#endif//FINITE_DIFFERENCE_SYSTEM_H