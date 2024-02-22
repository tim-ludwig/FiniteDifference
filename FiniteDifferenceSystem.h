#ifndef FINITE_DIFFERENCE_SYSTEM_H
#define FINITE_DIFFERENCE_SYSTEM_H

#include <Eigen/Sparse>
#include "Boundary.h"

class FiniteDifferenceSystem {
private:
    Boundary _boundary;
    double _xStep;
    double _yStep;
    int _rows;
    int _cols;
    int _firstInnerRow;
    int _lastInnerRow;
    int _firstInnerCol;
    int _lastInnerCol;

public:
    Eigen::SparseMatrix<double, Eigen::RowMajor> _a;
    Eigen::VectorXd _b;

public:
    FiniteDifferenceSystem(Boundary boundary, double xStep, double yStep, int rows, int cols, int firstInnerRow, int lastInnerRow, int firstInnerCol, int lastInnerCol, Eigen::SparseMatrix<double, Eigen::RowMajor> a, Eigen::VectorXd b)
    : _boundary{boundary},
      _xStep{xStep},
      _yStep{yStep},
      _rows{rows},
      _cols{cols},
      _firstInnerRow{firstInnerRow},
      _lastInnerRow{lastInnerRow},
      _firstInnerCol{firstInnerCol},
      _lastInnerCol{lastInnerCol},
      _a{a},
      _b{b}
    {

    }
    
    int rows() const { return _rows; }
    int cols() const { return _cols; }
    int unknowns() const { return _a.rows(); }

    Eigen::MatrixXd solve() const {
        #define NODE_IDX(row, col) ((((row) - _firstInnerRow) * (_lastInnerCol - _firstInnerCol + 1)) + (col) - _firstInnerCol)
        #define IS_INNER_NODE(row, col) ((_firstInnerRow <= (row) && (row) <= _lastInnerRow) && (_firstInnerCol <= (col) && (col) <= _lastInnerCol))
        #define X_COORD(col) (_boundary.xmin() + _xStep * (col))
        #define Y_COORD(row) (_boundary.ymin() + _yStep * (row))

        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::NaturalOrdering<int>> solver;
        solver.compute(_a);
        if (solver.info() != Eigen::Success) { std::cerr << "solver::compute failed" << std::endl; }
        Eigen::VectorXd x = solver.solve(_b);
        if (solver.info() != Eigen::Success) { std::cerr << "solver::solve failed" << std::endl; }
        Eigen::MatrixXd solution{_rows, _cols};

        for (int row = 0, nodeIdx = 0; row < _rows; ++row) {
            for (int col = 0; col < _cols; ++col, ++nodeIdx) {
                // add boundary values to solution matrix
                if (!IS_INNER_NODE(row, col))  {
                    solution(row, col) = _boundary(X_COORD(col), Y_COORD(row));
                } else {
                    solution(row, col) = x(NODE_IDX(row, col));
                }
            }
        }

        return solution;

        #undef NODE_IDX
        #undef IS_INNER_NODE
        #undef X_COORD
        #undef Y_COORD
    }
};

#endif//FINITE_DIFFERENCE_SYSTEM_H