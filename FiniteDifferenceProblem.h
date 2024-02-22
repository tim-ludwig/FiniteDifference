#ifndef FINITE_DIFFERENCE_PROBLEM_H
#define FINITE_DIFFERENCE_PROBLEM_H

#include <functional>
#include <Eigen/Sparse>
#include "Boundary.h"
#include "Stencil.h"
#include "FiniteDifferenceSystem.h"

class FiniteDifferenceProblem {
private:
    Boundary _boundary;
    std::function<Stencil(double, double)> _stencilGenerator;
    std::function<double(double, double)> _f;

public:
    FiniteDifferenceProblem(Boundary boundary, std::function<Stencil(double, double)> stencilGenerator, std::function<double(double, double)> f=[](double x, double y) { return 0; })
    : _boundary{boundary},
      _stencilGenerator{stencilGenerator},
      _f{f}
    {
    
    }

    FiniteDifferenceSystem prepareSystem(double xStep, double yStep) const {
        Stencil stencil = _stencilGenerator(xStep, yStep);
        
        int rows = (_boundary.ymax() - _boundary.ymin()) / yStep + 1;
        int cols = (_boundary.xmax() - _boundary.xmin()) / xStep + 1;

        int firstInnerRow = stencil.rowOffset() > 0 ? 1 : 0;
        int firstInnerCol = stencil.colOffset() > 0 ? 1 : 0;

        int lastInnerRow = stencil.rowOffset() < stencil.rows() - 1 ? rows - 2 : rows - 1;
        int lastInnerCol = stencil.colOffset() < stencil.cols() - 1 ? cols - 2 : cols - 1;

        int numInnerRows = lastInnerRow - firstInnerRow + 1;
        int numInnerCols = lastInnerCol - firstInnerCol + 1;

        int n = numInnerRows * numInnerCols;

        Eigen::SparseMatrix<double, Eigen::RowMajor> a{n, n};
        a.reserve(Eigen::VectorXi::Constant(n, stencil.nonZeros()));
        Eigen::VectorXd b{n};

        #define NODE_IDX(row, col) (((row) - firstInnerRow) * numInnerCols + (col) - firstInnerCol)
        #define IS_INNER_NODE(row, col) ((firstInnerRow <= (row) && (row) <= lastInnerRow) && (firstInnerCol <= (col) && (col) <= lastInnerCol))
        #define X_COORD(col) (_boundary.xmin() + xStep * (col))
        #define Y_COORD(row) (_boundary.ymin() + yStep * (row))

        // for each inner node
        for (int row = firstInnerRow; row <= lastInnerRow; ++row) {
            for (int col = firstInnerCol; col <= lastInnerCol; ++col) {
                b(NODE_IDX(row, col)) = stencil.commonCoeff() * _f(X_COORD(col), Y_COORD(row));

                // for each neighbour in the stencil
                for (int i = 0, otherRow = row - stencil.rowOffset(); i < stencil.rows(); ++i, ++otherRow) {
                    for (int j = 0, otherCol = col - stencil.colOffset(); j < stencil.cols(); ++j, ++otherCol) {
                        // eliminate boundary conditions
                        if (!IS_INNER_NODE(otherRow, otherCol)) {
                            b(NODE_IDX(row, col)) -= stencil(i, j) * _boundary(X_COORD(otherCol), Y_COORD(otherRow));
                        } else {
                            if (stencil(i, j) != 0)
                                a.insert(NODE_IDX(row, col), NODE_IDX(otherRow, otherCol)) = stencil(i, j);
                        }
                    }
                }
            }
        }

        a.makeCompressed();

        return FiniteDifferenceSystem{_boundary, xStep, yStep, rows, cols, firstInnerRow, lastInnerRow, firstInnerCol, lastInnerCol, a, b};

        #undef NODE_IDX
        #undef IS_INNER_NODE
        #undef X_COORD
        #undef Y_COORD
    }
};

#endif//FINITE_DIFFERENCE_PROBLEM_H