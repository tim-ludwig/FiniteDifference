#ifndef FINITE_DIFFERENCE_PROBLEM_H
#define FINITE_DIFFERENCE_PROBLEM_H

#include <functional>
#include <Eigen/Sparse>
#include "BoundaryCondition.h"
#include "Stencil.h"
#include "FiniteDifferenceSystem.h"

class FiniteDifferenceProblem {
private:
    Domain _domain;
    BoundaryConditions _boundaryConditions;
    Stencil::Generator _innerGenerator;

public:
    FiniteDifferenceProblem(Domain domain, BoundaryConditions boundaryConditions, Stencil::Generator innerGenerator)
    : _domain{domain},
      _boundaryConditions{boundaryConditions},
      _innerGenerator{innerGenerator}
    {
    
    }

    FiniteDifferenceSystem prepareSystem(double xStep, double yStep) const {
        Stencil innerStencil = _innerGenerator(xStep, yStep);
        std::optional<Stencil> leftStencil = _boundaryConditions.generateLeft(xStep, yStep);
        std::optional<Stencil> rightStencil = _boundaryConditions.generateRight(xStep, yStep);
        std::optional<Stencil> topStencil = _boundaryConditions.generateTop(xStep, yStep);
        std::optional<Stencil> bottomStencil = _boundaryConditions.generateBottom(xStep, yStep);
        
        int rows = (_domain.ymax - _domain.ymin) / yStep + 1;
        int cols = (_domain.xmax - _domain.xmin) / xStep + 1;

        int n = rows * cols;

        Eigen::SparseMatrix<double, Eigen::RowMajor> a{n, n};
        a.reserve(Eigen::VectorXi::Constant(n, innerStencil.nonZeros()));
        Eigen::VectorXd b{n};

        #define NODE_IDX(row, col) ((row) * cols + (col))
        #define X_COORD(col) (_domain.xmin + (col) * xStep)
        #define Y_COORD(row) (_domain.ymin + (row) * yStep)

        // for each node
        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                // TODO select appropriate stencil form problem statement / boundaries
                Stencil* stencil = &innerStencil;

                if (row == 0 && bottomStencil) {
                    stencil = &*bottomStencil;
                } else if (col == 0 && leftStencil) {
                    stencil = &*leftStencil;
                } else if (col == cols - 1 && rightStencil) {
                    stencil = &*rightStencil;
                } else if (row == rows - 1 && topStencil) {
                    stencil = &*topStencil;
                }

                // apply the stencil
                b(NODE_IDX(row, col)) = stencil->rhs(X_COORD(col), Y_COORD(row));
                for (int i = 0, otherRow = row - stencil->rowOffset(); i < stencil->rows(); ++i, ++otherRow) {
                    for (int j = 0, otherCol = col - stencil->colOffset(); j < stencil->cols(); ++j, ++otherCol) {
                        if ((*stencil)(i, j) != 0)
                            a.insert(NODE_IDX(row, col), NODE_IDX(otherRow, otherCol)) = (*stencil)(i, j);
                    }
                }
            }
        }

        a.makeCompressed();

        return FiniteDifferenceSystem{rows, cols, a, b};

        #undef NODE_IDX
        #undef X_COORD
        #undef Y_COORD
    }
};

#endif//FINITE_DIFFERENCE_PROBLEM_H