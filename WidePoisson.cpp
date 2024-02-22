#include <Eigen/Sparse>
#include <functional>
#include <iostream>

#include "Boundary.h"

typedef struct {
    Eigen::SparseMatrix<double, Eigen::RowMajor> a;
    Eigen::VectorXd b;
} System;

Eigen::Matrix<double, 5, 5> stencil{
    {0,   0,   1,   0, 0},
    {0,   0, -16,   0, 0},
    {1, -16,  60, -16, 1},
    {0,   0, -16,   0, 0},
    {0,   0,   1,   0, 0},
};

System generateSystem(std::function<double(double,double)> f, double h) {
    int rows = 1 / h + 1;
    int cols = 1 / h + 1;

    int firstInnerRow = 1;
    int firstInnerCol = 1;

    int lastInnerRow = rows - 2;
    int lastInnerCol = cols - 2;

    int numInnerRows = lastInnerRow - firstInnerRow + 1;
    int numInnerCols = lastInnerCol - firstInnerCol + 1;

    int n = numInnerRows * numInnerCols;

    Eigen::SparseMatrix<double, Eigen::RowMajor> a{n, n};
    a.reserve(Eigen::VectorXi::Constant(n, 9));
    Eigen::VectorXd b{n};

    #define NODE_IDX(row, col) (((row) - firstInnerRow) * numInnerCols + (col) - firstInnerCol)
    #define IS_INNER_NODE(row, col) ((firstInnerRow <= (row) && (row) <= lastInnerRow) && (firstInnerCol <= (col) && (col) <= lastInnerCol))
    #define X_COORD(col) (h * (col))
    #define Y_COORD(row) (h * (row))

    // for each inner node
    for (int row = firstInnerRow; row <= lastInnerRow; ++row) {
        for (int col = firstInnerCol; col <= lastInnerCol; ++col) {
            b(NODE_IDX(row, col)) = f(X_COORD(col), Y_COORD(row));

            // for each neighbour in the stencil
            for (int i = 0, otherRow = row - 2; i < 5; ++i, ++otherRow) {
                for (int j = 0, otherCol = col - 2; j < 5; ++j, ++otherCol) {
                    // eliminate boundary conditions
                    if (IS_INNER_NODE(otherRow, otherCol)) {
                        if (stencil(i, j) != 0)
                            a.insert(NODE_IDX(row, col), NODE_IDX(otherRow, otherCol)) = stencil(i, j);
                    }
                }
            }

            // change discretization for nodes close to borders
            if (row == firstInnerRow) {
                a.coeffRef(NODE_IDX(row, col), NODE_IDX(row, col)) += 2;
                a.coeffRef(NODE_IDX(row, col), NODE_IDX(3, col)) -= 1;
                b(NODE_IDX(row, col)) += (f(X_COORD(0), Y_COORD(col)) + f(X_COORD(1), Y_COORD(col)) + f(X_COORD(2), Y_COORD(col))) / 9;
            } else if (row == lastInnerRow) {

            }

            if (col == firstInnerCol) {
                a.coeffRef(NODE_IDX(row, col), NODE_IDX(row, col)) += 2;
                a.coeffRef(NODE_IDX(row, col), NODE_IDX(row, 3)) -= 1;
                b(NODE_IDX(row, col)) += (f(X_COORD(row), Y_COORD(0)) + f(X_COORD(row), Y_COORD(1)) + f(X_COORD(row), Y_COORD(2))) / 9;
            } else if (col == lastInnerCol) {

            }

            // apply common factor to rhs vector
            b(NODE_IDX(row, col)) *= 12 * h * h;
        }
    }

    a.makeCompressed();

    return System{a, b};

    #undef NODE_IDX
    #undef IS_INNER_NODE
    #undef X_COORD
    #undef Y_COORD
}

int main() {
    System s = generateSystem([](double x, double y) {
        return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
    }, 1.0 / 6);

    Eigen::IOFormat cleanFormat(0, 0, " ", "\n");
    std::cout << s.a.toDense().format(cleanFormat) << std::endl;
}