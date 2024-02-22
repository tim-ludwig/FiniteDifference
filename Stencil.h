#ifndef STENCIL_H
#define STENCIL_H

#include <cstddef>
#include <Eigen/Dense>

class Stencil {
private:
    int _rowOffset;
    int _colOffset;
    double _commonCoeff;
    Eigen::MatrixXd _coeff;

public:
    Stencil(int rowOffset, int colOffset, Eigen::MatrixXd coeff, double commonCoeff=1)
    : _rowOffset{rowOffset},
      _colOffset{colOffset},
      _coeff{coeff},
      _commonCoeff{commonCoeff}
    {

    }

    int rows() const { return _coeff.rows(); }
    int cols() const { return _coeff.cols(); }
    int rowOffset() const { return _rowOffset; }
    int colOffset() const { return _colOffset; }
    double commonCoeff() const { return _commonCoeff; }

    double operator()(int row, int col) const { return _coeff(row, col); }

    int nonZeros() const { return _coeff.nonZeros(); }
};

#endif//STENCIL_H