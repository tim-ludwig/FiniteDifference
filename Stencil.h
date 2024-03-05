#ifndef STENCIL_H
#define STENCIL_H

#include <cstddef>
#include <functional>
#include <Eigen/Dense>

class Stencil {
public:
    using Generator = std::function<Stencil(double, double)>;

private:
    int _rowOffset;
    int _colOffset;
    Eigen::MatrixXd _coeff;
    std::function<double(double, double)> _rhs;

public:
    Stencil(int rowOffset, int colOffset, Eigen::MatrixXd coeff, std::function<double(double, double)> rhs)
    : _rowOffset{rowOffset},
      _colOffset{colOffset},
      _coeff{coeff},
      _rhs{rhs}
    {

    }

    int rows() const { return _coeff.rows(); }
    int cols() const { return _coeff.cols(); }
    
    int rowOffset() const { return _rowOffset; }
    int colOffset() const { return _colOffset; }

    double operator()(int row, int col) const { return _coeff(row, col); }

    int nonZeros() const { return _coeff.nonZeros(); }

    double rhs(double x, double y) const { return _rhs(x, y); }
};

#endif//STENCIL_H