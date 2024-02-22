#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <cstddef>
#include <functional>

class Boundary {
private:
    double _xmin;
    double _xmax;
    double _ymin;
    double _ymax;
    std::function<double(double, double)> _boundaryFunction;

public:
    Boundary(double xmin, double xmax, double ymin, double ymax, std::function<double(double, double)> boundaryFunction)
    : _xmin{xmin},
      _xmax{xmax},
      _ymin{ymin},
      _ymax{ymax},
      _boundaryFunction{boundaryFunction}
    {
    
    }

    double xmin() const { return _xmin; }
    double xmax() const { return _xmax; }
    double ymin() const { return _ymin; }
    double ymax() const { return _ymax; }
    
    double operator()(double x, double y) const { return _boundaryFunction(x, y); }
};

#endif//BOUNDARY_H