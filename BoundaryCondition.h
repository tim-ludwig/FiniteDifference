#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include <cstddef>
#include <functional>
#include <optional>

#include "Stencil.h"

class BoundaryConditions {
private:
    std::optional<Stencil::Generator> _left;
    std::optional<Stencil::Generator> _right;
    std::optional<Stencil::Generator> _bottom;
    std::optional<Stencil::Generator> _top;

public:
    BoundaryConditions(std::optional<Stencil::Generator> left, std::optional<Stencil::Generator> right, std::optional<Stencil::Generator> bottom, std::optional<Stencil::Generator> top)
    : _left{left},
      _right{right},
      _bottom{bottom},
      _top{top}
    {
    
    }
    
    std::optional<Stencil> generateLeft(double xStep, double yStep) const { return _left ? (*_left)(xStep, yStep) : std::optional<Stencil>{}; }
    std::optional<Stencil> generateRight(double xStep, double yStep) const { return _right ? (*_right)(xStep, yStep) : std::optional<Stencil>{}; }
    std::optional<Stencil> generateBottom(double xStep, double yStep) const { return _bottom ? (*_bottom)(xStep, yStep) : std::optional<Stencil>{}; }
    std::optional<Stencil> generateTop(double xStep, double yStep) const { return _top ? (*_top)(xStep, yStep) : std::optional<Stencil>{}; }
};

#endif//BOUNDARY_CONDITIONS_H