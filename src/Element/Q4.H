#ifndef ELEMENT_Q4_H
#define ELEMENT_Q4_H

#include "Element.H"

namespace Element
{
class Q4 : public Element<2,4,4>
{
public:
    using Element<2,4,4>::Element;

    std::array<double,4> phi(Set::Vector Y) override
    {
        std::array<double,4> ret;

        double xi = Y(0), eta = Y(1);
        // Q4 bilinear shape functions: N_i = (1/4)(1 +/- xi)(1 +/- eta)
        // Node 0 (-1,-1): N0 = (1/4)(1-xi)(1-eta)
        ret[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
        // Node 1 (1,-1): N1 = (1/4)(1+xi)(1-eta)
        ret[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
        // Node 2 (1,1): N2 = (1/4)(1+xi)(1+eta)
        ret[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
        // Node 3 (-1,1): N3 = (1/4)(1-xi)(1+eta)
        ret[3] = 0.25 * (1.0 - xi) * (1.0 + eta);

        return ret;
    }

    std::array<Set::Vector,4> Dphi(Set::Vector Y) override
    {
        std::array<Set::Vector,4> Dphi;
        for (int i = 0; i < 4; i++) Dphi[i] = Set::Vector::Zero();

        double xi = Y(0), eta = Y(1);
        // Derivatives: dN/dxi, dN/deta
        // N0 = 0.25(1-xi)(1-eta): dN0/dxi = -0.25(1-eta), dN0/deta = -0.25(1-xi)
        Dphi[0](0) = -0.25 * (1.0 - eta);
        Dphi[0](1) = -0.25 * (1.0 - xi);
        // N1 = 0.25(1+xi)(1-eta): dN1/dxi = 0.25(1-eta), dN1/deta = -0.25(1+xi)
        Dphi[1](0) =  0.25 * (1.0 - eta);
        Dphi[1](1) = -0.25 * (1.0 + xi);
        // N2 = 0.25(1+xi)(1+eta): dN2/dxi = 0.25(1+eta), dN2/deta = 0.25(1+xi)
        Dphi[2](0) =  0.25 * (1.0 + eta);
        Dphi[2](1) =  0.25 * (1.0 + xi);
        // N3 = 0.25(1-xi)(1+eta): dN3/dxi = -0.25(1+eta), dN3/deta = 0.25(1-xi)
        Dphi[3](0) = -0.25 * (1.0 + eta);
        Dphi[3](1) =  0.25 * (1.0 - xi);

        return Dphi;
    }

public:
    virtual const std::array<Set::Vector,4> & getY0() {return Y0;}

private:
    //
    // Hard-coded values common to all elements
    //

    const std::array<Set::Vector,4> Y0 = {
        Set::Vector(-1.0,-1.0),
        Set::Vector(1.0,-1.0),
        Set::Vector(1.0,1.0),
        Set::Vector(-1.0,1.0)
    };

};

}
#endif
