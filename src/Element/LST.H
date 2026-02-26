#ifndef ELEMENT_LST_H
#define ELEMENT_LST_H

#include "Element.H"

namespace Element
{
class LST : public Element<2,6,4>
{
public:
    using Element<2,6,4>::Element;

    std::array<double,6> phi(Set::Vector Y) override
    {
        std::array<double,6> ret;

        // Area coordinates for right triangle: corners at (0,0), (1,0), (0,1)
        double x = Y(0), y = Y(1);
        double L1 = 1.0 - x - y;
        double L2 = x;
        double L3 = y;

        // LST shape functions: corners N_i = L_i(2L_i-1), mids N = 4*L_i*L_j
        // Node 0 (0,0): corner 1 -> N0 = L1(2L1-1)
        ret[0] = L1 * (2.0*L1 - 1.0);
        // Node 1 (0.5,0): mid 1-2 -> N1 = 4*L1*L2
        ret[1] = 4.0 * L1 * L2;
        // Node 2 (1,0): corner 2 -> N2 = L2(2L2-1)
        ret[2] = L2 * (2.0*L2 - 1.0);
        // Node 3 (0.5,0.5): mid 2-3 -> N3 = 4*L2*L3
        ret[3] = 4.0 * L2 * L3;
        // Node 4 (0,1): corner 3 -> N4 = L3(2L3-1)
        ret[4] = L3 * (2.0*L3 - 1.0);
        // Node 5 (0,0.5): mid 3-1 -> N5 = 4*L3*L1
        ret[5] = 4.0 * L3 * L1;

        return ret;
    }

    std::array<Set::Vector,6> Dphi(Set::Vector Y) override
    {
        std::array<Set::Vector,6> ret;
        for (int i = 0; i < 6; i++) ret[i] = Set::Vector::Zero();

        double x = Y(0), y = Y(1);
        double L1 = 1.0 - x - y;
        double L2 = x;
        double L3 = y;

        // dL1/dx=-1, dL1/dy=-1; dL2/dx=1, dL2/dy=0; dL3/dx=0, dL3/dy=1
        // dN0/dx = (4L1-1)*dL1/dx = -(4L1-1), dN0/dy = -(4L1-1)
        ret[0](0) = -(4.0*L1 - 1.0);
        ret[0](1) = -(4.0*L1 - 1.0);
        // dN1/dx = 4(L1*dL2/dx + L2*dL1/dx) = 4(L1 - L2), dN1/dy = 4*L2*dL1/dy = -4*L2
        ret[1](0) = 4.0 * (L1 - L2);
        ret[1](1) = -4.0 * L2;
        // dN2/dx = (4L2-1), dN2/dy = 0
        ret[2](0) = 4.0*L2 - 1.0;
        ret[2](1) = 0.0;
        // dN3/dx = 4*L3, dN3/dy = 4*L2
        ret[3](0) = 4.0 * L3;
        ret[3](1) = 4.0 * L2;
        // dN4/dx = 0, dN4/dy = (4L3-1)
        ret[4](0) = 0.0;
        ret[4](1) = 4.0*L3 - 1.0;
        // dN5/dx = 4*(L3*dL1/dx + L1*dL3/dx) = -4*L3, dN5/dy = 4*(L1 + L3*dL1/dy) = 4*(L1 - L3)
        ret[5](0) = -4.0 * L3;
        ret[5](1) = 4.0 * (L1 - L3);

        return ret;
    }

public:
    virtual const std::array<Set::Vector,6> & getY0() {return Y0;}

private:
    //
    // Hard-coded values common to all elements
    //

    const std::array<Set::Vector,6> Y0 = {
        Set::Vector(0.0,0.0),
        Set::Vector(0.5,0.0),
        Set::Vector(1.0,0.0),
        Set::Vector(0.5,0.5),
        Set::Vector(0.0,1.0),
        Set::Vector(0.0,0.5)
    };

};

}
#endif
