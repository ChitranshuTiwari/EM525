#ifndef ELEMENT_Q9_H
#define ELEMENT_Q9_H

#include "Set/Set.H"
#include "Element.H"

namespace Element
{
class Q9 : public Element<2,9,9>
{
public:
    using Element<2,9,9>::Element;

    std::array<double,9> phi(Set::Vector Y) override
    {
        std::array<double,9> ret;

        double xi = Y(0), eta = Y(1);
        // Q9 Lagrange: 1D basis at -1,0,1: N_{-1}(s)=s(s-1)/2, N_0(s)=1-s^2, N_1(s)=s(s+1)/2
        double Nm1_xi = 0.5 * xi * (xi - 1.0);
        double N0_xi  = 1.0 - xi*xi;
        double N1_xi  = 0.5 * xi * (xi + 1.0);
        double Nm1_eta = 0.5 * eta * (eta - 1.0);
        double N0_eta  = 1.0 - eta*eta;
        double N1_eta  = 0.5 * eta * (eta + 1.0);

        // Nodes: 0(-1,-1), 1(0,-1), 2(1,-1), 3(1,0), 4(1,1), 5(0,1), 6(-1,1), 7(-1,0), 8(0,0)
        ret[0] = Nm1_xi * Nm1_eta;  // (-1,-1)
        ret[1] = N0_xi  * Nm1_eta;  // (0,-1)
        ret[2] = N1_xi  * Nm1_eta;  // (1,-1)
        ret[3] = N1_xi  * N0_eta;   // (1,0)
        ret[4] = N1_xi  * N1_eta;   // (1,1)
        ret[5] = N0_xi  * N1_eta;   // (0,1)
        ret[6] = Nm1_xi * N1_eta;   // (-1,1)
        ret[7] = Nm1_xi * N0_eta;   // (-1,0)
        ret[8] = N0_xi  * N0_eta;   // (0,0)

        return ret;
    }

    std::array<Set::Vector,9> Dphi(Set::Vector Y) override
    {
        std::array<Set::Vector,9> Dret;
        for (int i = 0; i < 9; i++) Dret[i] = Set::Vector::Zero();

        double xi = Y(0), eta = Y(1);
        // 1D basis derivatives: dNm1/ds=(2s-1)/2, dN0/ds=-2s, dN1/ds=(2s+1)/2
        double dNm1_xi = 0.5 * (2.0*xi - 1.0);
        double dN0_xi  = -2.0 * xi;
        double dN1_xi  = 0.5 * (2.0*xi + 1.0);
        double dNm1_eta = 0.5 * (2.0*eta - 1.0);
        double dN0_eta  = -2.0 * eta;
        double dN1_eta  = 0.5 * (2.0*eta + 1.0);

        double Nm1_xi = 0.5 * xi * (xi - 1.0);
        double N0_xi  = 1.0 - xi*xi;
        double N1_xi  = 0.5 * xi * (xi + 1.0);
        double Nm1_eta = 0.5 * eta * (eta - 1.0);
        double N0_eta  = 1.0 - eta*eta;
        double N1_eta  = 0.5 * eta * (eta + 1.0);

        // dN/dxi = (dN_xi/dxi)*N_eta, dN/deta = N_xi*(dN_eta/deta)
        Dret[0](0) = dNm1_xi * Nm1_eta;  Dret[0](1) = Nm1_xi * dNm1_eta;
        Dret[1](0) = dN0_xi  * Nm1_eta;  Dret[1](1) = N0_xi  * dNm1_eta;
        Dret[2](0) = dN1_xi  * Nm1_eta;  Dret[2](1) = N1_xi  * dNm1_eta;
        Dret[3](0) = dN1_xi  * N0_eta;   Dret[3](1) = N1_xi  * dN0_eta;
        Dret[4](0) = dN1_xi  * N1_eta;   Dret[4](1) = N1_xi  * dN1_eta;
        Dret[5](0) = dN0_xi  * N1_eta;   Dret[5](1) = N0_xi  * dN1_eta;
        Dret[6](0) = dNm1_xi * N1_eta;   Dret[6](1) = Nm1_xi * dN1_eta;
        Dret[7](0) = dNm1_xi * N0_eta;   Dret[7](1) = Nm1_xi * dN0_eta;
        Dret[8](0) = dN0_xi  * N0_eta;   Dret[8](1) = N0_xi  * dN0_eta;

        return Dret;
    }

public:
    virtual const std::array<Set::Vector,9> & getY0() {return Y0;}

private:
    //
    // Hard-coded values common to all elements
    //

    const std::array<Set::Vector,9> Y0 = {
        Set::Vector(-1.0,-1.0),
        Set::Vector( 0.0,-1.0),
        Set::Vector( 1.0,-1.0),
        Set::Vector( 1.0, 0.0),
        Set::Vector( 1.0, 1.0),
        Set::Vector( 0.0, 1.0),
        Set::Vector(-1.0, 1.0),
        Set::Vector(-1.0, 0.0),
        Set::Vector( 0.0, 0.0),
    };
};

}
#endif
