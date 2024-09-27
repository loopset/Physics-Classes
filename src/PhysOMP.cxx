#include "PhysOMP.h"

#include "PhysColors.h"

#include <cmath>
#include <iostream>

void PhysOMP::OMP::Print() const
{
    std::cout << BOLDYELLOW << "····· " << fName << " ·····" << '\n';
    std::cout << "   For target (A, Z) : (" << fA << ", " << fZ << ")" << '\n';
    std::cout << "   and projectile (A, Z) : (" << fTargetA << ", " << fTargetZ << ") @ " << fEnergy << " MeV" << '\n';
    std::cout << "-> Coulomb :" << '\n';
    std::cout << "   rc : " << frc << '\n';
    std::cout << "-> Real volume : " << '\n';
    std::cout << "   Vr : " << fVr << '\n';
    std::cout << "   rv : " << frv << '\n';
    std::cout << "   av : " << fav << '\n';
    std::cout << "-> Imaginary volume : " << '\n';
    std::cout << "   Wv : " << fWv << '\n';
    std::cout << "   rw : " << frw << '\n';
    std::cout << "   aw : " << faw << '\n';
    std::cout << "-> Imaginary surface : " << '\n';
    std::cout << "   Ws : " << fWs << '\n';
    std::cout << "   rs : " << frs << '\n';
    std::cout << "   as : " << fas << '\n';
    std::cout << "-> Real spin-orbit : " << '\n';
    std::cout << "   Vso : " << fVso << '\n';
    std::cout << "   FRESCO Vso / 2 : " << fVso / 2 << '\n';
    std::cout << "   rso : " << frvso << '\n';
    std::cout << "   aso : " << favso << '\n';
    std::cout << "-> Imaginary spin-orbit : " << '\n';
    std::cout << "   Wso : " << fWso << '\n';
    std::cout << "   FRESCO Wso / 2 : " << fWso / 2 << '\n';
    std::cout << "   rwso : " << frwso << '\n';
    std::cout << "   awso : " << fawso << RESET << '\n';
}

PhysOMP::Haixia::Haixia(int Z, int A, double energy) : OMP(Z, A, energy, 1, 2, "Haixia")
{
    // Coulomb
    frc = 1.303;
    // Real volume
    fVr = 91.85 - 0.249 * fEnergy + 0.000115 * std::pow(fEnergy, 2) + 0.642 * Z / std::pow(fA, 1. / 3);
    frv = 1.152 - 0.00776 * std::pow(fA, -1. / 3);
    fav = 0.719 + 0.0126 * std::pow(fA, 1. / 3);
    // Imaginary volume
    fWv = 1.104 + 0.0622 * fEnergy;
    frw = 1.305 + 0.0997 * std::pow(fA, -1. / 3);
    faw = 0.855 - 0.100 * std::pow(fA, 1. / 3);
    // Imaginary surface
    fWs = 10.83 - 0.0306 * fEnergy;
    frs = 1.334 + 0.152 * std::pow(fA, -1. / 3);
    fas = 0.531 + 0.062 * std::pow(fA, 1. / 3);
    // Real spin-orbit
    fVso = 3.557;
    frvso = 0.972;
    favso = 1.011;
}

PhysOMP::Pang::Pang(int Z, int A, double energy, bool isTriton) : OMP(Z, A, energy, 2, 3, "Pang"), fIsTriton(isTriton)
{
    // Change Z of target if is triton
    if(fIsTriton)
        fTargetZ = 1;
    // Set parameters
    frc = 1.24 + 0.12 * std::pow(fA, -1. / 3);
    auto Rc {frc * std::pow(fA, 1. / 3)};
    auto Ec {6 * fTargetZ * 1.44 / (5 * Rc)}; // Correction from Coulomb energy

    auto V0 {118.3};
    auto Ve {-0.13};
    fVr = V0 + Ve * (fEnergy - Ec);
    auto r0 {1.3};
    auto r00 {-0.48};
    frv = r0 + r00 * std::pow(fA, -1. / 3);
    fav = 0.82;

    auto Wv0 {38.5};
    auto Wve0 {156.1};
    auto Wvew {52.4};
    fWv = Wv0 / (1 + std::exp((Wve0 - (fEnergy - Ec)) / Wvew));

    auto Ws0 {35.0};
    auto Wst {34.2};
    auto epsilon {(double)(fN - fZ) / fA};
    auto Wse0 {30.8};
    auto Wsew {106.4};
    fWs = (Ws0 + (fIsTriton ? -1 : 1) * Wst * epsilon) / (1 + std::exp(((fEnergy - Ec) - Wse0) / Wsew));
    auto rww {1.31};
    auto rw0 {-0.13};
    frw = rww + rw0 * std::pow(fA, -1. / 3);
    frs = frw;
    faw = 0.84;
    fas = faw;
}
