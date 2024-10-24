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

PhysOMP::Daehnick::Daehnick(int Z, int A, double energy) : OMP(Z, A, energy, 1, 2, "Daehnick")
{
    // Using option L in paper: NON-RELATIVISTIC KINEMATICS
    double beta {-std::pow(fEnergy / 100, 2)};
    auto mu {[](double N)
             {
                 double ret;
                 for(const auto& m : {8, 20, 28, 50, 82, 126})
                     ret += std::exp(-1 * std::pow((m - N) / 2, 2));
                 return ret;
             }};
    // Coulomb
    frc = 1.3;
    // Real volume
    fVr = 88.5 - 0.26 * fEnergy + 0.88 * fZ * std::pow(fA, -1. / 3);
    frv = 1.17;
    fav = 0.709 + 0.0017 * fEnergy;
    // Imaginary volume
    fWv = (12.2 + 0.026 * fEnergy) * (1 - std::exp(beta));
    frw = 1.325;
    faw = 0.53 + 0.07 * std::pow(fA, 1. / 3) - 0.04 * mu(fN);
    // Imaginary surface
    fWs = (12.2 + 0.026 * fEnergy) * std::exp(beta);
    frs = frw;
    fas = faw;
    // Real spin-orbit
    fVso = 7.33 - 0.029 * fEnergy;
    frvso = 1.07;
    favso = 0.66;
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


PhysOMP::HT1p::HT1p(int Z, int A, double energy, bool isTriton) : OMP(Z, A, energy, 2, 3, "HT1p"), fIsTriton(isTriton)
{
    // Change Z of target if is triton
    if(fIsTriton)
        fTargetZ = 1;
    // Coulomb
    frc = 1.3;
    // Correction to energy
    double Rc {frc * std::pow(fA, 1. / 3)};
    double Ec {6 * (fIsTriton ? 1 : 2) * fZ * 1.44 / (5 * Rc)};

    // Real volume
    double V0 {155.1};
    double Ve {-0.678};
    fVr = V0 + Ve * (fEnergy - Ec);
    double r0 {0.920};
    double r00 {0.108};
    double r0e {0.0031};
    frv = r0 + (r00 + r0e * (fEnergy - Ec)) * std::pow(fA, -1. / 3);
    fav = 0.792;

    // Imaginary volume
    double Wv0 {33.1};
    double Wve0 {156.1};
    double Wvew {52.4};
    fWv = Wv0 / (1 + std::exp((Wve0 - (fEnergy - Ec)) / Wvew));
    double rww {1.43};
    double rw0 {-0.16};
    frw = rww + rw0 * std::pow(fA, -1. / 3);
    faw = 0.801;

    // Imaginary surface
    double Ws0 {21.8};
    double Wst {13.1};
    double epsilon {(double)(fN - fZ) / fA};
    double Wse0 {30.8};
    double Wsew {106.4};
    fWs = (Ws0 + (fIsTriton ? -1 : 1) * Wst * epsilon) / (1 + std::exp(((fEnergy - Ec) - Wse0) / Wsew));
    frs = frw;
    fas = faw;
}

PhysOMP::KoningDelaroche::KoningDelaroche(int Z, int A, double energy) : OMP(Z, A, energy, 1, 1, "KonigDelaroche")
{
    // Some parameters
    double Efermi {-8.4075 + 0.01378 * fA};
    // Substract fermi energy to proton's
    fEnergy = fEnergy - Efermi;

    // Coulomb
    frc = 1.198 + 0.697 * std::pow(fA, -2. / 3) + 12.994 * std::pow(fA, -5. / 3);

    // Declare parameters
    double v1p {59.30 + 21 * (double)(fN - fZ) / fA - 0.024 * fA};
    double v2p {0.007067 + 4.23e-6 * fA};
    double v3p {1.729e-5 + 1.136e-8 * fA};
    double v4p {7e-9};
    double vbarc {1.73 / frc * fZ * std::pow(fA, -1. / 3)}; // Coulomb correction

    double w1p {14.667 + 0.009629 * fA};
    double w2p {73.55 + 0.0795 * fA};

    double d1p {16 - 16 * (double)(fN - fZ) / fA};
    double d2p {0.018 + 0.003802 / (1 + std::exp((double)(fA - 156.) / 8))};
    double d3p {11.5};

    double vso1 {5.922 + 0.0030 * fA};
    double vso2 {0.0040};

    double wso1 {-3.1};
    double wso2 {160};

    // 1-> Volume real
    fVr = v1p * (1 - v2p * fEnergy + v3p * std::pow(fEnergy, 2) - v4p * std::pow(fEnergy, 3)) +
          vbarc * v1p * (v2p - 2 * v3p * fEnergy + 3 * v4p * std::pow(fEnergy, 2));
    frv = 1.3039 - 0.4054 * std::pow(fA, -1. / 3);
    fav = 0.6778 - 1.487e-4 * fA;

    // 2-> Volume imaginary
    fWv = w1p * std::pow(fEnergy, 2) / (std::pow(fEnergy, 2) + std::pow(w2p, 2));
    frw = frv;
    faw = fav;

    // 3-> Surface imaginary
    fWs = d1p * std::pow(fEnergy, 2) * std::exp(-d2p * fEnergy) / (std::pow(fEnergy, 2) + std::pow(d3p, 2));
    frs = 1.3424 - 0.01585 * std::pow(fA, 1. / 3);
    fas = 0.5187 + 5.205e-4 * fA;

    // 4-> SO real
    fVso = vso1 * std::exp(-vso2 * fEnergy);
    frvso = 1.1854 - 0.647 * std::pow(fA, -1. / 3);
    favso = 0.59;

    // 5-> SO imaginary
    fWso = wso1 * std::pow(fEnergy, 2) / (std::pow(fEnergy, 2) + std::pow(wso2, 2));
    frwso = frvso;
    fawso = favso;
}
