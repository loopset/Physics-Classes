#include "PhysSF.h"

#include "PhysColors.h"

#include <iostream>

void PhysUtils::SpectroscopicFactor::Print() const
{
    std::cout << BOLDGREEN << "····· Spectroscopic factor ·····" << '\n';
    std::cout << "  SF      : " << fSF << " +- " << fUSF << '\n';
    std::cout << "  Chi2red : " << fChi2Red << '\n';
    std::cout << "  Ndf     : " << fNdf << RESET << '\n';
}
