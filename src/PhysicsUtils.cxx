#ifndef PhysicsUtils_cxx
#define PhysicsUtils_cxx

#include "PhysicsUtils.h"
#include "TMath.h"
#include <iostream>
#include <sstream>
#include <string>

double PhysicsUtils::ExperimentInfo::GetUNb() const
{
    if(fDiv == 0)
        std::cout<<"Warning: CATS division factor = 0! Set it!"<<'\n';
    double uct {TMath::Sqrt(fNtriggers) * fDiv *fDiv};
    return uct;
}

std::string PhysicsUtils::Uncertainty::GetStr() const
{
    return {};
}
#endif
