#include "PhysSF.h"

#include "TString.h"

#include "PhysColors.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

void PhysUtils::SpectroscopicFactor::Print() const
{
    std::cout << BOLDGREEN << "····· Spectroscopic factor ·····" << '\n';
    std::cout << "  SF      : " << fSF << " +- " << fUSF << '\n';
    std::cout << "  Chi2red : " << fChi2Red << '\n';
    std::cout << "  Ndf     : " << fNdf << RESET << '\n';
}

void PhysUtils::SFCollection::Add(const std::string& model, const SpectroscopicFactor& sf)
{
    fModels.push_back(model);
    fSFs.push_back(sf);
}

PhysUtils::SpectroscopicFactor* PhysUtils::SFCollection::Get(const std::string& model)
{
    auto it {std::find(fModels.begin(), fModels.end(), model)};
    if(it != fModels.end())
        return &(fSFs[std::distance(fModels.begin(), it)]);
    else
        throw std::invalid_argument("SFCollection::Get(): model not listed. You may try GetApprox if you want an "
                                    "approximative finding in vector");
}

PhysUtils::SpectroscopicFactor* PhysUtils::SFCollection::GetApprox(const std::string& model)
{
    auto it {std::find_if(fModels.begin(), fModels.end(),
                          [&](const std::string& str)
                          {
                              TString tstr {str};
                              return tstr.Contains(model);
                          })};
    if(it != fModels.end())
        return &fSFs[std::distance(fModels.begin(), it)];
    else
        throw std::invalid_argument("SFCollection::GetApprox(): cannot find any model with regex " + model);
}

PhysUtils::SpectroscopicFactor* PhysUtils::SFCollection::GetBestChi2()
{
    auto it {std::min_element(fSFs.begin(), fSFs.end(), [](const SpectroscopicFactor& a, const SpectroscopicFactor& b)
                              { return a.GetChi2Red() < b.GetChi2Red(); })};
    if(it != fSFs.end())
        return &(*it);
    else
        return nullptr; // not fitted so no minimum
}

void PhysUtils::SFCollection::Print() const
{
    std::cout << BOLDYELLOW << "----- SFCollection -----" << RESET << '\n';
    for(int i = 0; i < fModels.size(); i++)
    {
        std::cout << BOLDYELLOW << "-> Model : " << fModels[i] << RESET << '\n';
        fSFs[i].Print();
        std::cout << BOLDYELLOW << "--------------------" << RESET << '\n';
    }
}
