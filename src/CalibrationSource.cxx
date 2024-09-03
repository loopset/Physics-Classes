#ifndef CalibrationSource_cxx
#define CalibrationSource_cxx

#include "CalibrationSource.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>

Calibration::Source::Source(const std::string& name) : fName(name)
{
    if(name == "3AlphaGanil")
        Add3AlphaGanil();
    else if(name == "3AlphaUSC")
        Add3AlphaUSC();
    else
        throw std::runtime_error("Calibration::Source(): source not recognized");
}

void Calibration::Source::Add3AlphaGanil()
{
    double commonSigma {0.1};
    // 239Pu
    fLabels.push_back("239Pu");
    fEnergies.push_back({});
    fSigmas.push_back({});
    fBR.push_back({});
    // 1--
    fEnergies.back().push_back(5.1055);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(11.94);
    // 2--
    fEnergies.back().push_back(5.1443);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(17.11);
    // 3--
    fEnergies.back().push_back(5.15659);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(70.77);

    // 241Am
    fLabels.push_back("241Am");
    fEnergies.push_back({});
    fSigmas.push_back({});
    fBR.push_back({});
    // 1--
    fEnergies.back().push_back(5.338);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(1.66);
    // 2--
    fEnergies.back().push_back(5.4428);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(14.76);
    // 3--Alpha
    fEnergies.back().push_back(5.48556);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(84.8);

    // 244Cm
    fLabels.push_back("244Cm");
    fEnergies.push_back({});
    fSigmas.push_back({});
    fBR.push_back({});
    //  1--
    fEnergies.back().push_back(5.76264);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(23.1);
    // 2--
    fEnergies.back().push_back(5.80477);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(76.9);

    // Set limits for fits!
    fLimits = {{"239Pu", {4.7, 5.3}}, {"241Am", {5.3, 5.65}}, {"244Cm", {5.65, 6.0}}};
}

void Calibration::Source::Add3AlphaUSC()
{
    // Same as GANIL one but changing 244Cm by 233U
    double commonSigma {0.1};
    // 233U
    fLabels.push_back("233U");
    fEnergies.push_back({});
    fSigmas.push_back({});
    fBR.push_back({});
    //  1--
    fEnergies.back().push_back(4.729);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(1.61);
    //  2--
    fEnergies.back().push_back(4.783);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(13.2);
    //  3--
    fEnergies.back().push_back(4.824);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(84.3);

    // 239Pu
    fLabels.push_back("239Pu");
    fEnergies.push_back({});
    fSigmas.push_back({});
    fBR.push_back({});
    // 1--
    fEnergies.back().push_back(5.1055);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(11.94);
    // 2--
    fEnergies.back().push_back(5.1443);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(17.11);
    // 3--
    fEnergies.back().push_back(5.15659);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(70.77);

    // 241Am
    fLabels.push_back("241Am");
    fEnergies.push_back({});
    fSigmas.push_back({});
    fBR.push_back({});
    // 1--
    fEnergies.back().push_back(5.338);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(1.66);
    // 2--
    fEnergies.back().push_back(5.4428);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(14.76);
    // 3--Alpha
    fEnergies.back().push_back(5.48556);
    fSigmas.back().push_back(commonSigma);
    fBR.back().push_back(84.8);

    // Set limits for fits!
    fLimits = {{"233U", {4.5, 4.95}}, {"239Pu", {4.95, 5.3}}, {"241Am", {5.3, 5.65}}};
}
std::unordered_map<std::string, double> Calibration::Source::GetMajorPeaks() const
{
    std::unordered_map<std::string, double> ret {};
    for(int iso = 0; iso < fLabels.size(); iso++)
    {
        auto name {fLabels[iso]};
        auto maxBR {std::distance(fBR[iso].begin(), std::max_element(fBR[iso].begin(), fBR[iso].end()))};
        ret[name] = fEnergies[iso][maxBR];
    }
    return ret;
}

void Calibration::Source::Print() const
{
    std::cout << "--- Calibration::Source: " << fName << " ----" << '\n';
    for(int s = 0; s < fLabels.size(); s++)
    {
        std::cout << "-> Source : " << fLabels[s] << '\n';
        for(int p = 0; p < fEnergies[s].size(); p++)
        {
            std::cout << "   Peak " << p << " -> E : " << fEnergies[s][p] << " MeV; sigma : " << fSigmas[s][p]
                      << " MeV; BR : " << fBR[s][p] << '\n';
        }
        std::cout << "   Limits : [" << fLimits.at(fLabels[s]).first << ", " << fLimits.at(fLabels[s]).second << "] MeV"
                  << '\n';
    }
}
#endif
