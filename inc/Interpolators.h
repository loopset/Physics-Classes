#ifndef Interpolators_h
#define Interpolators_h

#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

#include <string>
namespace Interpolators
{
class Efficiency
{
private:
    TEfficiency* fEff {};
    TGraphAsymmErrors* fGraph {};

public:
    Efficiency(const std::string& file, const std::string& key = "eff");

    // Getters
    double GetPointEff(double thetaCM) const { return fGraph->Eval(thetaCM); }
    double GetMeanEff(double min, double max) const;
    TEfficiency* GetTEfficiency() const { return fEff; }
    TGraphAsymmErrors* GetGraph() const { return fGraph; }
};
} // namespace Interpolators

#endif
