#ifndef Interpolators_h
#define Interpolators_h

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

#include <string>
#include <unordered_map>
namespace Interpolators
{
class Efficiency
{
private:
    std::unordered_map<std::string, TEfficiency*> fEff {};
    std::unordered_map<std::string, TGraphAsymmErrors*> fGraph {};

public:
    Efficiency() = default;
    Efficiency(const std::string& peak, const std::string& file, const std::string& name = "eff")
    {
        Add(peak, file, name);
    };
    // Main method to add more peaks
    void Add(const std::string& peak, const std::string& file, const std::string& name = "eff");

    // Getters
    double GetPointEff(const std::string& peak, double thetaCM) const { return fGraph.at(peak)->Eval(thetaCM); }
    double GetMeanEff(const std::string& peak, double min, double max) const;
    TEfficiency* GetTEfficiency(const std::string& peak) const { return fEff.at(peak); }
    TGraphAsymmErrors* GetGraph(const std::string& peak) const { return fGraph.at(peak); }

    // Others
    TCanvas* Draw();

    void SaveAs(const std::string& file);
};
} // namespace Interpolators

#endif
