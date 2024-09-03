#ifndef Interpolators_h
#define Interpolators_h

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include <string>
#include <unordered_map>
namespace Interpolators
{
class Efficiency
{
private:
    std::unordered_map<std::string, TEfficiency*> fEff {};
    std::unordered_map<std::string, TGraphAsymmErrors*> fGraph {};
    double fMeanStep {0.5}; // step in deg to compute mean efficiency

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
    double GetPointUEff(const std::string& peak, double thetaCM) const;

    // Others
    TCanvas* Draw(bool multigraph = true, const TString& title = "");

    void SaveAs(const std::string& file);
};

class Sigmas
{
private:
    TGraphErrors* fGraph {};

public:
    Sigmas() = default;
    Sigmas(TGraphErrors* gs) : fGraph(gs) {}
    Sigmas(const std::string& file, const std::string& name = "gsigmas") { Read(file, name); }

    // Write
    void SaveAs(const std::string& file, const std::string& name = "gsigmas");
    // Read
    void Read(const std::string& file, const std::string& name = "gsigmas");

    // Main functions
    double Eval(double Ex) { return fGraph->Eval(Ex); }

    void Draw();

    TGraphErrors* GetGraph() const { return fGraph; }
};
} // namespace Interpolators

#endif
