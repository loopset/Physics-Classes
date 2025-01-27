#ifndef Interpolators_h
#define Interpolators_h

#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TVirtualPad.h"

#include <string>
#include <unordered_map>
namespace Interpolators
{
class Efficiency
{
private:
    std::unordered_map<std::string, TEfficiency*> fEff {};
    std::unordered_map<std::string, TGraphAsymmErrors*> fGraph {};
    int fMeanDiv {5};
    bool fIsLab {}; //!< Infer from efficiency name if we're in Lab frame instead of standard CM

public:
    Efficiency() = default;
    Efficiency(const std::string& peak, const std::string& file, const std::string& name = "eff")
    {
        Add(peak, file, name);
    };
    // Main method to add more peaks
    void Add(const std::string& peak, const std::string& file, const std::string& name = "eff");
    void Add(const std::string& peak, TEfficiency* eff);

    // Getters
    double GetPointEff(const std::string& peak, double thetaCM) const { return fGraph.at(peak)->Eval(thetaCM); }
    double GetMeanEff(const std::string& peak, double min, double max) const;
    TEfficiency* GetTEfficiency(const std::string& peak) const { return fEff.at(peak); }
    TGraphAsymmErrors* GetGraph(const std::string& peak) const { return fGraph.at(peak); }
    double GetPointUEff(const std::string& peak, double thetaCM) const;

    // Setters
    void SetMeanDiv(int div) { fMeanDiv = div; }

    // Others
    TVirtualPad* Draw(bool multigraph = true, const TString& title = "", TVirtualPad* par = nullptr);

    void SaveAs(const std::string& file);

private:
    void AddGraph(const std::string& peak, TEfficiency* eff);
    std::pair<double, double> GetXaxisRange(TMultiGraph* mg);
    double GetYaxisRange(TMultiGraph* mg);
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
