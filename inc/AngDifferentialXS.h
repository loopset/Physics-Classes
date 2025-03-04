#ifndef AngDifferentialXS_h
#define AngDifferentialXS_h

#include "TCanvas.h"
#include "TGraphErrors.h"

#include "AngFitter.h"
#include "AngIntervals.h"
#include "Interpolators.h"
#include "PhysExperiment.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace Angular
{
class DifferentialXS
{
private:
    // Pointers to objs defined outside
    Intervals* fIvs {};
    Fitter* fFitter {};
    Interpolators::Efficiency* fEff {};
    PhysUtils::Experiment* fExp {};

    // Results
    std::unordered_map<std::string, TGraphErrors*> fXS {};

    // Counts threshold
    double fNThresh {2};

public:
    DifferentialXS(Intervals* ivs, Interpolators::Efficiency* eff, PhysUtils::Experiment* exp)
        : fIvs(ivs),
          fEff(eff),
          fExp(exp)
    {
    }
    DifferentialXS(Intervals* ivs, Fitter* fits, Interpolators::Efficiency* eff, PhysUtils::Experiment* exp)
        : fIvs(ivs),
          fFitter(fits),
          fEff(eff),
          fExp(exp)
    {
    }

    // Main method
    void DoFor(const std::vector<std::string>& peaks = {});
    void DoFor(TGraphErrors* gexp, const std::string& peak);

    // Getters
    const std::unordered_map<std::string, TGraphErrors*>& GetAll() const { return fXS; }
    TGraphErrors* Get(const std::string& peak) const;
    double GetNThresh() const { return fNThresh; }

    // Setters
    void SetNThresh(double t) { fNThresh = t; }

    // Draw method
    TCanvas* Draw(const TString& title = "") const;

    // Write to file
    void Write(const std::string& dir, const std::string& name = "xs") const;

    // Function to trim some points in X, since in some cases the thetaCM range cannot be the same for all
    void TrimX(const std::string& peak, double xok, bool low = true);

private:
    void Do(const std::vector<double>& N, const std::string& peak);
    double
    Uncertainty(const std::string& peak, double N, double Nt, double Nb, double Omega, double eps, double thetaCM);
};
} // namespace Angular

#endif // !AngDifferentialXS_h
