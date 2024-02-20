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

public:
    DifferentialXS(Intervals* ivs, Fitter* fits, Interpolators::Efficiency* eff, PhysUtils::Experiment* exp)
        : fIvs(ivs),
          fFitter(fits),
          fEff(eff),
          fExp(exp)
    {
    }

    // Main method
    void DoFor(const std::vector<std::string>& peaks = {});

    // Getters
    const std::unordered_map<std::string, TGraphErrors*>& GetAll() const { return fXS; }
    TGraphErrors* Get(const std::string& peak) const { return ((fXS.count(peak) ? fXS.at(peak) : nullptr)); }

    // Draw method
    TCanvas* Draw() const;

private:
    void Do(const std::string& peak);
    double Uncertainty(double N, double Nt, double Nb, double Omega, double eps, double thetaCM);
};
} // namespace Angular

#endif // !AngDifferentialXS_h
