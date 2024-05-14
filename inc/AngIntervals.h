#ifndef AngIntervals_h
#define AngIntervals_h
#include "ROOT/RDF/HistoModels.hxx"

#include "TCanvas.h"

#include <mutex>
#include <utility>
#include <vector>
namespace Angular
{
class Intervals
{
private:
    std::vector<std::pair<double, double>> fRanges {};
    std::vector<double> fOmegas {};
    std::vector<TH1D*> fHs {};
    std::mutex fMutex {}; //!< Mutex to ensure thread-safety of Fill function

public:
    Intervals(double xmin, double xmax, const ROOT::RDF::TH1DModel& model, double step = -1);

    // Setters
    void Fill(double thetaCM, double Ex);

    // Getters
    const std::vector<TH1D*> GetHistos() const { return fHs; }
    double GetLow(int i) const { return fRanges[i].first; }
    double GetCenter(int i) const { return (fRanges[i].first + fRanges[i].second) / 2; }
    double GetUp(int i) const { return fRanges[i].second; }
    double GetOmega(int i) const { return fOmegas[i]; }
    double GetMin() const { return fRanges.front().first; }
    double GetMax() const { return fRanges.back().second; }

    // Others
    TCanvas* Draw(const TString& title = "") const;

private:
    double ComputeSolidAngle(double min, double max);
};
} // namespace Angular

#endif // !AngInterval_h
