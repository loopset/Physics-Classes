#ifndef AngIntervals_h
#define AngIntervals_h
#include "ROOT/RDF/HistoModels.hxx"

#include "TCanvas.h"
#include "TH1.h"

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
    std::vector<std::vector<TH1D*>> fHsPS {};
    std::mutex fMutex {}; //!< Mutex to ensure thread-safety of Fill function

public:
    Intervals(double xmin, double xmax, const ROOT::RDF::TH1DModel& model, double step = -1, int nps = 0);

    // Setters
    void Fill(double thetaCM, double Ex);
    void FillPS(int idx, double thetaCM, double Ex, double weight);
    // Templates are defined in class header
    // Otherwise we get a linking error!
    template <typename T>
    void FillConstPS(int ps, T* hf)
    {
        for(int i = 0; i < fRanges.size(); i++)
            fHsPS[ps][i]->Add(hf);
    }
    template <typename T>
    void FillIntervalPS(int ps, int iv, T* hf)
    {
        fHsPS[ps][iv]->Add(hf);
    }
    void TreatPS(int nsmooth = 10, double scale = 0.2);

    // Getters
    std::vector<TH1D*> GetHistos() const { return fHs; }
    std::vector<std::vector<TH1D*>> GetHistosPS() const { return fHsPS; }
    double GetLow(int i) const { return fRanges[i].first; }
    double GetCenter(int i) const { return (fRanges[i].first + fRanges[i].second) / 2; }
    double GetUp(int i) const { return fRanges[i].second; }
    double GetOmega(int i) const { return fOmegas[i]; }
    double GetMin() const { return fRanges.front().first; }
    double GetMax() const { return fRanges.back().second; }
    int GetSize() const { return fRanges.size(); }

    // Others
    TCanvas* Draw(const TString& title = "") const;

private:
    double ComputeSolidAngle(double min, double max);
};
} // namespace Angular

#endif // !AngInterval_h
