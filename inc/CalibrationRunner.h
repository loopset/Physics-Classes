#ifndef CalibrationRunner_h
#define CalibrationRunner_h

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TList.h"
#include "TSpectrum.h"

#include "CalibrationSource.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Calibration
{
class Runner
{
    typedef std::unordered_map<std::string, std::vector<std::shared_ptr<TF1>>> SatelliteCont;
    typedef std::unordered_map<std::string, std::shared_ptr<TF1>> GaussCont;

private:
    Source* fSource {}; // Pointer to previously defined Source obj
    TH1D* fData {};     // Pointer to hist containing channel data (could be rebinned)
    TH1D* fRawData {};  // Same but always unrebinned!
    bool fDebug {};     // Print info of calibration procedure

    // Settings for data
    std::pair<double, double> fRange {}; // Channel range to search for peaks
    // Settings for TSpectrum
    double fSpeSigma {2};
    double fSpeThresh {0.1};
    // Pedestal
    GaussCont fPedestals {};
    // Precalibration
    double fPreGaussWidth {15};
    GaussCont fGaussPre {};
    std::shared_ptr<TGraphErrors> fGraphPre {};
    std::shared_ptr<TF1> fCalibPre {};
    std::pair<double, double> fHistOpts {10., 0.015};
    std::shared_ptr<TH1D> fHistPre {};
    // Fine calibration
    GaussCont fGaussFine {};
    std::shared_ptr<TGraphErrors> fGraphFine {};
    std::shared_ptr<TF1> fCalibFine {};
    std::unordered_map<std::string, std::vector<std::string>> fSatelliteStr {};
    SatelliteCont fFineSat {};
    // Final: combination of both calibrations
    std::shared_ptr<TF1> fCalibFinal {};
    std::shared_ptr<TH1D> fHistFinal {};
    GaussCont fGaussFinal {};
    SatelliteCont fFinalSat {};
    // Fit settings
    // Graphs
    // Range is not important as we use the whole graph range
    // (had to delete R bc it caused wrong plotting in some case, altough fit was done perfectly)
    std::string fFitOptsGraph {"0QM+"}; // range determined by graph itself
    std::string fFitOptsGraphDebug {"0M+"};
    // WARNING: in some cases, the TGraphErrors fails when fitting (Etheo, Channel) in the precalibration
    // due to the errors in the X axis! This flag disables them ONLY for the precalibration
    bool fDisableXErrors {};
    // Histograms
    std::string fFitOpts {"0QRM+"};
    std::string fFitOptsDebug {"0RM+"};
    // Fitting ranges
    double fMinSigma {0.005};
    double fMaxSigma {2}; // MeV

public:
    Runner(Source* source, TH1D* data, TH1D* originalData = nullptr, bool debug = false)
        : fSource(source),
          fData(data),
          fRawData(originalData),
          fDebug(debug)
    {
        TH1::AddDirectory(false);
    }

    // Setters
    void SetRange(double min, double max)
    {
        fRange = {min, max};
        ApplyRange();
    }
    void SetGaussPreWidth(double w) { fPreGaussWidth = w; }
    void SertMinSigma(double sigma) { fMinSigma = sigma; }
    void SetMaxSigma(double sigma) { fMaxSigma = sigma; }

    void DisableXErrors() { fDisableXErrors = true; }

    // Main execution function
    bool DoIt();

    // Draw
    void Draw(TCanvas* c = nullptr) const { Debug(c); }

    // Methods to retrieve resolution
    void PrintRes() const;
    double GetRes(const std::string& peak) const
    {
        if(fGaussFinal.count(peak))
            return fGaussFinal.at(peak)->GetParameter(2);
        else
            return -1;
    }
    double GetURes(const std::string& peak) const
    {
        if(fGaussFinal.count(peak))
            return fGaussFinal.at(peak)->GetParError(2);
        else
            return -1;
    }

    // Get calibration parameters
    std::pair<double, double> GetParameters();
    std::pair<double, double> GetPedestal();
    std::shared_ptr<TH1D> GetHistFinal() const { return fHistFinal; }
    inline void DrawSatFinal() const { DrawSat(fFinalSat); }
    const SatelliteCont& GetFinalSat() const { return fFinalSat; }

private:
    void ApplyRange();
    void FindPedestal();
    bool DoPreCalibration();
    std::vector<std::pair<double, double>> FilterPeaks(const TSpectrum& spe);
    void FillHistPre();
    void DoFineCalibration();
    std::pair<double, double> GetAmpMeanInRange(TH1D* h, double min, double max);
    void AddSatellite(SatelliteCont& map, const GaussCont& funcs);
    void InitFinalCalib();
    void FillHistFinal();
    void DoFinalPlots();
    void Debug(TCanvas* c = nullptr) const;
    inline void DrawAll(TList* l, bool clone = false) const
    {
        for(auto* o : *l)
            if(o)
            {
                if(clone)
                    o->DrawClone("same");
                else
                    o->Draw("same");
            }
    }
    inline void DrawSat(const SatelliteCont& sats) const
    {
        for(const auto& [_, funcs] : sats)
            for(const auto& f : funcs)
            {
                f->DrawClone("same");
            }
    }
};
} // namespace Calibration

#endif // !CalibrationRunner_h
