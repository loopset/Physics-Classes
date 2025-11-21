#ifndef AngFitter_h
#define AngFitter_h

#include "TCanvas.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TH1.h"

#include "AngIntervals.h"
#include "FitData.h"
#include "FitModel.h"
#include "FitRunner.h"

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
namespace Angular
{
class Fitter
{
public:
    using CountsIv = std::vector<double>;
    using Counts = std::unordered_map<std::string, CountsIv>;
    using WhichFree = std::vector<std::string>;

private:
    // Pointer to Intervals if called through this ctor
    Intervals* fIvs {};
    std::vector<TH1D*> fHistos {};
    std::vector<Fitters::Data> fData {};          //!
    std::vector<std::vector<TH1D*>> fHistosPS {}; //!
    std::vector<Fitters::Model> fModels {};       //!
    std::vector<TFitResult> fRes;
    // Saving the results of the fits
    std::vector<double> fResIvs {};
    std::vector<TGraph*> fResGlobal {};
    std::vector<std::vector<TH1D*>> fResHistos {};
    std::vector<std::vector<std::string>> fResNames {};

    // Declare integral vectors
    Counts fIgCounts {};  //!
    Counts fSumCounts {}; //!

    // Global fit read from file
    std::vector<std::string> fParNames {};
    TFitResult fGlobalFit {}; //!

    // Settings to be sent to fitter
    bool fUseDivisions {true};
    bool fUseIntegral {};
    // Allow variation of mean of gaussians during interval fit
    bool fAllowFreeMean {};
    double fFreeMeanRange {0.5};                // MeV
    std::vector<std::string> fWhichFreeMean {}; // For which states apply a free mean. If empty, for all!
    // Allow free sigma
    bool fAllowFreeSigma {};
    double fFreeSigmaRange {0.25}; // MeV
    WhichFree fWhichFreeSigma {};  // For which states apply free sigma. If empty, for all!
    // Allow free gamma
    bool fAllowFreeGamma {};
    double fFreeGammaRange {0.4}; // MeV
    WhichFree fWhichFreeGamma {};
    // Ignore PSs in fit by intervals
    // (it is interesing in some cases)
    bool fIgnorePS {false};
    // Or fix amplitude to ivs PS by given values
    std::map<int, std::vector<double>> fPSFixAmps {}; //!
    // Specify a manual range, do not use one used in config file
    std::pair<double, double> fManualRange {-1, -1};

public:
    Fitter() = default;
    Fitter(Intervals* ivs) : fIvs(ivs) {}

    // Setters
    void Configure(const std::string& file);
    void SetUseDivisions(bool use) { fUseDivisions = use; }
    void SetUseIntegral(bool use) { fUseIntegral = use; }
    // Free mean settings
    void SetAllowFreeMean(bool allow, const WhichFree& which = {})
    {
        fAllowFreeMean = allow;
        fWhichFreeMean = which;
    }
    void SetFreeMeanRange(double range) { fFreeMeanRange = range; }
    // Free sigma settings
    void SetAllowFreeSigma(bool allow, const WhichFree& which = {})
    {
        fAllowFreeSigma = allow;
        fWhichFreeSigma = which;
    }
    void SetFreeSigmaRange(double range) { fFreeSigmaRange = range; }
    // Gamma settings
    void SetAllowFreeGamma(bool allow, const WhichFree& which = {})
    {
        fAllowFreeGamma = allow;
        fWhichFreeGamma = which;
    }
    void SetGreeGammaRange(double range) { fFreeGammaRange = range; }
    // PS settings
    void SetIgnorePS(bool ignore) { fIgnorePS = ignore; }
    void SetFixAmpPS(int ips, const std::vector<double>& vals) { fPSFixAmps[ips] = vals; }
    // Specify a manual fitting range different than for the global fit
    void SetManualRange(double min, double max) { fManualRange = {min, max}; };

    // Getters
    CountsIv GetIgCountsFor(const std::string& peak) const;
    CountsIv GetSumCountsFor(const std::string& peak) const;
    TGraphErrors* GetIgCountsGraph(const std::string& peak) const;
    TGraphErrors* GetSumCountsGraph(const std::string& peak) const;
    const std::vector<std::string>& GetParNames() const { return fParNames; }
    const std::vector<TFitResult>& GetTFitResults() const { return fRes; }
    std::vector<double> GetResIvs() const { return fResIvs; }
    std::vector<TGraph*> GetResGlobal() const { return fResGlobal; }
    std::vector<std::vector<TH1D*>> GetResHistos() const { return fResHistos; }
    std::vector<std::vector<std::string>> GetResNames() const { return fResNames; }
    Intervals* GetIvs() const { return fIvs; }

    // Main methods
    void Run(bool print = false);
    void ComputeIntegrals(int nsigma = 2);

    TCanvas* Draw(const TString& title = "");
    TCanvas* DrawCounts(bool both = true, const TString& title = "");

    void Write(const std::string& file);
    void Read(const std::string& file);

    void Print() const;

private:
    void AddData(double exmin, double exmax);
    void AddModels();
    void ConfigRunner(int iv, Fitters::Runner& runner);
    void DoCounts(unsigned int iv, int nsigma);
    void CountsBySum(const std::string& key, unsigned int iv, int nsigma, TF1* f);
    void FillResHistos();
};
} // namespace Angular

#endif // !AngFitter_h
