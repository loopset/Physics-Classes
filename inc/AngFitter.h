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
    typedef std::vector<double> CountsIv;
    typedef std::unordered_map<std::string, CountsIv> Counts;

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
    std::vector<std::string> fParNames {}; //!
    TFitResult fGlobalFit {};              //!

    // Settings to be sent to fitter
    bool fUseDivisions {true};
    bool fUseIntegral {};
    // Allow variation of mean of gaussians during interval fit
    bool fAllowFreeMean {};
    double fFreeMeanRange {0.5}; // MeV
    bool fAllowFreeSigma {};
    double fFreeSigmaRange {0.25}; // MeV
    // Ignore PSs in fit by intervals
    // (it is interesing in some cases)
    bool fIgnorePS {false};
    // Specify a manual range, do not use one used in config file
    std::pair<double, double> fManualRange {-1, -1};

public:
    Fitter() = default;
    Fitter(Intervals* ivs) : fIvs(ivs) {}

    // Setters
    void Configure(const std::string& file);
    void SetUseDivisions(bool use) { fUseDivisions = use; }
    void SetUseIntegral(bool use) { fUseIntegral = use; }
    void SetAllowFreeMean(bool allow) { fAllowFreeMean = allow; }
    void SetFreeMeanRange(double range) { fFreeMeanRange = range; }
    void SetIgnorePS(bool ignore) { fIgnorePS = ignore; }
    void SetManualRange(double min, double max) { fManualRange = {min, max}; };
    void SetAllowFreeSigma(bool allow) { fAllowFreeSigma = allow; }
    void SetFreeSigmaRange(double range) { fFreeSigmaRange = range; }

    // Getters
    CountsIv GetIgCountsFor(const std::string& peak) const;
    CountsIv GetSumCountsFor(const std::string& peak) const;
    TGraphErrors* GetIgCountsGraph(const std::string& peak) const;
    TGraphErrors* GetSumCountsGraph(const std::string& peak) const;
    std::vector<double> GetResIvs() const { return fResIvs; }
    std::vector<TGraph*> GetResGlobal() const { return fResGlobal; }
    std::vector<std::vector<TH1D*>> GetResHistos() const { return fResHistos; }
    std::vector<std::vector<std::string>> GetResNames() const { return fResNames; }
    Intervals* GetIvs() const { return fIvs; }

    // Main methods
    void Run();
    void ComputeIntegrals(int nsigma = 2);

    TCanvas* Draw(const TString& title = "");
    TCanvas* DrawCounts(bool both = true, const TString& title = "");

    void Write(const std::string& file);
    void Read(const std::string& file);

private:
    void AddData(double exmin, double exmax);
    void AddModels();
    void ConfigRunner(Fitters::Runner& runner);
    void DoCounts(unsigned int iv, int nsigma);
    void CountsBySum(const std::string& key, unsigned int iv, int nsigma, TF1* f);
    void FillResHistos();
};
} // namespace Angular

#endif // !AngFitter_h
