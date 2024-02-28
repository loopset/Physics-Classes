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
    std::vector<Fitters::Data> fData {};
    std::vector<TH1D> fPS {};
    std::vector<Fitters::Model> fModels {};
    std::vector<TFitResult> fRes;

    // Declare integral vectors
    Counts fIgCounts {};
    Counts fSumCounts {};

    // Global fit read from file
    std::vector<std::string> fParNames {};
    TFitResult fGlobalFit {};

    // Settings to be sent to fitter
    bool fUseDivisions {true};
    bool fUseIntegral {};

public:
    Fitter(Intervals* ivs) : fIvs(ivs) {}

    // Setters
    void Configure(const std::string& file, const std::vector<TH1D>& ps = {});
    void SetUseDivisions(bool use) { fUseDivisions = use; }
    void SetUseIntegral(bool use) { fUseIntegral = use; }

    // Getters
    CountsIv GetIgCountsFor(const std::string& peak) const;
    CountsIv GetSumCountsFor(const std::string& peak) const;
    TGraphErrors* GetIgCountsGraph(const std::string& peak) const;
    TGraphErrors* GetSumCountsGraph(const std::string& peak) const;

    // Main methods
    void Run();
    void ComputeIntegrals(int nsigma = 2);

    TCanvas* Draw(const TString& title = "");
    TCanvas* DrawCounts(const TString& title = "");

private:
    void AddData(double exmin, double exmax);
    void AddModels();
    void ConfigRunner(Fitters::Runner& runner);
    void DoCounts(unsigned int iv, int nsigma);
    void CountsBySum(const std::string& key, unsigned int iv, int nsigma, TF1* f);
};
} // namespace Angular

#endif // !AngFitter_h
