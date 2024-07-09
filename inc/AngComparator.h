#ifndef AngComparator_h
#define AngComparator_h

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFitResultPtr.h"
#include "TGraphErrors.h"
#include "TLegend.h"

#include "PhysExperiment.h"

#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
namespace Angular
{
class Comparator
{
private:
    std::string fName {};
    TGraphErrors* fExp {};
    // Vector to hold keys (so they are plotted in the order given by Add())
    std::vector<std::string> fKeys {};
    std::unordered_map<std::string, TGraphErrors*> fTheo {};
    // Graph fit results from fit
    std::unordered_map<std::string, TGraphErrors*> fFit {};
    // Fit results, containing SF
    std::unordered_map<std::string, TFitResultPtr> fRes {};
    // Store also fitting range
    std::pair<double, double> fFitRange {-1, -1};
    // To format lines...
    std::unordered_map<std::string, std::tuple<int, int, int>> fStyles {};

public:
    Comparator(const std::string& name, TGraphErrors* exp) : fName(name), fExp((TGraphErrors*)exp->Clone()) {}

    // Main method to add theoretical models
    void Add(const std::string& name, const std::string& file, int lc = -1, int ls = 1, int lw = 3);

    // Fit all models to exp in given range
    void Fit(double xmin = -1, double xmax = -1);

    // Print SFs
    void Print() const;

    // Main draw method
    TCanvas* Draw(const TString& title = "", bool logy = false, bool withSF = true, double offset = 3);

    // Draw theoretical and fits
    TCanvas* DrawTheo();

    // Compare to experimental counts to check efficiency
    TCanvas* ScaleToExp(const std::string& model, PhysUtils::Experiment* exp, TGraphErrors* gcounts,
                        TEfficiency* teff = nullptr, double SF = -1);

    // Canvas with ratio per point
    TCanvas* QuotientPerPoint();

    // Retrieve SF
    double GetSF(const std::string& model);

private:
    TGraphErrors* ReadTwoFNR(const std::string& file);
    TGraphErrors* GetFitGraph(TGraphErrors* g, TF1* f);
    TLegend* BuildLegend(double width = 0.25, double height = 0.2);
};
} // namespace Angular

#endif // !AngComparator_h
