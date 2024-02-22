#ifndef AngComparator_h
#define AngComparator_h

#include "TCanvas.h"
#include "TFitResultPtr.h"
#include "TGraphErrors.h"
#include "TLegend.h"

#include <string>
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
    std::pair<double, double> fFitRange {};

public:
    Comparator(const std::string& name, TGraphErrors* exp) : fName(name), fExp((TGraphErrors*)exp->Clone()) {}

    // Main method to add theoretical models
    void Add(const std::string& name, const std::string& file);

    // Fit all models to exp in given range
    void Fit(double xmin = -1, double xmax = -1);

    // Print SFs
    void Print() const;

    // Main draw method
    TCanvas* Draw(bool withSF = true, double offset = 3);

    // Draw theoretical and fits
    TCanvas* DrawTheo();

private:
    TGraphErrors* ReadTwoFNR(const std::string& file);
    TGraphErrors* GetFitGraph(TGraphErrors* g, TF1* f);
    TLegend* BuildLegend(double width = 0.25, double height = 0.2);
};
} // namespace Angular

#endif // !AngComparator_h
