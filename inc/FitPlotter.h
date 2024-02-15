#ifndef FitPlotter_h
#define FitPlotter_h

#include "TF1.h"
#include "TFitResult.h"
#include "TGraph.h"

#include "FitData.h"
#include "FitModel.h"

#include <unordered_map>
namespace Fitters
{
class Plotter
{
private:
    // Just hold pointers to variables already defined in the macro
    Data* fData {};
    Model* fModel {};
    TFitResult* fRes {};

public:
    Plotter(Data* data, Model* model, TFitResult* res) : fData(data), fModel(model), fRes(res) {}

    // Getters
    TGraph* GetGlobalFit();

    std::unordered_map<std::string, TH1D*> GetIndividualHists();

private:
    void FillHistoFromFunc(TH1D* h, TF1* f);
};
} // namespace Fitters

#endif // !FitPlotter_h
