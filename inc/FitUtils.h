#ifndef FitUtils_h
#define FitUtils_h

#include "TGraph.h"
#include "TH1.h"

#include "FitModel.h"
#include "FitRunner.h"

namespace Fitters
{
void TreatPS(TH1D* hEx, TH1D* hPS);

void DrawGlobalFit(TGraph* g, const std::unordered_map<std::string, TH1D*>& hs);

void RunFit(TH1D* h, double exmin, double exmax, Fitters::Model& model, const Fitters::Runner::Init& initial,
            const Fitters::Runner::Bounds& bounds, const Fitters::Runner::Fixed& fixed, const std::string& outfile,
            const std::string& title = "");
} // namespace Fitters
#endif
