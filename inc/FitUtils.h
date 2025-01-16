#ifndef FitUtils_h
#define FitUtils_h

#include "TGraph.h"
#include "TH1.h"
#include "TLegend.h"

#include "FitModel.h"
#include "FitRunner.h"

#include <string>
#include <unordered_map>

namespace Fitters
{
void TreatPS(TH1D* hEx, TH1D* hPS, int nsmooth = 10, double scale = 0.2);

void DrawGlobalFit(TGraph* g, const std::unordered_map<std::string, TH1D*>& hs, TLegend* leg,
                   const std::unordered_map<std::string, std::string>& labels = {});

void SaveGlobalFit(const std::string& file, TH1D* h, TGraph* g, const std::unordered_map<std::string, TH1D*>& hs,
                   TLegend* leg);

void RunFit(TH1D* h, double exmin, double exmax, Fitters::Model& model, const Fitters::Runner::Init& initial,
            const Fitters::Runner::Bounds& bounds, const Fitters::Runner::Fixed& fixed, const std::string& outfile,
            const std::string& title = "", const std::unordered_map<std::string, std::string>& labels = {},
            const Fitters::Runner::Step& steps = {});

std::pair<Fitters::Runner::Init, Fitters::Runner::Init> ReadInit(const std::string& name);

void ReadDrawGlobalFit(const std::string& file);
} // namespace Fitters
#endif
