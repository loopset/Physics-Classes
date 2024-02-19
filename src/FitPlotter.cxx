#include "FitPlotter.h"

#include "Rtypes.h"

#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"

#include <string>
#include <unordered_map>

TGraph* Fitters::Plotter::GetGlobalFit()
{
    // Init return object
    auto* ret {new TGraph};
    // Set use spline to not have gaps in plotting
    fModel->SetUseSpline(true);
    for(auto x = fData->GetXLow(); x < fData->GetXUp(); x += fData->GetBinWidth() / 10)
    {
        ret->SetPoint(ret->GetN(), x, (*fModel)(&x, fRes->GetParams()));
    }
    // Disable use spline
    fModel->SetUseSpline(false);
    // A few default settings
    ret->SetLineWidth(2);
    ret->SetLineColor(kRed);
    return ret;
}

void Fitters::Plotter::FillHistoFromFunc(TH1D* h, TF1* f)
{
    for(int bin = 1; bin <= h->GetNbinsX(); bin++)
    {
        auto x {h->GetXaxis()->GetBinCenter(bin)};
        auto y {f->Eval(x)};
        h->SetBinContent(bin, y);
    }
}

std::unordered_map<std::string, TH1D*> Fitters::Plotter::GetIndividualHists()
{
    std::unordered_map<std::string, TH1D*> ret;
    // Unpack parameters
    auto pack = fModel->UnpackParameters(fRes->GetParams());
    auto gaus = pack[0];
    auto voigt = pack[1];
    auto phase = pack[2];
    auto cte = pack[3];
    // Gauss
    int idx = 0;
    for(const auto& pars : gaus)
    {
        std::string key {"g" + std::to_string(idx)};
        // Function
        auto f = new TF1(key.c_str(), "gaus", fData->GetXLow(), fData->GetXUp());
        f->SetParameters(&pars[0]);
        // Histogram
        ret[key] = new TH1D(("h" + key).c_str(), key.c_str(), fData->GetSize(), fData->GetXLow(), fData->GetXUp());
        FillHistoFromFunc(ret[key], f);
        delete f;
        idx++;
    }
    // Voigt
    idx = 0;
    for(const auto& pars : voigt)
    {
        std::string key {"v" + std::to_string(idx)};
        // Function
        auto f = new TF1(key.c_str(), "[0] * TMath::Voigt(x - [1], [2], [3])", fData->GetXLow(), fData->GetXUp());
        f->SetParameters(&pars[0]);
        // Histogram
        ret[key] = new TH1D(("h" + key).c_str(), key.c_str(), fData->GetSize(), fData->GetXLow(), fData->GetXUp());
        FillHistoFromFunc(ret[key], f);
        delete f;
        idx++;
    }
    // Phase space
    idx = 0;
    for(const auto& pars : phase)
    {
        std::string key {"ps" + std::to_string(idx)};
        ret[key] = new TH1D(("h" + key).c_str(), key.c_str(), fData->GetSize(), fData->GetXLow(), fData->GetXUp());
        // Manually fill it
        for(int bin = 1; bin <= ret[key]->GetNbinsX(); bin++)
        {
            auto x {ret[key]->GetXaxis()->GetBinCenter(bin)};
            auto y {pars.front() * fModel->EvalPS(idx, x)};
            ret[key]->SetBinContent(bin, y);
        }
        idx++;
    }
    // Cte
    idx = 0;
    for(const auto& pars : cte)
    {
        std::string key {"cte" + std::to_string(idx)};
        // Function
        auto f = new TF1(key.c_str(), "[0]", fData->GetXLow(), fData->GetXUp());
        f->SetParameters(&pars[0]);
        // Histogram
        ret[key] = new TH1D(("h" + key).c_str(), key.c_str(), fData->GetSize(), fData->GetXLow(), fData->GetXUp());
        FillHistoFromFunc(ret[key], f);
        delete f;
        idx++;
    }
    // Set directory to null
    for(auto& [_, h] : ret)
        h->SetDirectory(nullptr);
    return std::move(ret);
}
