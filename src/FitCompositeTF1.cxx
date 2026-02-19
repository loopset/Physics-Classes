#include "FitCompositeTF1.h"

#include "TGraphErrors.h"
#include "TH1.h"
#include "TString.h"
#include "TVirtualFitter.h"

#include <initializer_list>
#include <string>
#include <vector>

Fitters::CompositeTF1::CompositeTF1(std::initializer_list<TF1*> ptrs)
{
    if(!ptrs.size())
        return;
    for(auto ptr : ptrs)
    {
        fKeys.push_back(ptr->GetName());
        fFuncs.push_back(ptr);
        fNPars += ptr->GetNpar();
        fRange = {ptr->GetXmin(), ptr->GetXmax()};
    }
    // Init composite model
    fComposite =
        new TF1 {"composite", [this](double* x, double* p) { return Eval(x, p); }, fRange.first, fRange.second, fNPars};
    // Init parameters
    Init();
}

Fitters::CompositeTF1::ParVec Fitters::CompositeTF1::UnpackPars(double* p)
{
    ParVec ret {};
    int idx {};
    for(auto func : fFuncs)
    {
        ret.push_back({});
        auto npar {func->GetNpar()};
        for(int i = 0; i < npar; i++)
        {
            ret.back().push_back(p[idx]);
            idx++;
        }
    }
    return ret;
}

double Fitters::CompositeTF1::Eval(double* x, double* p)
{
    double ret {};
    // Unpack parameters
    auto pars {UnpackPars(p)};
    // Call each function
    int idx {};
    for(auto func : fFuncs)
    {
        func->SetParameters(pars[idx].data());
        ret += func->Eval(x[0]);
        idx++;
    }
    return ret;
}

double Fitters::CompositeTF1::Eval(double x)
{
    auto pars {PackPars()};
    double xx[1] {x};
    return Eval(xx, pars.data());
}

std::vector<double> Fitters::CompositeTF1::PackPars()
{
    std::vector<double> pars {};
    for(auto func : fFuncs)
    {
        for(int i = 0; i < func->GetNpar(); i++)
            pars.push_back(func->GetParameter(i));
    }
    return pars;
}

void Fitters::CompositeTF1::Init()
{
    fComposite->SetParameters(PackPars().data());
    int idx {};
    for(auto func : fFuncs)
    {
        for(auto p = 0; p < func->GetNpar(); p++)
        {
            fComposite->SetParName(idx, func->GetParName(p));
            idx++;
        }
    }
}

void Fitters::CompositeTF1::InitDraw(TH1* hmodel)
{
    // Settings
    auto xmin {fComposite->GetXmin()};
    auto xmax {fComposite->GetXmax()};
    auto step {(xmax - xmin) / 300};
    // Global
    fGlobal = new TGraphErrors;
    for(double x = xmin; x <= xmax; x += step)
        fGlobal->SetPointX(fGlobal->GetN(), x);
    fGlobal->SetLineWidth(2);
    fGlobal->SetLineColor(kRed);
    // Add also 1 sigma band
    TVirtualFitter::GetFitter()->GetConfidenceIntervals(fGlobal, 0.68);
    // Functions
    if(!hmodel)
        return;
    for(auto func : fFuncs)
    {
        auto h {(TH1D*)hmodel->Clone(TString::Format("h%s", func->GetName()))};
        h->SetTitle(TString::Format("%s fit", func->GetName()));
        h->Reset();
        for(int b = 1; b <= h->GetNbinsX(); b++)
            h->SetBinContent(b, func->Eval(h->GetBinCenter(b)));
        h->SetLineWidth(2);
        fHs.push_back(h);
    }
}

void Fitters::CompositeTF1::Draw(const std::string& opts)
{
    if(!fGlobal)
        InitDraw();
    fGlobal->Draw("xl same");
    for(auto* h : fHs)
    {
        h->Draw("same plc pmc");
    }
}

std::vector<double> Fitters::CompositeTF1::Integral()
{
    std::vector<double> ret;
    for(auto func : fFuncs)
        ret.push_back(func->Integral(func->GetXmin(), func->GetXmax()));
    return ret;
}
