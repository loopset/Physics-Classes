#ifndef TwoFNR_cxx
#define TwoFNR_cxx

#include "TwoFNR.h"
#include "Rtypes.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TSpline.h"
#include "TString.h"
#include "TMultiGraph.h"
#include "TVirtualPad.h"
#include "TFitResult.h"
#include "TAxis.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

TGraphErrors* TheoreticalUtils::TwoFNR::ReadDiffXS(const std::string& file)
{
    auto g = new TGraphErrors(file.c_str(), "%lg %lg");
    return g;
}

TGraph* TheoreticalUtils::TwoFNR::ReadAsym(const std::string& file)
{
    std::ifstream streamer {file.c_str()};
    if(!streamer)
        throw std::runtime_error("File does not exist when reading TWOFNR Asymmetry!");
    double x {}; double y {}; double asym {};
    auto* ret {new TGraph()};
    while(streamer >> x >> y >> asym)
        ret->SetPoint(ret->GetN(), x, asym);
    return ret;
}

void TheoreticalUtils::TwoFNR::Add(const std::string &name, const std::string &file)
{
    //Graph
    fTheo[name] = ReadDiffXS(file);
    fTheo[name]->SetTitle(TString::Format("Diff. xs from %s;#theta_{CM} [degree];#frac{d#sigma}{d#Omega} [mb/sr]", name.c_str()));
    //Asymmetry graph
    fAsym[name] = ReadAsym(file);
    fAsym[name]->SetTitle(TString::Format("Asymmetry from %s;#theta_{CM} [degree];Asymmetry [mb/sr]", name.c_str()));
    //Save keys to keep ordering
    fKeys.push_back(name);
}

void TheoreticalUtils::TwoFNR::FitToExperimental(TGraphErrors *gexp, double xmin, double xmax)
{
    for(const auto& key : fKeys)
    {
        const auto& g {fTheo[key]};
        auto* sptheo {new TSpline3("sptheo", (TGraph*)g, "b2,e2", 0, 0)};
        auto* f1 {new TF1("fittheo", [&](double* x, double* p){return p[0] * sptheo->Eval(x[0]);}, 0, 180, 1)};
        f1->SetParameters(1);
        auto fitres = gexp->Fit(f1, "SQN", "", xmin, xmax);
        //Built output graph
        fFits[key] = new TGraphErrors();
        for(int p = 0; p < g->GetN(); p++)
            fFits[key]->SetPoint(p, g->GetPointX(p), g->GetPointY(p) * f1->GetParameter(0));
        fFits[key]->SetTitle(TString::Format("Fitted to exp. from %s", key.c_str()));
        //Store fit parameter
        fSF[key] = {f1->GetParameter(0), f1->GetParError(0)};
        //Print
        std::cout<<"============================"<<'\n';
        std::cout<<"-- Fit for "<<key<<'\n';
        std::cout<<"->C2S = "<<f1->GetParameter(0)<<" +/- "<<f1->GetParError(0)<<'\n';
        std::cout<<"->chisqr / dof = "<<fitres->Chi2() / fitres->Ndf()<<'\n';
        std::cout<<"============================"<<'\n';
        //Delete
        delete f1;
        delete sptheo;
    }
}

double TheoreticalUtils::TwoFNR::Integral(const std::string& key, double thetamin, double thetamax)
{
    double ret {};
    auto* g {new TGraph};
    for(int p = 0; p < fTheo.at(key)->GetN(); p++)
    {
        //Convert deg to rad in X axis of function
        auto x {fTheo[key]->GetPointX(p) * TMath::DegToRad()};
        auto y {fTheo[key]->GetPointY(p)};
        g->SetPoint(p, x, y);
    }
    // //To check that ROOT integration works
    // for(int p = 0; p < g->GetN(); p++)
    // {
    //     auto x {g->GetPointX(p)}; auto y {g->GetPointY(p)};
    //     if(x >= thetamin * TMath::DegToRad() && x < thetamax * TMath::DegToRad())
    //     {
    //         ret += TMath::TwoPi() * y * TMath::Sin(x) * 1 * TMath::DegToRad();
    //     }
    // }
    //Build TF1
    auto* spline {new TSpline3("spline", g, "b2,e2", 0, 0)};
    auto* func {new TF1("func", [&](double* x, double* p)
    {
        //One needs to take into account the solid angle element
        return TMath::TwoPi() * spline->Eval(x[0]) * TMath::Sin(x[0]);
    },
            0, TMath::TwoPi(), 1)};
    ret = func->Integral(thetamin * TMath::DegToRad(), thetamax * TMath::DegToRad());
    //deletes
    delete func;
    delete spline;
    delete g;
    //return and exit
    return ret;
}

void TheoreticalUtils::TwoFNR::IntegralAll(double thetamin, double thetamax, double exp, double uexp)
{
    std::cout<<"++++++++ TwoFNR abs xs computation "<<fName<<" ++++++++"<<'\n';
    for(const auto& key : fKeys)
    {
        auto res {Integral(key, thetamin, thetamax)};
        std::cout<<"Key = "<<key<<", integral = "<<res<<" mb"<<'\n';
        double ures {};
        if(uexp > 0)
            ures = 1. / res * uexp;
        if(exp > 0)
            std::cout<<"  Ratio to exp: C2S = "<<exp / res<<" +/- "<<ures<<'\n';
    }
    std::cout<<"+++++++++++++++++++++++++++++++++"<<'\n';
}

void TheoreticalUtils::TwoFNR::DrawTheoretical(const std::string& opts)
{
    int idx {24};//for marker style
    for(auto& [key, g] : fTheo)
    {
        g->SetLineWidth(2);
        g->SetMarkerStyle(0);
        if(opts.length() == 0)
            g->Draw("pl plc pmc same");
        else if(idx == 24)
            g->Draw(opts.c_str());
        else
            g->Draw((opts + " same").c_str());
        idx++;
    }
    //Fill colors
    gPad->Update();
    for(auto& [key, g] : fTheo)
        fStyle[key] = {g->GetLineColor(), g->GetMarkerStyle()};
}

void TheoreticalUtils::TwoFNR::DrawFitted()
{
    int idx {24};//for marker style
    for(auto& [key, g] : fFits)
    {
        g->SetLineWidth(2);
        g->SetLineStyle(7);
        if(fStyle.count(key))
        {
            g->SetLineColor(fStyle[key].first);
            g->SetMarkerColor(fStyle[key].first);
            g->SetMarkerStyle(fStyle[key].second);
            g->Draw("pl same");
        }
        else
        {
            g->SetMarkerStyle(0);
            g->Draw("pl plc pmc same");
        }
        idx++;
    }
}

TLegend* TheoreticalUtils::TwoFNR::DrawLegend(bool fancy)
{
    auto* leg {new TLegend(0.4, 0.2)};
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->SetNColumns(2);
    leg->SetColumnSeparation(0.25);
    //Add entries
    for(const auto& key : fKeys)
    {
        //Theoretical
        leg->AddEntry(fTheo[key], key.c_str(), "lp");
        //Fitted
        if(fSF.count(key))
        {
            TString name {};
            if(fancy)
                name = TString::Format("C^{2}S(%s) = %.3f \\pm %.3f", key.c_str(), fSF[key].first, fSF[key].second);
            else
                name = TString::Format("(%s) #times %.2f", key.c_str(), fSF[key].first);
            leg->AddEntry(fFits[key], name, "lp");
        }
        else
            leg->AddEntry(fFits[key], ("Fit. " + key).c_str(), "lp");
    }
    return leg;
}

TCanvas* TheoreticalUtils::TwoFNR::GetCanvas(TGraphErrors *gexp, double xmin, double xmax, bool asym)
{
    static int counterfnr {0};
    auto* cret {new TCanvas(TString::Format("cTheo%d", counterfnr), "TwoFNR canvas")};
    counterfnr++;
    cret->DivideSquare(1 + asym);
    //create multigraphs
    auto* mg {new TMultiGraph()};
    mg->SetTitle(TString::Format("%s;#theta_{CM} [#circ];d#sigma / d#Omega [mb/sr]", (fName.length() > 0) ? fName.c_str() : ""));
    auto* mas {new TMultiGraph()};
    mas->SetTitle(";#theta_{CM} [#circ];Asymmetry");
    //Legend
    auto* leg {new TLegend(0.3, 0.3)};
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    //Experimental
    if(gexp)
    {
        mg->Add(gexp);
        leg->AddEntry(gexp, "\\mathrm{Exp.}", "pe");
    }
    //Theoretical
    for(auto& [key, theo] : fTheo)
    {
        theo->SetLineWidth(2);
        theo->SetLineStyle(1);
        mg->Add(theo, "c");
        leg->AddEntry(theo, key.c_str(), "l");
    }
    //Asymmetry
    for(auto& [key, as] : fAsym)
    {
        as->SetLineWidth(2);
        as->SetLineStyle(7);
        mas->Add(as);
    }
    //Draw
    cret->cd(1);
    mg->Draw("apc plc pmc");
    leg->Draw();
    if(xmin != -1 && xmax != -1)
        mg->GetXaxis()->SetLimits(xmin, xmax);
    cret->cd(1)->Modified();
    cret->cd(1)->Update();
    
    if(asym)
    {
        cret->cd(2);
        mas->Draw("apl plc pmc");
    }
    return cret;
}

TCanvas* TheoreticalUtils::TwoFNR::GetCanvasPublication(TGraphErrors *gexp, double xmin, double xmax, const std::vector<int>& ls, const std::vector<int>& colors, bool sf)
{
    static int cpubcounter {0};
    auto* cret {new TCanvas(TString::Format("cPub%d", cpubcounter), "TwoFNR canvas for publication")};
    cpubcounter++;
    //create multigraphs
    auto* mg {new TMultiGraph()};
    mg->SetMinimum(1e-2);//to avoid any issues with log scale (bug of root)
    mg->SetTitle(TString::Format("%s;#theta_{CM} [#circ];d#sigma / d#Omega [mb/sr]", (fName.length() > 0) ? fName.c_str() : ""));
    //Legend
    auto* leg {new TLegend(0.3, 0.3)};
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    //Experimental
    if(gexp)
    {
        //Clone to avoid interferences with other functions
        auto* clone {(TGraphErrors*)gexp->Clone()};
        clone->SetLineColor(kBlack);
        mg->Add(clone, "pe");
        leg->AddEntry(clone, "\\mathrm{Exp.}", "pe");
    }
    //Style settings
    bool customColors {false};
    //Fitted
    int idx {1};
    for(const auto& key : fKeys)
    {
        //Clone to avoid interferences with Draw() functions
        auto* clone {(TGraphErrors*)fFits[key]->Clone()};
        clone->SetLineWidth(3);
        //Line style
        if(ls.size() > 0)
        {
            if(ls.size() != fKeys.size())
            {
                std::cout<<"Warning: Defaulting to automatic line style since size of specified ls does not match database!"<<'\n';
                clone->SetLineStyle(idx);
            }
            else
                clone->SetLineStyle(ls[idx - 1]);
        }
        else
            clone->SetLineStyle(idx);
        //Colors
        if(colors.size() == fKeys.size())
        {
            clone->SetLineColor(colors[idx - 1]);
            customColors = true;
        }
        mg->Add(clone, "c");//c line to smooth the angular intervals from the theoretical files
        //Append to legend
        TString entry;
        if(sf)
            entry = TString::Format("\\mathit{%s} \\Rightarrow %.2f", key.c_str(), fSF[key].first);
        else
            entry = key;
        leg->AddEntry(clone, entry, "l");
        idx++;
    }    
    //Draw
    cret->cd();
    mg->Draw((customColors) ? "apl" : "apc plc pmc");
    if(xmin != -1 && xmax != -1)
        mg->GetXaxis()->SetLimits(xmin, xmax);
    leg->Draw();
    cret->Modified();
    cret->Update();
    return cret;
}
#endif
