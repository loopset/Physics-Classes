#ifndef TwoFNR_cxx
#define TwoFNR_cxx

#include "TwoFNR.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TString.h"
#include "TMultiGraph.h"
#include "TVirtualPad.h"
#include "TFitResult.h"
#include "TAxis.h"

#include <fstream>
#include <iostream>
#include <string>

TheoreticalUtils::TwoFNR::TwoFNR(const std::string& name, const std::string& file)
{
    Add(name, file);
}

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
    for(auto& [key, g] : fTheo)
    {
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
        fPar[key] = {f1->GetParameter(0), f1->GetParError(0)};
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

void TheoreticalUtils::TwoFNR::DrawTheoretical()
{
    int idx {24};//for marker style
    for(auto& [key, g] : fTheo)
    {
        g->SetLineWidth(2);
        g->SetMarkerStyle(idx);
        g->Draw("pl plc pmc same");
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
            g->SetMarkerStyle(idx);
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
    leg->SetTextSize(0.04);
    leg->SetNColumns(2);
    //Add entries
    for(const auto& key : fKeys)
    {
        //Theoretical
        leg->AddEntry(fTheo[key], key.c_str(), "lp");
        //Fitted
        if(fPar.count(key))
        {
            TString name {};
            if(fancy)
                name = TString::Format("C^{2}S(%s) = %.3f \\pm %.3f", key.c_str(), fPar[key].first, fPar[key].second);
            else
                name = TString::Format("(%s) \\times %.2f", key.c_str(), fPar[key].first);
            leg->AddEntry(fFits[key], name, "lp");
        }
        else
            leg->AddEntry(fFits[key], ("Fit. " + key).c_str(), "lp");
    }
    return leg;
}

TCanvas* TheoreticalUtils::TwoFNR::GetCanvas(TGraphErrors *gexp, bool asym)
{
    auto* cret {new TCanvas("cFNR", "TwoFNR canvas")};
    cret->DivideSquare(1 + asym);
    //create multigraphs
    auto* mg {new TMultiGraph()};
    mg->SetTitle(";#theta_{CM} [#circ];#frac{d#sigma}{d#Omega} [mb/sr]");
    auto* mas {new TMultiGraph()};
    mas->SetTitle(";#theta_{CM} [#circ];Asymmetry");
    //Experimental
    gexp->SetLineWidth(2);
    gexp->SetLineStyle(1);
    gexp->SetMarkerStyle(24);
    mg->Add(gexp);
    //Theoretical
    int idx {25};//for marker style
    for(auto& [key, theo] : fTheo)
    {
        theo->SetLineWidth(2);
        theo->SetMarkerStyle(idx);
        theo->SetLineStyle(2);
        mg->Add(theo);
        idx++;
    }
    //Fitted
    idx = 25;
    for(auto& [key, fit] : fFits)
    {
        fit->SetLineWidth(2);
        fit->SetMarkerStyle(idx);
        fit->SetLineStyle(4);
        mg->Add(fit);
        idx++;
    }
    //Asymmetry
    idx = 25;
    for(auto& [key, as] : fAsym)
    {
        as->SetLineWidth(2);
        as->SetMarkerStyle(idx);
        as->SetLineStyle(7);
        mas->Add(as);
        idx++;
    }
    //Draw
    cret->cd(1);
    mg->Draw("apl plc pmc");
    if(asym)
    {
        cret->cd(2);
        mas->Draw("apl plc pmc");
    }
    return cret;
}

TCanvas* TheoreticalUtils::TwoFNR::GetCanvasPublication(TGraphErrors *gexp, double xmin, double xmax)
{
    auto* cret {new TCanvas("cFNRPub", "TwoFNR canvas for publication")};
    //create multigraphs
    auto* mg {new TMultiGraph()};
    mg->SetTitle(";#theta_{CM} [#circ];#frac{d#sigma}{d#Omega} [mb/sr]");
    //Legend
    auto* leg {new TLegend(0.3, 0.3)};
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    //Experimental
    mg->Add(gexp);
    leg->AddEntry(gexp, "\\mathrm{Exp.}", "lpe");
    //Fitted
    for(auto& [key, fit] : fFits)
    {
        //Clone to avoid interferences with Draw() functions
        auto* clone {(TGraphErrors*)fit->Clone()};
        clone->SetLineWidth(2);
        clone->SetLineStyle(1);
        leg->AddEntry(clone, key.c_str(), "l");
        mg->Add(clone, "l");
    }
    //Draw
    cret->cd();
    mg->Draw("apl plc pmc");
    if(xmin != -1 && xmax != -1)
        mg->GetXaxis()->SetLimits(xmin, xmax);
    leg->Draw();
    cret->Modified();
    cret->Update();
    return cret;
}
#endif
