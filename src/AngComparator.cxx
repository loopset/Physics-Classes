#include "AngComparator.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TSpline.h"
#include "TString.h"
#include "TVirtualPad.h"

#include "PhysColors.h"

#include <iostream>
#include <memory>
#include <string>

void Angular::Comparator::Add(const std::string& name, const std::string& file)
{
    // So far, assume a twofnr file format
    fTheo[name] = ReadTwoFNR(file);
    // Set title
    fTheo[name]->SetTitle(name.c_str());
}

TGraphErrors* Angular::Comparator::ReadTwoFNR(const std::string& file)
{
    return new TGraphErrors {file.c_str(), "%lg %lg"};
}

TGraphErrors* Angular::Comparator::GetFitGraph(TGraphErrors* g, TF1* f)
{
    auto* ret {new TGraphErrors};
    for(int p = 0; p < g->GetN(); p++)
    {
        auto x {g->GetPointX(p)};
        auto y {g->GetPointY(p)};
        ret->SetPoint(ret->GetN(), x, y * f->GetParameter(0));
    }
    return ret;
}

void Angular::Comparator::Fit(double xmin, double xmax)
{
    // Set fit range for later
    fFitRange = {xmin, xmax};
    // Fit is based on a TSpline
    for(const auto& [name, gt] : fTheo)
    {
        auto spline {std::make_unique<TSpline3>("spline", (TGraph*)gt, "b2,e2", 0, 0)};
        auto func {std::make_unique<TF1>(
            "func", [&](double* x, double* p) { return p[0] * spline->Eval(x[0]); }, 0, 180, 1)};
        // Set parameters
        func->SetParameters(1);
        func->SetParName(0, "SF");
        // And fit exp to theo!
        auto res {fExp->Fit(func.get(), "SQN", "", xmin, xmax)};
        // Add fit graph
        fFit[name] = GetFitGraph(gt, func.get());
        fFit[name]->SetTitle(gt->GetTitle());
        // And store results in map
        fRes[name] = res;
    }
    Print();
}

void Angular::Comparator::Print() const
{
    std::cout << BOLDYELLOW << "···· Comparator for " << fName << " ····" << '\n';
    for(const auto& [name, res] : fRes)
    {
        std::cout << "-> Model : " << name << '\n';
        std::cout << "   SF    : " << res->Value(0) << " +/- " << res->Error(0) << '\n';
        auto chi2 {res->Chi2()};
        auto ndf {res->Ndf()};
        std::cout << "   chi2  : " << chi2 << '\n';
        std::cout << "   ndf   : " << ndf << '\n';
        std::cout << "   chi2 / ndf : " << chi2 / ndf << '\n';
    }
    std::cout << "······························" << RESET << '\n';
}

TLegend* Angular::Comparator::BuildLegend(double width, double height)
{
    auto* l {new TLegend {width, height}};
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(42);
    l->SetTextSize(0.04);
    return l;
}

TCanvas* Angular::Comparator::Draw(bool withSF)
{
    // Draw all using a TMultiGraph
    auto* mg {new TMultiGraph};
    mg->SetTitle((fName + " ;#theta_{CM} [#circ];d#sigma / d#Omega [mb / sr]").c_str());
    // Create a legend
    auto* leg {BuildLegend()};
    // 1-> Add experimental
    fExp->SetMarkerStyle(25);
    fExp->SetLineWidth(2);
    leg->AddEntry(fExp, "Exp", "pe");
    mg->Add(fExp, "p");
    // 2-> Add all fitted
    for(const auto& [name, g] : fFit)
    {
        g->SetLineWidth(2);
        TString desc {name};
        if(withSF)
            desc += TString::Format(" #Rightarrow SF = %.2f", fRes[name]->Value(0));
        leg->AddEntry(g, desc, "l");
        mg->Add(g, "c");
    }
    // Plot
    auto* c {new TCanvas {"cComp", "Theo to exp comp"}};
    mg->Draw("a plc pmc");
    mg->GetXaxis()->SetLimits(fFitRange.first, fFitRange.second);
    c->cd()->Update();
    leg->Draw();
    // Somehow GetXaxis sets the selected pad and this causes
    // all posterior DrawClones to be drawn on this canvas
    // reset it!
    gROOT->SetSelectedPad(nullptr);
    return c;
}

TCanvas* Angular::Comparator::DrawTheo()
{
    // Create multigraph
    auto* mg {new TMultiGraph};
    mg->SetTitle((fName + " models;#theta_{CM} [#circ];d#sigma / d#Omega [mb / sr]").c_str());
    // Legend
    auto* leg {BuildLegend()};
    // Add graphs
    for(const auto& [name, g] : fTheo)
    {
        g->SetLineWidth(2);
        leg->AddEntry(g, name.c_str(), "l");
        mg->Add(g, "c");
    }
    // Plot
    auto* c {new TCanvas {"cCompTheo", "Theo models"}};
    mg->Draw("a plc");
    leg->Draw();
    return c;
}
