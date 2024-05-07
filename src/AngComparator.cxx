#include "AngComparator.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TSpline.h"
#include "TString.h"
#include "TVirtualPad.h"

#include "PhysColors.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

void Angular::Comparator::Add(const std::string& name, const std::string& file)
{
    // So far, assume a twofnr file format
    fTheo[name] = ReadTwoFNR(file);
    // Set title
    fTheo[name]->SetTitle(name.c_str());
    // And add key
    fKeys.push_back(name);
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
    // Set auto range if not passed
    if(xmin == -1 || xmax == -1)
    {
        xmin = fExp->GetPointX(0);
        xmax = fExp->GetPointX(fExp->GetN() - 1);
    }
    // Set fit range for later
    fFitRange = {xmin, xmax};
    // Fit is based on a TSpline
    for(const auto& name : fKeys)
    {
        auto& gt {fTheo[name]};
        auto spline {std::make_unique<TSpline3>("spline", (TGraph*)gt, "b2,e2", 0, 0)};
        auto func {
            std::make_unique<TF1>("func", [&](double* x, double* p) { return p[0] * spline->Eval(x[0]); }, 0, 180, 1)};
        // Set parameters
        func->SetParameters(1);
        func->SetParName(0, "SF");
        // And fit exp to theo!
        auto res {fExp->Fit(func.get(), "RSQN", "", xmin, xmax)};
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
    std::cout << "·· Fitted in range [" << fFitRange.first << ", " << fFitRange.second << "] deg" << '\n';
    for(const auto& name : fKeys)
    {
        auto& res {fRes.at(name)};
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

TCanvas* Angular::Comparator::Draw(const TString& title, bool withSF, double offset)
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
    bool isTheo {false};
    for(const auto& name : fKeys)
    {
        auto& g {fFit[name]};
        if(!g)
        {
            std::cout << "Comparator::Draw(): Fit() has not been called before Draw()!" << '\n';
            std::cout << "-> plotting theoretical model!" << '\n';
            // continue;
            g = fTheo[name];
            withSF = false;
            isTheo = true;
        }
        g->SetLineWidth(2);
        TString desc {name};
        if(withSF)
            desc += TString::Format(" #Rightarrow SF = %.2f", fRes[name]->Value(0));
        leg->AddEntry(g, desc, "l");
        mg->Add(g, "c");
    }
    // Plot
    // Canvas counter
    static int cCompIdx {};
    auto* c {new TCanvas {TString::Format("cComp%d", cCompIdx), (title.Length()) ? title : "Angular::Comparator"}};
    cCompIdx++;
    if(isTheo)
    {
        // Range in Y
        auto ymin {TMath::MinElement(fExp->GetN(), fExp->GetY())};
        auto ymax {TMath::MaxElement(fExp->GetN(), fExp->GetY())};
        mg->SetMinimum(ymin * 0.8);
        mg->SetMaximum(ymax * 1.2);
    }
    mg->Draw("a plc pmc");
    if(fFitRange.first > 0 && fFitRange.second > 0)
        mg->GetXaxis()->SetLimits(fFitRange.first - offset, fFitRange.second + offset);
    else
    {
        // Get range from exp
        auto xmin {TMath::MinElement(fExp->GetN(), fExp->GetX())};
        auto xmax {TMath::MaxElement(fExp->GetN(), fExp->GetX())};
        mg->GetXaxis()->SetLimits(xmin - offset, xmax + offset);
    }
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
    for(const auto& name : fKeys)
    {
        auto& g {fTheo[name]};
        g->SetLineWidth(2);
        leg->AddEntry(g, name.c_str(), "l");
        mg->Add(g, "c");
    }
    // Plot
    // Canvas counter
    static int cTheoIdx {};
    auto* c {new TCanvas {TString::Format("cTheo%d", cTheoIdx), "Theo models"}};
    cTheoIdx++;
    mg->Draw("a plc");
    leg->Draw();
    return c;
}

TCanvas*
Angular::Comparator::ScaleToExp(const std::string& model, double theoSF, TGraphErrors* gcounts, TEfficiency* eff)
{
    if(!fTheo.count(model))
        throw std::runtime_error("Comparator::ScaleToExp(): could not locate model " + model);
    // Clone to avoid any change in model
    auto theo {(TGraphErrors*)fTheo[model]->Clone()};
    // Scale to theoSF
    theo->Scale(theoSF);
    // Create TMultiGraph to draw it
    auto* mg {new TMultiGraph};
    mg->SetTitle((model + ";#theta_{CM} [#circ];d#sigma / d#Omega [mb /sr]").c_str());
    fTheo[model]->SetLineWidth(2);
    mg->Add(fTheo[model]);
    theo->SetLineStyle(2);
    theo->SetLineWidth(2);
    mg->Add(theo);
    // Theoretical function
    double theoMin {theo->GetPointX(0)};
    double theoMax {theo->GetPointX(theo->GetN() - 1)};
    auto* funcTheo {new TF1 {"funcTheo", [=](double* x, double* p) { return theo->Eval(x[0], nullptr, "S"); }, theoMin,
                             theoMax, 1}};
    // Counts function
    double cMin {gcounts->GetPointX(0)};
    double cMax {gcounts->GetPointX(gcounts->GetN() - 1)};
    auto* funcC {new TF1 {"funcC",
                          [=](double* x, double* p)
                          {
                              if(x[0] < cMin)
                                  return 0.;
                              if(x[0] > cMax)
                                  return 0.;
                              return gcounts->Eval(x[0], nullptr, "S");
                          },
                          cMin, cMax, 1}};
    // Divide function!
    auto* funcDiv {new TF1 {"funcDiv", [=](double* x, double* p) { return funcC->Eval(x[0]) / funcTheo->Eval(x[0]); },
                            theoMin, theoMax, 1}};

    // Plot!
    static int cScaleIdx {};
    auto* c {new TCanvas {TString::Format("cScale%d", cScaleIdx), "Scaling to exp"}};
    cScaleIdx++;
    c->DivideSquare(2);
    c->cd(1);
    mg->Draw("al plc");
    c->cd(2);
    if(eff)
        eff->Draw("apl");
    funcDiv->SetNormalized(true);
    funcDiv->SetTitle("Ratio exp to theo;#theta_{CM} [#circ];N_{exp} / (SF #upoint N_{theo})");
    funcDiv->SetNpx(1000);
    funcDiv->Draw((eff) ? "same" : "");
    return c;
}

TCanvas* Angular::Comparator::QuotientPerPoint()
{
    // Create quotients
    auto* mg {new TMultiGraph};
    mg->SetTitle("Quotient exp / theo;#theta_{CM} [#circ];Quotient exp / theo");
    for(const auto& [name, gtheo] : fTheo)
    {
        auto* gq {new TGraphErrors};
        for(int p = 0; p < fExp->GetN(); p++)
        {
            auto xexp {fExp->GetPointX(p)};
            auto yexp {fExp->GetPointY(p)};
            // Eval theoretical
            auto ytheo {gtheo->Eval(xexp, nullptr, "S")};
            auto uq {1. / ytheo * fExp->GetErrorY(p)};
            // Fill with quotient
            gq->SetPoint(p, xexp, yexp / ytheo);
            gq->SetPointError(p, 0, uq);
        }
        gq->SetTitle(name.c_str());
        gq->SetMarkerStyle(24);
        gq->SetLineWidth(2);
        mg->Add(gq);
    }
    static int cQIdx {};
    auto* c {new TCanvas {TString::Format("cQuot%d", cQIdx), "Quotient theo to exp"}};
    cQIdx++;
    mg->Draw("apl plc pmc");
    c->BuildLegend();
    return c;
}
