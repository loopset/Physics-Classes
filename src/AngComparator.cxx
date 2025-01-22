#include "AngComparator.h"

#include "Rtypes.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TSpline.h"
#include "TString.h"
#include "TVirtualPad.h"

#include "PhysColors.h"
#include "PhysExperiment.h"
#include "PhysSF.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

void Angular::Comparator::Add(const std::string& name, const std::string& file, int lc, int ls, int lw)
{
    // So far, assume a twofnr file format
    fTheo[name] = ReadTwoFNR(file);
    // Set title
    fTheo[name]->SetTitle(name.c_str());
    // And add key
    fKeys.push_back(name);
    // And line style
    fStyles[name] = {lc, ls, lw};
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
        fRes[name] = *res.Get();
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
        std::cout << "   SF    : " << res.Value(0) << " +/- " << res.Error(0) << '\n';
        auto chi2 {res.Chi2()};
        auto ndf {res.Ndf()};
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
    return l;
}

TVirtualPad* Angular::Comparator::Draw(const TString& title, bool logy, bool withSF, double offset, TVirtualPad* pad)
{
    // Draw all using a TMultiGraph
    fMulti = new TMultiGraph;
    fMulti->SetTitle((fName + " ;#theta_{CM} [#circ];d#sigma / d#Omega [mb / sr]").c_str());
    // Create a legend
    auto* leg {BuildLegend()};
    // 1-> Add experimental
    fExp->SetMarkerStyle(25);
    fExp->SetLineWidth(2);
    leg->AddEntry(fExp, "Exp", "pe");
    fMulti->Add(fExp, "p");
    // 2-> Add all fitted
    // In case nothing was fitted, add theoretical lines
    bool isTheo {fFit.size() == 0};
    bool haveCustomColor {};
    for(const auto& name : fKeys)
    {
        TGraphErrors* g {};
        if(isTheo)
        {
            g = fTheo[name];
            withSF = false;
        }
        else
            g = fFit[name];

        // Set styles!
        auto [lc, ls, lw] {fStyles[name]};
        if(lc != -1)
        {
            g->SetLineColor(lc);
            haveCustomColor = true;
        }
        g->SetLineStyle(ls);
        g->SetLineWidth(lw);
        TString desc {name};
        if(withSF)
            desc += TString::Format(" #Rightarrow SF = %.2f", fRes[name].Value(0));
        leg->AddEntry(g, desc, "l");
        fMulti->Add(g, "c");
    }
    // Plot
    if(!pad)
    {
        // Canvas counter
        static int cCompIdx {};
        pad = new TCanvas {TString::Format("cComp%d", cCompIdx), (title.Length()) ? title : "Angular::Comparator"};
        cCompIdx++;
    }
    if(logy)
        pad->SetLogy();
    fMulti->Draw((haveCustomColor) ? "a" : "a plc pmc");
    // Set ranges
    if(fFitRange.first > 0 && fFitRange.second > 0)
        fMulti->GetXaxis()->SetLimits(fFitRange.first - offset, fFitRange.second + offset);
    if(isTheo || logy)
    {
        if(isTheo)
        {
            // Range in X
            auto xmin {TMath::MinElement(fExp->GetN(), fExp->GetX())};
            auto xmax {TMath::MaxElement(fExp->GetN(), fExp->GetX())};
            fMulti->GetXaxis()->SetLimits(xmin - offset, xmax + offset);
        }
        // Range in Y
        auto ymin {TMath::MinElement(fExp->GetN(), fExp->GetY())};
        auto ymax {TMath::MaxElement(fExp->GetN(), fExp->GetY())};
        double scaleMin {(logy) ? 0.1 : 0.6};
        double scaleMax {(logy) ? 10. : 1.3};
        fMulti->SetMinimum(ymin * scaleMin);
        fMulti->SetMaximum(ymax * scaleMax);
    }
    pad->Update();
    leg->Draw();
    // Somehow GetXaxis sets the selected pad and this causes
    // all posterior DrawClones to be drawn on this canvas
    // reset it!
    gROOT->SetSelectedPad(nullptr);
    return pad;
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

TCanvas* Angular::Comparator::ScaleToExp(const std::string& model, PhysUtils::Experiment* exp, TGraphErrors* gcounts,
                                         TEfficiency* teff, double SF)
{
    if(!fTheo.count(model))
        throw std::runtime_error("Comparator::ScaleToExp(): could not locate model " + model);
    // Get SF
    double theoSF {1};
    if(SF != -1)
        theoSF = SF;
    else
    {
        auto fitSF {GetSF(model)};
        if(fitSF != -1)
            theoSF = fitSF;
    }
    // Compute binning
    double bwidth {gcounts->GetPointX(1) - gcounts->GetPointX(0)}; // deg
    auto* gEff {new TGraphErrors};
    gEff->SetTitle("Computed");

    // And manually set contents
    for(int p = 0; p < gcounts->GetN(); p++)
    {
        auto theta {gcounts->GetPointX(p)};
        auto thetaMin {theta - bwidth / 2};
        auto thetaMax {theta + bwidth / 2};
        auto Omega {TMath::TwoPi() *
                    (TMath::Cos(thetaMin * TMath::DegToRad()) - TMath::Cos(thetaMax * TMath::DegToRad()))};
        auto xstheo {fTheo[model]->Eval(theta)};
        auto N {gcounts->GetPointY(p)};
        auto denom {exp->GetNb() * exp->GetNt() * theoSF * xstheo * Omega};
        auto eff {N / denom};
        eff *= 1e27; // properly convert mb units to cm2 in xs * Nt
        // Compute uncertainty
        auto coeffN {theoSF / (exp->GetNb() * exp->GetNt() * xstheo * Omega)};
        auto uN {TMath::Sqrt(N)};
        auto coeffNb {-theoSF * N / (exp->GetNt() * Omega * xstheo) / TMath::Power(exp->GetNb(), 2)};
        auto uNb {exp->GetUNb()};
        auto ueff {TMath::Sqrt(coeffN * coeffN * uN * uN + coeffNb * coeffNb * uNb * uNb)};
        ueff *= 1e27;
        // Fill
        gEff->SetPoint(p, theta, eff);
        gEff->SetPointError(p, 0, ueff);
    }

    // Draw
    // Counter
    static int cScaleIdx {};
    auto* c {new TCanvas {TString::Format("cScale%d", cScaleIdx), "Scaling to exp"}};
    // Increase counter
    cScaleIdx++;
    // Create multigraph to draw all together
    TGraphAsymmErrors* gteff {};
    if(teff)
    {
        gteff = teff->CreateGraph();
        gteff->SetTitle("Simulated");
    }

    auto* mg {new TMultiGraph};
    mg->SetTitle(TString::Format("Eff. for %s;#theta_{CM} [#circ];#epsilon", model.c_str()));
    gEff->SetMarkerStyle(24);
    gEff->SetLineWidth(2);
    gEff->SetLineColor(kMagenta);
    mg->Add(gEff);
    if(gteff)
    {
        gteff->SetLineWidth(2);
        gteff->SetLineColor(8);
        gteff->SetFillStyle(0);
        mg->Add(gteff);
    }
    mg->Draw("apl");
    c->BuildLegend();

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

double Angular::Comparator::GetSF(const std::string& model)
{
    if(fRes.count(model))
        return fRes[model].Parameter(0);
    else
    {
        std::cout << BOLDRED << "Angular::Comparator::GetSF(): cannot get SF because model " << model
                  << " hasn't been fitted yet" << RESET << '\n';
        return -1;
    }
}

double Angular::Comparator::GetuSF(const std::string& model)
{
    if(fRes.count(model))
        return fRes[model].Error(0);
    else
    {
        std::cout << BOLDRED << "Angular::Comparator::GetuSF(): cannot get u(SF) because model " << model
                  << " hasn't been fitted yet" << RESET << '\n';
        return -1;
    }
}

void Angular::Comparator::Write(const std::string& file)
{
    auto f {std::make_unique<TFile>(file.c_str(), "recreate")};
    // Save in two vectors (std::map not supported wo dict in root file)
    std::vector<std::string> names;
    std::vector<PhysUtils::SpectroscopicFactor> sfs;
    for(const auto& [model, res] : fRes)
    {
        names.push_back(model);
        sfs.emplace_back(res.Parameter(0), res.ParError(0), res.Chi2() / res.Ndf(), res.Ndf());
    }
    f->WriteObject(&names, "Models");
    f->WriteObject(&sfs, "SFs");
    f->WriteObject(fMulti, "MultiGraph");
}
