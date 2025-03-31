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

#include "AngGlobals.h"
#include "PhysColors.h"
#include "PhysExperiment.h"
#include "PhysSF.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

Angular::Comparator::Comparator(const std::string& name, TGraphErrors* exp) : fName(name)
{
    if(exp)
        fExp = (TGraphErrors*)exp->Clone();
}

void Angular::Comparator::Add(const std::string& name, const std::string& file, int lc, int ls, int lw)
{
    // Read any two column part of a text file
    fTheo[name] = ReadFile(file);
    // Set title
    fTheo[name]->SetTitle(name.c_str());
    // And add key
    fKeys.push_back(name);
    // And line style
    fStyles[name] = {lc, ls, lw};
}

void Angular::Comparator::Replace(const std::string& name, TGraphErrors* gnew)
{
    if(fTheo.count(name))
        delete fTheo[name];
    // Clone to avoid modifications of passed graph
    auto* clone {(TGraphErrors*)gnew->Clone()};
    fTheo[name] = ProcessTheo(clone);
    // And delete clone once processed
    delete clone;
    // Set title
    fTheo[name]->SetTitle(name.c_str());
    // And add key in case is was not there
    auto it {std::find(fKeys.begin(), fKeys.end(), name)};
    if(it == fKeys.end())
        fKeys.push_back(name);
    // And default line styles
    if(!fStyles.count(name))
        fStyles[name] = {-1, 1, 3};
}

TGraphErrors* Angular::Comparator::ReadFile(const std::string& file)
{
    TGraphErrors theo {file.c_str(), "%lg %lg"};
    return ProcessTheo(&theo);
}

TGraphErrors* Angular::Comparator::ProcessTheo(TGraphErrors* theo)
{
    if(!fExp)
        return new TGraphErrors(*theo);
    // INFO: 24/01/2025. Theoretical graph is preprocessed
    //  so as it has the same binning as the experimental
    //  and the bin content is the INTEGRAL in that bin width averaged by its width

    // Read the content
    // Convert x axis to rad
    theo->Scale(TMath::DegToRad(), "x");
    // Build function to integrate in bin!
    TF1 func {"func", [&](double* x, double* p) { return theo->Eval(x[0], nullptr, "S"); }, 0, TMath::TwoPi(), 0};
    // Minimum and maximum for theoretical graph
    auto tMin {theo->GetPointX(0)};
    auto tMax {theo->GetPointX(theo->GetN() - 1)};
    // Min and max of experimental xs in rad units
    auto eMin {fExp->GetPointX(0) * TMath::DegToRad()};
    auto eMax {fExp->GetPointX(fExp->GetN() - 1) * TMath::DegToRad()};
    double eBW {};
    if(fExp->GetN() > 1)
        eBW = (fExp->GetPointX(1) - fExp->GetPointX(0));
    else
        eBW = 1;
    eBW *= TMath::DegToRad();
    // Get starting bin centre
    double start {};
    for(double x = eMin; x > tMin; x -= eBW)
        start = x;
    // And create theoretical graph with the same binning!
    auto* ret {new TGraphErrors};
    for(double x = start; x < tMax; x += eBW)
    {
        auto low {x - eBW / 2};
        auto up {x + eBW / 2};
        auto integral {func.Integral(low, up)};
        ret->AddPoint(x, integral / eBW);
    }
    // And convert back X axis to deg units
    ret->Scale(TMath::RadToDeg(), "x");
    return ret;
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
    // Check if exp data is empty; in that case, continue
    if(!fExp->GetN())
        return;
    if(!fTheo.size())
        return;
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
        fRes[name] = TFitResult(*res.Get());
    }
    Print();
}

void Angular::Comparator::Print() const
{
    if(!fExp->GetN())
        return;
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

TVirtualPad*
Angular::Comparator::Draw(const TString& title, bool logy, bool withSF, double offset, TVirtualPad* pad, bool withChi)
{
    // Draw all using a TMultiGraph
    fMulti = new TMultiGraph;
    fMulti->SetTitle((fName + ";#theta_{" + (gIsLab ? "Lab" : "CM") + "} [#circ];d#sigma / d#Omega [mb / sr]").c_str());
    // Create a legend
    auto* leg {BuildLegend()};
    // 1-> Add experimental
    if(fExp->GetN())
    {
        fExp->SetMarkerStyle(25);
        fExp->SetLineWidth(2);
        leg->AddEntry(fExp, "Exp", "pe");
        fMulti->Add(fExp, "p");
    }
    // 2-> Add all fitted
    // In case nothing was fitted, add theoretical lines
    bool isTheo {fFit.size() == 0 || fExp->GetN() == 0};
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
        if(fRes.count(name))
        {
            if(withSF)
                desc += TString::Format(" #Rightarrow SF = %.2f", fRes[name].Value(0));
            if(withChi)
                desc += TString::Format(" #cbar #chi^{2}_{#nu} = %.2f", fRes[name].Chi2() / fRes[name].Ndf());
        }
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
    {
        auto xmin {fFitRange.first - offset};
        fMulti->GetXaxis()->SetLimits((xmin < 0) ? 0 : xmin, fFitRange.second + offset);
    }
    if((isTheo || logy) && fExp->GetN())
    {
        if(isTheo)
        {
            // Range in X
            auto xmin {TMath::MinElement(fExp->GetN(), fExp->GetX())};
            auto xmax {TMath::MaxElement(fExp->GetN(), fExp->GetX())};
            xmin -= offset;
            xmax += offset;
            fMulti->GetXaxis()->SetLimits((xmin < 0) ? 0 : xmin, xmax);
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
    mg->SetTitle(
        (fName + " models;#theta_{" + (gIsLab ? "Lab" : "CM") + "} [#circ];d#sigma / d#Omega [mb / sr]").c_str());
    // Legend
    auto* leg {BuildLegend()};
    // Custom color?
    bool haveCustomColor {false};
    // Add graphs
    for(const auto& name : fKeys)
    {
        auto& g {fTheo[name]};
        // Set styles!
        auto [lc, ls, lw] {fStyles[name]};
        if(lc != -1)
        {
            g->SetLineColor(lc);
            haveCustomColor = true;
        }
        g->SetLineStyle(ls);
        g->SetLineWidth(lw);
        leg->AddEntry(g, name.c_str(), "l");
        mg->Add(g, "c");
    }
    // Plot
    // Canvas counter
    static int cTheoIdx {};
    auto* c {new TCanvas {TString::Format("cTheo%d", cTheoIdx), "Theo models"}};
    cTheoIdx++;
    mg->SetMinimum(1.e-3);
    mg->Draw((haveCustomColor) ? "a" : "a plc pmc");
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
    double uTheoSF {0};
    if(SF != -1)
        theoSF = SF;
    else
    {
        auto fitSF {GetSF(model)};
        if(fitSF != -1)
        {
            theoSF = fitSF;
            uTheoSF = GetuSF(model);
        }
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
        auto coeffN {1. / (exp->GetNb() * exp->GetNt() * theoSF * xstheo * Omega)};
        auto uN {TMath::Sqrt(N)};
        auto coeffNb {-N / (exp->GetNt() * Omega * theoSF * xstheo) / TMath::Power(exp->GetNb(), 2)};
        auto uNb {exp->GetUNb()};
        auto coeffSF {-N / (exp->GetNt() * exp->GetNb() * xstheo * Omega) / TMath::Power(theoSF, 2)};
        auto ueff {TMath::Sqrt(coeffN * coeffN * uN * uN + coeffNb * coeffNb * uNb * uNb +
                               coeffSF * coeffSF * uTheoSF * uTheoSF)};
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
    mg->SetTitle(TString::Format("Eff. for %s;#theta_{%s} [#circ];#epsilon", model.c_str(), gIsLab ? "Lab" : "CM"));
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
    mg->SetTitle(TString::Format("%s;#theta_{%s} [#circ];Quotient exp / theo", fName.c_str(), gIsLab ? "Lab" : "CM"));
    for(const auto& key : fKeys)
    {
        auto& gtheo {fTheo[key]};
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
        gq->SetTitle(key.c_str());
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

void Angular::Comparator::Write(const std::string& key, const std::string& file)
{
    bool weOwnFile {};
    TFile* f {};
    if(file.length())
    {
        f = new TFile {file.c_str(), "recreate"};
        weOwnFile = true;
    }
    else
        f = gFile;
    // Create collection
    PhysUtils::SFCollection sfcol;
    for(const auto& [model, res] : fRes)
        sfcol.Add(model, {res.Parameter(0), res.ParError(0), res.Chi2() / res.Ndf(), (int)res.Ndf()});
    // Write!
    f->cd();
    f->WriteObject(&sfcol, (key + "_sfs").c_str());
    if(fMulti)
        f->WriteObject(fMulti, (key + "_mg").c_str());
    if(weOwnFile)
    {
        f->Close();
        delete f;
    }
}
