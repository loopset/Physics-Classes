#include "FitUtils.h"

#include "Rtypes.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TList.h"
#include "TRatioPlot.h"
#include "TString.h"
#include "TStyle.h"

#include "FitModel.h"
#include "FitPlotter.h"
#include "FitRunner.h"
#include "PhysColors.h"

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

void Fitters::TreatPS(TH1D* hEx, TH1D* hPS, int nsmooth, double scale)
{
    // 1-> Smooth it
    hPS->Smooth(nsmooth);
    // 2-> Scale it to have a reasonable height
    auto intEx {hEx->Integral()};
    auto intPS {hPS->Integral()};
    if(intPS == 0 || intEx == 0)
        return;
    hPS->Scale(scale * intEx / intPS);
}

void Fitters::FitPS(TH1D* hPS, const std::string& pol, bool replace, bool draw)
{
    // Find fitting range
    auto bmax {hPS->GetMaximumBin()};
    auto max {hPS->GetBinContent(bmax)};
    auto thresh {0.01};
    auto blow {hPS->FindFirstBinAbove(max * thresh)};
    auto bup {hPS->FindLastBinAbove(max * thresh)};
    hPS->Fit(pol.c_str(), "0QM", "", hPS->GetBinCenter(blow), hPS->GetBinCenter(bup));
    auto* func {hPS->GetFunction(pol.c_str())};
    auto* hclone {(TH1D*)hPS->Clone()};
    if(func)
    {
        func->SetLineColor(kMagenta);
        func->ResetBit(TF1::kNotDraw);
        // And replace if requested
        if(replace)
        {
            auto* clone {(TF1*)func->Clone()};
            hPS->Reset();
            hPS->Add(clone);
            delete clone;
        }
    }
    // If draw
    if(draw)
    {
        static int counter {0};
        auto* c {new TCanvas {TString::Format("cFitPS%d", counter), TString::Format("Fit PS %d", counter)}};
        hclone->SetLineColor(kRed);
        hclone->DrawClone();
        hPS->DrawClone("same");
    }
    // deletes
    delete hclone;
}

void Fitters::DrawGlobalFit(TGraph* g, const std::unordered_map<std::string, TH1D*>& hs, TLegend* leg,
                            const std::unordered_map<std::string, std::string>& labels)
{
    // Draw in current gPad
    if(g)
    {
        g->Draw("same");
        leg->AddEntry(g, "Global fit", "l");
    }
    auto* stack {new THStack};
    // Set fill styles
    std::vector<int> fills {3345, 3354, 3344};
    int idx {};
    for(auto& [key, h] : hs)
    {
        TString rstr {key};
        bool isPS {rstr.Contains("ps")};
        h->SetLineWidth(2);
        h->SetFillStyle(0);
        if(idx < fills.size() && isPS)
        {
            h->SetFillStyle(fills[idx]);
            idx++;
        }
        leg->AddEntry(h, (labels.count(key) ? labels.at(key) : key).c_str(), isPS ? "lf" : "l");
        stack->Add(h);
    }
    stack->Draw("nostack plc pfc same");
    if(hs.size() > 4)
        leg->SetNColumns((int)hs.size() / 4 + 1);
    leg->Draw();
}

void Fitters::RunFit(TH1D* h, double exmin, double exmax, Fitters::Model& model, const Fitters::Runner::Init& initial,
                     const Fitters::Runner::Bounds& bounds, const Fitters::Runner::Fixed& fixed,
                     const std::string& outfile, const std::string& title,
                     const std::unordered_map<std::string, std::string>& labels, bool residuals, bool minos)
{
    std::cout << BOLDCYAN << "++++ Global fit " << title << " ++++" << RESET << '\n';
    // Init data
    Fitters::Data data {*h, exmin, exmax};

    // Init runner
    Fitters::Runner runner {data, model};
    runner.GetObjective().SetUseDivisions(true);
    // And initial parameters
    runner.SetInitial(initial);
    runner.SetBounds(bounds);
    runner.SetFixed(fixed);
    // Run
    runner.Fit(true, minos);
    // Save
    runner.Write(outfile);

    // Get fit result
    auto res {runner.GetFitResult()};

    // Plotter
    Fitters::Plotter plot {&data, &model, &res};
    auto* gfit {plot.GetGlobalFit()};
    auto hfits {plot.GetIndividualHists()};

    // Draw
    static int cCount {};
    auto* c {new TCanvas {TString::Format("cGF%d", cCount), title.c_str()}};
    cCount++;
    // Create a legend
    auto* leg {new TLegend {0.5, 0.6, 0.9, 0.9}};
    // Draw depending on mode
    if(!residuals)
    {
        // Configure global Ex histogram
        h->SetTitle("");
        h->SetStats(false);
        h->GetXaxis()->SetRangeUser(exmin, exmax);
        h->SetLineWidth(2);
        auto* clone = h->DrawClone("e");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.045);
        leg->AddEntry(clone, "Experimental", "le");
        // Draw all the other histograms and legend
        DrawGlobalFit(gfit, hfits, leg, labels);
    }
    else
    {
        auto* clone {(TH1D*)h->Clone()};
        clone->SetLineWidth(2);
        clone->GetXaxis()->SetRangeUser(exmin, exmax);
        clone->SetStats(false);
        // Model pointer
        auto* ptr {model.Clone()};
        // Wrap into TF1
        auto* tf1 {new TF1 {"fitfunc", [=](double* x, double* p) { return (*ptr)(x, p); }, exmin, exmax,
                            static_cast<int>(ptr->NPar())}};
        tf1->SetParameters(res.Parameters().data());
        tf1->SetNpx(2000);
        // And append to funtion list
        clone->GetListOfFunctions()->Clear();
        clone->GetListOfFunctions()->Add(tf1);
        // Results pointer just in case...
        auto* rptr {new TFitResult {res}};
        // And ratio plot now!
        auto* ratio {new TRatioPlot {clone, "", rptr}};
        ratio->SetH1DrawOpt("e");    // visualize errors for histogram
        ratio->SetGraphDrawOpt("p"); // same for residuals
        ratio->Draw();
        ratio->GetLowerRefYaxis()->SetTitle("Residuals");
        ratio->GetLowerRefYaxis()->SetLabelSize(0.5 * gStyle->GetLabelSize("Y"));
        // std::cout << "Canvas : " << c << '\n';
        // std::cout << "up : " << ratio->GetUpperPad() << '\n';
        // std::cout << "low : " << ratio->GetLowerPad() << '\n';
        // std::cout << "gPad : " << gPad << '\n';
        // And draw the other things
        ratio->GetUpperPad()->cd();
        DrawGlobalFit(nullptr, hfits, leg, labels);
        c->cd();
    }
    // Save to file
    SaveGlobalFit(outfile, h, gfit, hfits, leg);
    // End :)
    std::cout << BOLDCYAN << "++++++++++++++++++++++++++++++" << RESET << '\n';
}

std::pair<Fitters::Runner::Init, Fitters::Runner::Init> Fitters::ReadInit(const std::string& name)
{
    auto file {std::make_unique<TFile>(name.c_str())};
    if(file->IsZombie())
    {
        std::cout << BOLDRED << "Fitters::ReadInit(): file" << name << " does not exist!" << RESET << '\n';
        return {};
    }
    // Read parameter names
    auto* names {file->Get<std::vector<std::string>>("ParNames")};
    auto* res {file->Get<TFitResult>("FitResult")};

    Fitters::Runner::Init vals, uncs;
    for(int i = 0; i < names->size(); i++)
    {
        // Split
        auto it {names->at(i).find_first_of("_")};
        auto key {names->at(i).substr(0, it)};
        auto val {res->Parameter(i)};
        auto unc {res->Error(i)};
        vals[key].push_back(val);
        uncs[key].push_back(unc);
    }
    return {vals, uncs};
}


void Fitters::SaveGlobalFit(const std::string& file, TH1D* h, TGraph* g,
                            const std::unordered_map<std::string, TH1D*>& hs, TLegend* leg)
{
    auto f {std::make_unique<TFile>(file.c_str(), "update")};
    // Save Ex spectrum
    h->Write("HistoEx");
    // Global fit
    g->Write("GraphGlobal");
    // Individual fits
    std::vector<std::string> keys;
    TList vhs; // Cannot be saved in std::vector
    for(const auto& [key, h] : hs)
    {
        keys.push_back(key);
        vhs.Add(h);
    }
    f->WriteObject(&keys, "NamePeaks");
    f->WriteObject(&vhs, "HistoPeaks");
    f->WriteObject(leg, "Legend");
}

void Fitters::ReadDrawGlobalFit(const std::string& file)
{
    auto f {std::make_unique<TFile>(file.c_str())};
    // Ex
    auto* h {f->Get<TH1D>("HistoEx")};
    h->SetDirectory(nullptr);
    // Global fit
    auto* g {f->Get<TGraph>("GraphGlobal")};
    // Individual fits
    auto* vhs {f->Get<TList>("HistoPeaks")};
    // Names
    auto* names {f->Get<std::vector<std::string>>("NamePeaks")};
    // Draw
    auto* leg {new TLegend {0.5, 0.6, 0.9, 0.9}};
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h, "Experimental", "le");
    leg->AddEntry(g, "Global fit", "l");
    h->Draw("e");
    g->Draw("l");
    auto* stack {new THStack};
    int idx {};
    for(auto* o : *vhs)
    {
        auto* h {dynamic_cast<TH1D*>(o)};
        h->SetDirectory(nullptr);
        leg->AddEntry(h, names->at(idx).c_str(), "l");
        stack->Add(h, "hist");
        idx++;
    }
    stack->Draw("nostack plc same");
    auto size {leg->GetListOfPrimitives()->GetSize()};
    auto ncols {size / 4 + (size % 4 ? 1 : 0)};
    leg->SetNColumns(ncols);
    leg->Draw();
}
