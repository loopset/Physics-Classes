#include "FitUtils.h"

#include "Rtypes.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TFitResult.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"

#include "FitPlotter.h"
#include "FitRunner.h"
#include "PhysColors.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

void Fitters::TreatPS(TH1D* hEx, TH1D* hPS)
{
    // 1-> Smooth it
    hPS->Smooth(20);
    // 2-> Scale it to have a reasonable height
    auto intEx {hEx->Integral()};
    auto intPS {hPS->Integral()};
    double factor {0.1};
    hPS->Scale(factor * intEx / intPS);
}

void Fitters::DrawGlobalFit(TGraph* g, const std::unordered_map<std::string, TH1D*>& hs, TLegend* leg,
                            const std::unordered_map<std::string, std::string>& labels)
{
    // Draw in current gPad
    g->Draw("same");
    leg->AddEntry(g, "Global fit", "l");
    auto* stack {new THStack};
    // Set fill styles
    std::vector<int> fills {3345, 3354, 3344};
    int idx {};
    for(auto& [key, h] : hs)
    {
        TString rstr {key};
        bool isPeak {rstr.Contains("g") || rstr.Contains("v")};
        h->SetLineWidth(2);
        h->SetFillStyle(0);
        if(idx < fills.size() && isPeak)
            h->SetFillStyle(fills[idx]);
        leg->AddEntry(h, (labels.count(key) ? labels.at(key) : key).c_str(), (isPeak) ? "lf" : "l");
        stack->Add(h);
        idx++;
    }
    stack->Draw("nostack plc pfc same");
    if(hs.size() > 4)
        leg->SetNColumns((int)hs.size() / 4 + 1);
    leg->Draw();
}

void Fitters::RunFit(TH1D* h, double exmin, double exmax, Fitters::Model& model, const Fitters::Runner::Init& initial,
                     const Fitters::Runner::Bounds& bounds, const Fitters::Runner::Fixed& fixed,
                     const std::string& outfile, const std::string& title,
                     const std::unordered_map<std::string, std::string>& labels, const Fitters::Runner::Step& steps)
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
    if(fixed.size() > 0)
        runner.SetStep(steps);
    // Run
    runner.Fit();
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
    // Configure global Ex histogram
    h->SetTitle("");
    h->SetStats(false);
    h->GetXaxis()->SetRangeUser(exmin, exmax);
    h->SetLineWidth(2);
    auto* clone = h->DrawClone("e");
    // Create a legend
    auto* leg {new TLegend {0.5, 0.6, 0.9, 0.9}};
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(clone, "Experimental", "le");
    // Draw all the other histograms and legend
    DrawGlobalFit(gfit, hfits, leg, labels);
    // End :)
    std::cout << BOLDCYAN << "++++++++++++++++++++++++++++++" << RESET << '\n';
}

Fitters::Runner::Init Fitters::ReadInit(const std::string& name)
{
    auto file {std::make_unique<TFile>(name.c_str())};
    std::cout << BOLDCYAN << "Fitters::ReadInit() from file : " << name << RESET << '\n';
    // Read parameter names
    auto* names {file->Get<std::vector<std::string>>("ParNames")};
    auto* res {file->Get<TFitResult>("FitResult")};

    Fitters::Runner::Init ret;
    for(int i = 0; i < names->size(); i++)
    {
        // Split
        auto it {names->at(i).find_first_of("_")};
        auto key {names->at(i).substr(0, it)};
        auto val {res->Parameter(i)};
        ret[key].push_back(val);
    }
    return ret;
}
