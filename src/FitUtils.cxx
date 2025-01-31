#include "FitUtils.h"

#include "Rtypes.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TList.h"
#include "TString.h"

#include "FitPlotter.h"
#include "FitRunner.h"
#include "PhysColors.h"

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
