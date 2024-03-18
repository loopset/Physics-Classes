#include "FitUtils.h"

#include "TCanvas.h"
#include "THStack.h"

#include "FitPlotter.h"
#include "PhysColors.h"

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

void Fitters::DrawGlobalFit(TGraph* g, const std::unordered_map<std::string, TH1D*>& hs)
{
    // Draw in current gPad
    g->Draw("same");
    auto* stack {new THStack};
    for(auto& [key, h] : hs)
    {
        h->SetLineWidth(2);
        stack->Add(h);
    }
    stack->Draw("nostack plc same");
}

void Fitters::RunFit(TH1D* h, double exmin, double exmax, Fitters::Model& model, const Fitters::Runner::Init& initial,
                     const Fitters::Runner::Bounds& bounds, const Fitters::Runner::Fixed& fixed,
                     const std::string& outfile, const std::string& title)
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
    h->GetXaxis()->SetRangeUser(exmin, exmax);
    h->SetLineWidth(2);
    h->DrawClone("e");
    DrawGlobalFit(gfit, hfits);
    // End :)
    std::cout << BOLDCYAN << "++++++++++++++++++++++++++++++" << RESET << '\n';
}
