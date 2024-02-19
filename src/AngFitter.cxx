#include "AngFitter.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TString.h"

#include "FitData.h"
#include "FitModel.h"
#include "FitPlotter.h"
#include "FitRunner.h"
#include "PhysColors.h"

#include <ios>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

Angular::Fitter::Fitter(const std::vector<TH1D*>& data, double exmin, double exmax)
{
    fHistos = data;
    for(auto& h : data)
        AddData(h, exmin, exmax);
}

Angular::Fitter::Fitter(Intervals& ivs, double exmin, double exmax)
{
    for(auto& h : ivs.GetHistos())
        AddData(h, exmin, exmax);
    fHistos = ivs.GetHistos();
}

void Angular::Fitter::AddData(TH1D* data, double xmin, double xmax)
{
    fData.push_back(Fitters::Data {*data, xmin, xmax});
}

void Angular::Fitter::AddModels()
{
    // count unique values
    std::set<std::string> counts;
    for(const auto& name : fParNames)
    {
        auto str {name.substr(0, name.find_first_of('_'))};
        counts.insert(str);
    }
    // Count number of functions
    int ngauss {};
    int nvoigt {};
    int nps {};
    int ncte {};
    for(const auto& name : counts)
    {
        TString str {name};
        if(str.Contains("g"))
            ngauss++;
        if(str.Contains("v"))
            nvoigt++;
        if(str.Contains("ps"))
            nps++;
        if(str.Contains("cte"))
            ncte++;
    }
    // Assert nps matches size of passed data
    if(fPS.size() != nps)
        throw std::runtime_error("Fitter::AddModels(): parsed nps does not match size of std::vector<TH1D>");
    // Assert ncte = 0 or 1
    if(ncte > 1)
        throw std::runtime_error("Fitter::AddModels(): parsed ncte is greater than 1");
    for(int i = 0; i < fData.size(); i++)
        fModels.push_back(Fitters::Model {ngauss, nvoigt, fPS, (bool)ncte});
    // Print
    std::cout << BOLDGREEN << "-> NGauss : " << ngauss << RESET << '\n';
    std::cout << BOLDGREEN << "-> NVoigt : " << nvoigt << RESET << '\n';
    std::cout << BOLDGREEN << "-> NPS    : " << nps << RESET << '\n';
    std::cout << BOLDGREEN << "-> Cte    ? " << std::boolalpha << (bool)ncte << RESET << '\n';
}

void Angular::Fitter::Configure(const std::string& file, const std::vector<TH1D>& ps)
{
    std::cout << BOLDGREEN << "---- Angular::Fitter ----" << '\n';
    std::cout << "-> Config : " << file << RESET << '\n';
    // Open
    auto f {std::make_unique<TFile>(file.c_str())};
    // Read par names
    std::vector<std::string>* names {};
    f->GetObject("ParNames", names);
    fParNames = *names;
    // And fit result
    auto* res {f->Get<TFitResult>("FitResult")};
    fGlobalFit = *res;
    // Set PS data also
    fPS = ps;
    // Init models
    AddModels();
}

void Angular::Fitter::ConfigRunner(Fitters::Runner& runner)
{
    // Fitter
    auto& f {runner.GetFitter()};
    // Number of parameters
    auto npars {runner.GetObjective().NDim()};
    for(unsigned int p = 0; p < npars; p++)
    {
        // 1-> Set parameter name
        f.Config().ParSettings(p).SetName(fParNames[p]);
        // 2-> Initial value
        f.Config().ParSettings(p).SetValue(fGlobalFit.Value(p));
        // 3-> Bounds
        double min {};
        double max {};
        fGlobalFit.ParameterBounds(p, min, max);
        if(TString {fParNames[p]}.Contains("_Amp")) // only set bounds for amp
            f.Config().ParSettings(p).SetLimits(min, max);
        // 4-> Fix: only amp is free
        if(!(TString {fParNames[p]}.Contains("_Amp")))
            f.Config().ParSettings(p).Fix();
    }
}

void Angular::Fitter::Run()
{
    for(int i = 0; i < fData.size(); i++)
    {
        // Initialize runnner
        Fitters::Runner runner {fData[i], fModels[i]};
        // Pass integral opts to fitter
        runner.GetObjective().SetUseDivisions(fUseDivisions);
        runner.GetObjective().SetUseIntegral(fUseIntegral);
        // Config it
        ConfigRunner(runner);
        // Fit!
        runner.Fit();
        fRes.push_back(runner.GetFitResult());
    }
    // Implicitly compute integrals
    ComputeIntegrals();
}

TCanvas* Angular::Fitter::Draw()
{
    auto* c {new TCanvas {"cFitter", "Angular::Fitter"}};
    c->DivideSquare(fData.size());
    for(int i = 0; i < fData.size(); i++)
    {
        c->cd(i + 1);
        Fitters::Plotter plot {&fData[i], &fModels[i], &fRes[i]};
        auto* gfit {plot.GetGlobalFit()};
        auto hfits {plot.GetIndividualHists()};
        // Draw with fitted range
        fHistos[i]->GetXaxis()->SetRangeUser(fData[i].GetXLow(), fData[i].GetXUp());
        fHistos[i]->Draw();
        gfit->Draw("same");
        // Create stack
        auto* hs {new THStack};
        for(auto& [kye, h] : hfits)
        {
            h->SetLineWidth(2);
            hs->Add(h);
        }
        hs->Draw("plc nostack same");
    }
    return c;
}

void Angular::Fitter::ComputeIntegrals(int nsigma)
{
    // Clear vectors
    fIgCounts.clear();
    fSumCounts.clear();
    // Run for each interval
    for(int i = 0; i < fData.size(); i++)
        DoCounts(i, nsigma);
}

void Angular::Fitter::DoCounts(unsigned int iv, int nsigma)
{
    // Unpack parameters
    auto pack {fModels[iv].UnpackParameters(fRes[iv].GetParams())};
    auto gaus = pack[0];
    auto voigt = pack[1];
    // Range of numerical integration
    auto xmin {fData[iv].GetXLow()};
    auto xmax {fData[iv].GetXUp()};
    // And integrate!
    // 1-> Gauss
    int idx {};
    for(const auto& pars : gaus)
    {
        std::string key {"g" + std::to_string(idx)};
        auto* f {new TF1 {"gaus", "gaus", xmin, xmax}};
        f->SetParameters(pars.data());
        fIgCounts[key].push_back(f->Integral(xmin, xmax) / fData[iv].GetBinWidth());
        // By sum
        CountsBySum(key, iv, nsigma, f);
        delete f;
        idx++;
    }
    // 2-> Voigt
    idx = 0;
    for(const auto& pars : voigt)
    {

        std::string key {"v" + std::to_string(idx)};
        auto f = new TF1("voigt", "[0] * TMath::Voigt(x - [1], [2], [3])", xmin, xmax);
        f->SetParameters(pars.data());
        fIgCounts[key].push_back(f->Integral(xmin, xmax) / fData[iv].GetBinWidth());
        // By sum
        CountsBySum(key, iv, nsigma, f);
        delete f;
        idx++;
    }
}

void Angular::Fitter::CountsBySum(const std::string& key, unsigned int iv, int nsigma, TF1* f)
{
    // Sum must be done with fitted histogram! Not with raw data
    int nbins {static_cast<int>((fData[iv].GetXUp() - fData[iv].GetXLow()) / fData[iv].GetBinWidth())};
    auto* h {new TH1D("htemp", "h for sum", nbins, fData[iv].GetXLow(), fData[iv].GetXUp())};
    for(int bin = 1; bin <= h->GetNbinsX(); bin++)
    {
        auto x {h->GetXaxis()->GetBinCenter(bin)};
        auto y {f->Eval(x)};
        h->SetBinContent(bin, y);
    }
    // Set scale according to nsigma
    double scale {};
    if(nsigma == 1)
        scale = 0.68;
    else if(nsigma == 2)
        scale = 0.95;
    else
        throw std::runtime_error("Angular::Fitter::CountsBySum(): no nsigma correction factor implemented");
    // And do integral
    auto mean {f->GetParameter(1)};
    auto sigma {f->GetParameter(2)};
    // Get intervals of integration, dependent on nsigma around mean
    double low {mean - nsigma * sigma};
    auto binLow {h->FindBin(low)};
    double up {mean + nsigma * sigma};
    auto binUp {h->FindBin(up)};
    // Integral(bin1, bin2) just sums bin contents
    auto integral {h->Integral(binLow, binUp)};
    // Scale to full data
    fSumCounts[key].push_back(integral / scale);
    delete h;
}

TCanvas* Angular::Fitter::DrawCounts()
{
    auto* c {new TCanvas {"cCounts", "Angular::Fitter counts"}};
    auto* leg {new TLegend {0.2, 0.2}};
    leg->SetTextSize(0.035);
    leg->SetFillStyle(0);
    leg->SetNColumns(2);
    // Integral
    auto* mg {new TMultiGraph};
    mg->SetTitle("Counts per #theta_{CM} interval;#theta_{CM} index;Counts");
    int idx {0};
    for(const auto& [key, counts] : fIgCounts)
    {
        auto* g {new TGraphErrors};
        for(int iv = 0; iv < counts.size(); iv++)
        {
            g->SetPoint(iv, iv, counts[iv]);
            g->SetPointError(iv, 0, TMath::Sqrt(counts[iv]));
        }
        // Style
        g->SetLineWidth(2);
        // Legend
        leg->AddEntry(g, key.c_str(), "lp");
        mg->Add(g);
        idx++;
    }
    // Counts
    auto* mc {new TMultiGraph};
    mc->SetTitle("Counts per #theta_{CM} interval;#theta_{CM} index;Counts");
    idx = 0;
    for(const auto& [key, counts] : fSumCounts)
    {
        auto* g {new TGraphErrors};
        for(int iv = 0; iv < counts.size(); iv++)
        {
            g->SetPoint(iv, iv, counts[iv]);
            g->SetPointError(iv, 0, TMath::Sqrt(counts[iv]));
        }
        // Style
        g->SetLineWidth(2);
        g->SetLineStyle(2);
        mc->Add(g);
        idx++;
    }
    // Draw
    mc->Draw("alp plc pmc pfc");
    mg->Draw("lp plc pmc pfc");
    leg->Draw();
    return c;
}
