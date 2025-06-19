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
#include "TRegexp.h"
#include "TString.h"

#include "Math/IntegratorOptions.h"

#include "AngGlobals.h"
#include "FitData.h"
#include "FitModel.h"
#include "FitPlotter.h"
#include "FitRunner.h"
#include "PhysColors.h"

#include <algorithm>
#include <cmath>
#include <ios>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

void Angular::Fitter::AddData(double exmin, double exmax)
{
    fHistos = fIvs->GetHistos();
    if(!fIgnorePS)
        fHistosPS = fIvs->GetHistosPS();
    for(auto& h : fHistos)
        fData.push_back(Fitters::Data {*h, exmin, exmax});
    std::cout << BOLDGREEN << "-> Range  : [" << exmin << ", " << exmax << "] MeV" << RESET << '\n';
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
    if(!fIgnorePS && (fHistosPS.size() != nps))
        throw std::runtime_error(
            "Angular::Fitter::AddModels(): parsed nps does not match size of fHistosPS from Angular::Intervals!");
    // Assert ncte = 0 or 1
    if(ncte > 1)
        throw std::runtime_error("Angular::Fitter::AddModels(): parsed ncte is greater than 1");
    for(int i = 0; i < fData.size(); i++)
    {
        // Build PS vector
        std::vector<TH1D> vps;
        for(auto& hps : fHistosPS)
            vps.push_back(*hps[i]);
        fModels.push_back(Fitters::Model {ngauss, nvoigt, vps, (bool)ncte});
    }
    // Print
    std::cout << BOLDGREEN << "-> NGauss : " << ngauss << RESET << '\n';
    std::cout << BOLDGREEN << "-> NVoigt : " << nvoigt << RESET << '\n';
    std::cout << BOLDGREEN << "-> NPS    : " << nps << RESET << '\n';
    if(fIgnorePS)
        std::cout << BOLDRED << "    but set to ignore them!" << RESET << '\n';
    std::cout << BOLDGREEN << "-> Cte    ? " << std::boolalpha << (bool)ncte << RESET << '\n';
}

void Angular::Fitter::Configure(const std::string& file)
{
    std::cout << BOLDGREEN << "---- Angular::Fitter ----" << '\n';
    std::cout << "-> Config : " << file << '\n';
    std::cout << "-> FreeMean ? " << std::boolalpha << fAllowFreeMean << RESET << '\n';
    if(fAllowFreeMean)
        std::cout << BOLDGREEN << "-> MeanRange : " << fFreeMeanRange << " MeV" << RESET << '\n';
    std::cout << BOLDGREEN << "-> FreeSigma? " << std::boolalpha << fAllowFreeSigma << RESET << '\n';
    if(fAllowFreeSigma)
        std::cout << BOLDGREEN << "-> SigmaRange: " << fFreeSigmaRange << " MeV" << RESET << '\n';
    // Open
    auto f {std::make_unique<TFile>(file.c_str())};
    // Read par names
    std::vector<std::string>* names {};
    f->GetObject("ParNames", names);
    fParNames = *names;
    // And fit result
    auto* res {f->Get<TFitResult>("FitResult")};
    fGlobalFit = *res;
    // Init data in range!
    auto range {f->Get<std::pair<double, double>>("FitRange")};
    if(!range)
        throw std::runtime_error("Fitter::Configure(): could not find FitRange in file " + file);
    // Determine which range to use
    if(fManualRange.first != -1 && fManualRange.second != -1)
    {
        std::cout << BOLDRED << "Angular::Fitter::Configure(): using manual!" << RESET << '\n';
        AddData(fManualRange.first, fManualRange.second);
    }
    else
        AddData(range->first, range->second);
    // Init models (after initializing data, mandatory)
    AddModels();
}

void Angular::Fitter::ConfigRunner(int iv, Fitters::Runner& runner)
{
    // Fitter
    auto& f {runner.GetFitter()};
    // Number of parameters
    auto npars {runner.GetObjective().NDim()};
    for(unsigned int p = 0; p < npars; p++)
    {
        TString str {fParNames[p]};

        // 1-> Set parameter name
        f.Config().ParSettings(p).SetName(str.Data());
        // Extract information
        auto sep {str.Index("_")};     // Example below for "g1_Mean" input str
        TString funcIdx {str(0, sep)}; // g1
        auto num {funcIdx.First("0123456789")};
        TString nameFunc {funcIdx(0, num)};                           // g
        int idxFunc {TString(funcIdx(num, funcIdx.Length())).Atoi()}; // 1
        TString namePar {str(sep + 1, str.Length())};                 // Mean

        // 2-> Initial value
        auto value {fGlobalFit.Value(p)};
        f.Config().ParSettings(p).SetValue(value);
        // 3-> Bounds
        double min {};
        double max {};
        fGlobalFit.ParameterBounds(p, min, max);
        if(str.Contains("_Amp")) // only set bounds for amp
            f.Config().ParSettings(p).SetLimits(min, max);

        // 4-> Fix parameters
        if(!(str.Contains("_Amp")))
            f.Config().ParSettings(p).Fix();

        // If requested, allow a slight variation in mean
        if(fAllowFreeMean && str.Contains("_Mean"))
        {
            // If specified particular states for which allow free mean, apply only to those!
            if(fWhichFreeMean.size() &&
               std::find(fWhichFreeMean.begin(), fWhichFreeMean.end(), funcIdx) == fWhichFreeMean.end())
                continue;
            f.Config().ParSettings(p).Release();
            f.Config().ParSettings(p).SetLimits(value - fFreeMeanRange, value + fFreeMeanRange);
        }

        // If requested, allow a slight variation in sigma
        if(fAllowFreeSigma && str.Contains("_Sigma"))
        {
            // If specified free sigma per state, only for those!
            if(fWhichFreeSigma.size() &&
               std::find(fWhichFreeSigma.begin(), fWhichFreeSigma.end(), funcIdx) == fWhichFreeSigma.end())
                continue;
            f.Config().ParSettings(p).Release();
            f.Config().ParSettings(p).SetLimits(TMath::Max(0.1, value - fFreeSigmaRange), value + fFreeSigmaRange);
        }
        // If requested, allow a slight variation in gamma
        if(fAllowFreeGamma && str.Contains("_Lg"))
        {
            // If specified free sigma per state, only for those!
            if(fWhichFreeGamma.size() &&
               std::find(fWhichFreeGamma.begin(), fWhichFreeGamma.end(), funcIdx) == fWhichFreeGamma.end())
                continue;
            f.Config().ParSettings(p).Release();
            f.Config().ParSettings(p).SetLimits(TMath::Max(0., value - fFreeGammaRange), value + fFreeGammaRange);
        }

        // If requested, fix amplitude for PS in this interval
        if(str.Contains("ps") && str.Contains("_Amp"))
        {
            // Check if that PS has given amplitudes
            if(fPSFixAmps.count(idxFunc))
            {
                double fixVal {1};
                // And then if we have amplitude for that interval
                // if not, fallback to 1
                if(fPSFixAmps[idxFunc].size() > iv)
                    fixVal = fPSFixAmps[idxFunc][iv];
                f.Config().ParSettings(p).SetValue(fixVal);
                f.Config().ParSettings(p).Fix();
                std::cout << BOLDGREEN << "Fitter::ConfigRunner: fixing " << str << " in iv " << iv << " at " << fixVal
                          << RESET << '\n';
            }
        }
    }
}

void Angular::Fitter::Run(bool print)
{
    for(int i = 0; i < fData.size(); i++)
    {
        // Initialize runnner
        Fitters::Runner runner {fData[i], fModels[i]};
        // Pass integral opts to fitter
        runner.GetObjective().SetUseDivisions(fUseDivisions);
        runner.GetObjective().SetUseIntegral(fUseIntegral);
        // Config it
        ConfigRunner(i, runner);
        // Fit!
        runner.Fit(print, Angular::GetUseHessErrors());
        fRes.push_back(runner.GetFitResult());
    }
    // Implicitly compute integrals
    ComputeIntegrals();
    // And fill fitted histograms
    FillResHistos();
}

void Angular::Fitter::FillResHistos()
{
    fResIvs.clear();
    fResGlobal.clear();
    fResHistos.clear();
    fResNames.clear();
    for(int i = 0; i < fData.size(); i++)
    {
        Fitters::Plotter plot {&fData[i], &fModels[i], &fRes[i]};
        // Centre
        if(fIvs)
            fResIvs.push_back(fIvs->GetCenter(i));
        // Global fit
        fResGlobal.push_back(plot.GetGlobalFit());
        // Histograms
        auto hfits {plot.GetIndividualHists()};
        fResHistos.push_back({});
        fResNames.push_back({});
        for(const auto& [key, h] : hfits)
        {
            fResHistos.back().push_back(h);
            fResNames.back().push_back(key);
        }
    }
}

TCanvas* Angular::Fitter::Draw(const TString& title)
{
    static int cFitIdx {};
    auto* c {new TCanvas {TString::Format("cFitter%d", cFitIdx), (title.Length()) ? title : "Angular::Fitter"}};
    cFitIdx++;
    // Size
    auto size {fData.size() ? fData.size() : fResIvs.size()};
    c->DivideSquare(size);
    for(int i = 0; i < size; i++)
    {
        c->cd(i + 1);
        auto* gfit {fResGlobal[i]};
        auto& hfits {fResHistos[i]};
        // Keep drawn ones independent from fResGlobal and fResHistos
        // Fitters::Plotter plot {&fData[i], &fModels[i], &fRes[i]};
        // auto* gfit {plot.GetGlobalFit()};
        // auto hfits {plot.GetIndividualHists()};
        // Draw with fitted range if Data is available
        if(i < fData.size())
            fHistos[i]->GetXaxis()->SetRangeUser(fData[i].GetXLow(), fData[i].GetXUp());
        fHistos[i]->Draw();
        gfit->Draw("same");
        // Create stack
        auto* hs {new THStack};
        for(auto& h : hfits)
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
    // Issue with Voigts: the default
    // AdiptiveSingular integration algorithm always
    // yields errors related with the tolerance
    // Changing to Gauss provides the same results without warning
    if(voigt.size())
        ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
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
    // Create a histogram from the interval data to sum counts
    int nbins {static_cast<int>((fData[iv].GetXUp() - fData[iv].GetXLow()) / fData[iv].GetBinWidth())};
    auto* h {new TH1D("htemp", "h for sum", nbins, fData[iv].GetXLow(), fData[iv].GetXUp())};
    for(int histbin = 1; histbin <= h->GetNbinsX(); histbin++)
    {
        auto x {h->GetXaxis()->GetBinCenter(histbin)};
        // auto y {f->Eval(x)};
        // WARNING: Contents are filled from the raw data, not from the evaluation of the fitted function
        // because in that case of course there is a match!!
        // Must convert to inner binning of each data element
        auto inbin {fData[iv].GetBin(x)};
        auto y {fData[iv].GetY(inbin)};
        h->SetBinContent(histbin, y);
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

Angular::Fitter::CountsIv Angular::Fitter::GetIgCountsFor(const std::string& peak) const
{
    if(fIgCounts.count(peak))
        return fIgCounts.at(peak);
    else
        throw std::invalid_argument("Fitter::GetIgCountsFor(): received not listed peak");
}

TGraphErrors* Angular::Fitter::GetIgCountsGraph(const std::string& peak) const
{
    if(!fIgCounts.count(peak))
        throw std::runtime_error("Fitter::GetIgCountsGraph(): could not find peak " + peak);
    const auto& counts {fIgCounts.at(peak)};
    auto* g {new TGraphErrors};
    for(int iv = 0; iv < counts.size(); iv++)
    {
        double xlabel {static_cast<double>(iv)};
        if(fIvs)
            xlabel = fIvs->GetCenter(iv);
        g->SetPoint(iv, xlabel, counts[iv]);
        g->SetPointError(iv, 0, TMath::Sqrt(counts[iv]));
    }
    // Style
    g->SetLineWidth(2);
    return g;
}

Angular::Fitter::CountsIv Angular::Fitter::GetSumCountsFor(const std::string& peak) const
{
    if(fSumCounts.count(peak))
        return fSumCounts.at(peak);
    else
        throw std::invalid_argument("Fitter::GetSumCountsFor(): received not listed peak");
}

TGraphErrors* Angular::Fitter::GetSumCountsGraph(const std::string& peak) const
{
    if(!fSumCounts.count(peak))
        throw std::runtime_error("Fitter::GetSumCountsGraph(): could not find peak " + peak);
    const auto& counts {fSumCounts.at(peak)};
    auto* g {new TGraphErrors};
    for(int iv = 0; iv < counts.size(); iv++)
    {
        double xlabel {static_cast<double>(iv)};
        if(fIvs)
            xlabel = fIvs->GetCenter(iv);
        g->SetPoint(iv, xlabel, counts[iv]);
        g->SetPointError(iv, 0, TMath::Sqrt(counts[iv]));
    }
    // Style
    g->SetLineWidth(2);
    g->SetLineStyle(2);
    return g;
}

TCanvas* Angular::Fitter::DrawCounts(bool both, const TString& title)
{
    static int cCountsIdx {};
    auto* c {
        new TCanvas {TString::Format("cCounts%d", cCountsIdx), (title.Length()) ? title : "Angular::Fitter counts"}};
    cCountsIdx++;
    auto* leg {new TLegend {0.2, 0.2}};
    leg->SetTextSize(0.035);
    leg->SetFillStyle(0);
    leg->SetNColumns(2);
    // Integral
    auto* mg {new TMultiGraph};
    TString label {gIsLab ? "#theta_{Lab}" : "#theta_{CM}"};
    if(fIvs)
        mg->SetTitle(TString::Format("Counts per %s interval;%s [#circ];Counts", label.Data(), label.Data()));
    else
        mg->SetTitle(TString::Format("Counts per %s interval;%s interval idx;Counts", label.Data(), label.Data()));
    int idx {0};
    for(const auto& [key, counts] : fIgCounts)
    {
        auto* g {GetIgCountsGraph(key)};
        g->SetNameTitle(key.c_str(), key.c_str());
        // Legend
        leg->AddEntry(g, key.c_str(), "lp");
        mg->Add(g);
        idx++;
    }
    // Counts
    auto* mc {new TMultiGraph};
    if(fIvs)
        mc->SetTitle(TString::Format("Counts per %s interval;%s [#circ];Counts", label.Data(), label.Data()));
    else
        mc->SetTitle(TString::Format("Counts per %s interval;%s interval idx;Counts", label.Data(), label.Data()));

    idx = 0;
    for(const auto& [key, counts] : fSumCounts)
    {
        auto* g {GetSumCountsGraph(key)};
        mc->Add(g);
        idx++;
    }
    // Draw
    if(both)
    {
        mc->Draw("alp plc pmc");
        mg->Draw("lp plc pmc");
    }
    else
    {
        mg->Draw("alp plc pmc");
        delete mc;
    }
    leg->Draw();
    return c;
}

void Angular::Fitter::Write(const std::string& file)
{
    // Save as TGraphs
    auto* fout {new TFile {file.c_str(), "recreate"}};
    for(const auto& [peak, _] : fIgCounts)
    {
        auto* g {GetIgCountsGraph(peak)};
        g->SetNameTitle(("g" + peak).c_str(), ("Integral counts for " + peak).c_str());
        g->Write();
    }

    // Save this object
    fout->WriteObject(this, "Fitter");

    fout->Close();
}

void Angular::Fitter::Read(const std::string& file)
{
    auto f {std::make_unique<TFile>(file.c_str())};
    auto fitter {f->Get<Fitter>("Fitter")};
    *this = *fitter;
}
