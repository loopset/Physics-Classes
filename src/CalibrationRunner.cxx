#ifndef CalibrationRunner_cxx
#define CalibrationRunner_cxx
#include "CalibrationRunner.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TROOT.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TVirtualPad.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

void Calibration::Runner::ApplyRange()
{
    fData->GetXaxis()->SetRangeUser(fRange.first, fRange.second);
}

bool Calibration::Runner::DoIt()
{
    // Find pedestal
    FindPedestal();
    // Run precalibration
    if(!DoPreCalibration())
        return false;
    // And fine calibration
    DoFineCalibration();
    // Final plots
    DoFinalPlots();
    if(fDebug)
        Debug();
    return true;
}

void Calibration::Runner::FindPedestal()
{
    // Pedestal is the maximum in the histogram
    // from [0, fRange.first]
    auto [mean, amp] {GetAmpMeanInRange(fData, 0, fRange.first)};
    // Temporaly set range properly
    fData->GetXaxis()->SetRangeUser(0, fRange.first);
    // Fit to gaussian
    fPedestals["ped"] = std::make_unique<TF1>("ped", "gaus", 0, mean + fPreGaussWidth);
    fData->Fit(fPedestals["ped"].get(), (fDebug ? fFitOptsDebug : fFitOpts).c_str());
    // Reset range! (as in constructor)
    fData->GetXaxis()->SetRangeUser(fRange.first, fRange.second);
}

std::pair<double, double> Calibration::Runner::GetPedestal()
{
    auto& f {fPedestals["ped"]};
    return {f->GetParameter(1), f->GetParameter(2)};
}

std::vector<std::pair<double, double>> Calibration::Runner::FilterPeaks(const TSpectrum& spe)
{
    // Sort them in increasing order of X
    std::vector<std::pair<double, double>> ret;
    auto* legacyX {spe.GetPositionX()};
    auto* legacyY {spe.GetPositionY()};
    for(int i = 0; i < spe.GetNPeaks(); i++)
        ret.push_back({legacyX[i], legacyY[i]});
    // Sort
    std::sort(ret.begin(), ret.end(), [](const auto& l, const auto& r) { return l.first < r.first; });
    return ret;
}

void Calibration::Runner::FillHistPre()
{
    // Initialize histogram after preliminar calibration using major peaks
    fHistPre = std::make_shared<TH1D>("hPre", "Precalibration hist;E_{pre} [MeV];Counts",
                                      (int)(fHistOpts.first / fHistOpts.second), 0, fHistOpts.first);
    TH1D* data {(fRawData) ? fRawData : fData};
    for(int bin = 1; bin <= data->GetNbinsX(); bin++)
    {
        // Get X value in channel
        auto x {data->GetBinCenter(bin)};
        // Energy of that x
        auto energy {fCalibPre->Eval(x)};
        // Get content of bin and fill energy content times
        auto content {data->GetBinContent(bin)};
        for(int i = 0; i < content; i++)
            fHistPre->Fill(energy);
    }
    // And set new range
    auto min {fCalibPre->Eval(fRange.first)};
    auto max {fCalibPre->Eval(fRange.second)};
    fHistPre->GetXaxis()->SetRangeUser(min, max);
}

bool Calibration::Runner::DoPreCalibration()
{
    // Use TSpectrum class to look for peaks in histogram!
    TSpectrum spe {5};
    spe.Search(fData, fSpeSigma, "nodraw", fSpeThresh);
    auto peaks {FilterPeaks(spe)};
    if(peaks.size() < 3)
    {
        std::cout << "Runner:DoPreCalibration(): TSpectrum # peaks < 3, cannot proceed. Change range of histogram"
                  << '\n';
        return false;
    }
    // Fit peaks to gaussians
    int idx {0};
    for(const auto& s : fSource->GetLabels())
    {
        fGaussPre[s] = std::make_shared<TF1>(("pre" + s).c_str(), "gaus", peaks[idx].first - fPreGaussWidth,
                                             peaks[idx].first + fPreGaussWidth);
        fData->Fit(fGaussPre[s].get(), (fDebug) ? fFitOptsDebug.c_str() : fFitOpts.c_str());
        idx++;
    }
    // And now get precalibration using major peaks
    fGraphPre = std::make_shared<TGraphErrors>();
    fGraphPre->SetTitle("Precalibration graph;Channel;E_{major} [MeV]");
    fGraphPre->SetMarkerStyle(24);
    auto major {fSource->GetMajorPeaks()};
    for(const auto& s : fSource->GetLabels())
    {
        auto mean {fGaussPre[s]->GetParameter("Mean")};
        auto umean {fGaussPre[s]->GetParError(1)};
        fGraphPre->SetPoint(fGraphPre->GetN(), mean, major[s]);
        fGraphPre->SetPointError(fGraphPre->GetN() - 1, umean, 0);
    }
    fCalibPre = std::make_shared<TF1>("precalib", "pol1", fRange.first, fRange.second);
    fCalibPre->SetParameters(-10, 0.001); // Initial guess needed!
    // Fit options
    std::string opts {};
    if(fDebug)
        opts = fFitOptsGraphDebug;
    else
        opts = fFitOptsGraph;
    if(fDisableXErrors)
        opts += "EX0";
    fGraphPre->Fit(fCalibPre.get(), opts.c_str());
    // Fill new histogram
    FillHistPre();
    return true;
}

std::pair<double, double> Calibration::Runner::GetAmpMeanInRange(TH1D* h, double min, double max)
{
    // Enhanced determination of initial parameter for each peak
    auto* clone {(TH1D*)h->Clone("hclone")};
    clone->GetXaxis()->SetRangeUser(min, max);
    int bin {clone->GetMaximumBin()};
    double x {clone->GetXaxis()->GetBinCenter(bin)};
    double y {clone->GetBinContent(bin)};
    delete clone;
    // Some corrections
    if(x < min)
        x = min;
    if(x > max)
        x = max;
    return {x, y};
}

void Calibration::Runner::AddSatellite(SatelliteCont& map, const GaussCont& funcs)
{
    std::vector<int> colors {6, 8, 9, 46};
    std::vector<int> lines {2, 5, 7, 3};
    for(const auto& [name, f] : funcs)
    {
        for(int p = 0; p < fSatelliteStr[name].size(); p++)
        {
            map[name].push_back(
                std::make_shared<TF1>(TString::Format("sat%d%s", p, name.c_str()), fSatelliteStr[name][p].c_str()));
            // Range from fine gauss
            double min {};
            double max {};
            f->GetRange(min, max);
            map[name][p]->SetRange(min, max);
            map[name][p]->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
            // Set styles
            map[name][p]->SetLineColor(colors[p]);
            map[name][p]->SetLineStyle(lines[p]);
        }
    }
}

void Calibration::Runner::DoFineCalibration()
{
    // Get all data
    auto [energies, sigmas, brs] {fSource->GetComponents()};
    // Declare functions
    for(int s = 0; s < fSource->GetLabels().size(); s++)
    {
        auto name {fSource->GetLabels()[s]};
        TString sumf {};
        // Add all peaks (main + satellite)
        // Get maximum energy and br
        double maxE {*std::max_element(energies[s].begin(), energies[s].end())};
        double maxBR {*std::max_element(brs[s].begin(), brs[s].end())};
        for(int p = 0; p < energies[s].size(); p++)
        {
            auto satellite {TString::Format(" + [0] * %f * TMath::Exp(-0.5 * TMath::Power((x - ([1] - %f)) / [2], 2))",
                                            brs[s][p] / maxBR, TMath::Abs(energies[s][p] - maxE))};
            fSatelliteStr[name].push_back(satellite.Data());
            sumf += satellite;
        }
        fGaussFine[name] = std::make_shared<TF1>(("fine" + name).c_str(), sumf);
        // Increase number of plotting points
        fGaussFine[name]->SetNpx(1000);
        // Set range form source
        auto [min, max] {fSource->GetLimits(name)};
        fGaussFine[name]->SetRange(min, max);
        // Set some limits on parameters to avoid wrong fits
        fGaussFine[name]->SetParLimits(0, 0, 1e8);   // amp
        fGaussFine[name]->SetParLimits(1, min, max); // mean
        fGaussFine[name]->SetParLimits(2, 0.01, 5);  // sigma c [0, 5] MeV
        // And now initial parameters!
        auto [mean, amp] {GetAmpMeanInRange(fHistPre.get(), min, max)};
        fGaussFine[name]->SetParameters(amp, mean, sigmas[s].back());
    }
    // Fit!
    for(auto& [name, f] : fGaussFine)
        fHistPre->Fit(f.get(), (fDebug ? fFitOptsDebug : fFitOpts).c_str());

    // And get finer calibration
    fGraphFine = std::make_shared<TGraphErrors>();
    fGraphFine->SetTitle("Fine calibration graph;E_{pre} [MeV];E_{fine} [MeV]");
    fGraphFine->SetMarkerStyle(25);
    auto major {fSource->GetMajorPeaks()};
    for(const auto& [name, f] : fGaussFine)
    {
        auto mean {f->GetParameter(1)};
        auto umean {f->GetParError(1)};
        fGraphFine->SetPoint(fGraphFine->GetN(), major[name], mean);
        fGraphFine->SetPointError(fGraphFine->GetN() - 1, umean, 0);
    }
    fCalibFine = std::make_shared<TF1>("finecalib", "pol1", 0, fHistOpts.first);
    fCalibFine->SetParameters(0., 1);
    fGraphFine->Fit(fCalibFine.get(), (fDebug ? fFitOptsGraphDebug : fFitOptsGraph).c_str());

    // Add satellite peaks
    AddSatellite(fFineSat, fGaussFine);

    // Init final calibration
    InitFinalCalib();
    // and histogram
    FillHistFinal();
}

void Calibration::Runner::InitFinalCalib()
{
    fCalibFinal =
        std::make_shared<TF1>("finalcalib", "pol1", fData->GetXaxis()->GetXmin(), fData->GetXaxis()->GetXmax());
    double a {fCalibFine->GetParameter(0) + fCalibFine->GetParameter(1) * fCalibPre->GetParameter(0)};
    double b {fCalibFine->GetParameter(1) * fCalibPre->GetParameter(1)};
    fCalibFinal->SetParameters(a, b);
    fCalibFinal->SetTitle("Final calibration;Channel;E [MeV]");
}

void Calibration::Runner::FillHistFinal()
{
    // Initialize histogram after preliminar calibration using major peaks
    fHistFinal = std::make_shared<TH1D>("hFinal", "Calibrated hist;E [MeV];Counts",
                                        (int)(fHistOpts.first / fHistOpts.second), 0, fHistOpts.first);
    TH1D* data {(fRawData) ? fRawData : fData};
    for(int bin = 1; bin <= data->GetNbinsX(); bin++)
    {
        // Get X value (if rebinned, bin != channel!!!)
        auto x {data->GetBinCenter(bin)};
        // Energy
        auto energy {fCalibFinal->Eval(x)};
        // Fill energy content times
        auto content {data->GetBinContent(bin)};
        for(int i = 0; i < content; i++)
            fHistFinal->Fill(energy);
    }
    // And set new range
    auto min {fCalibFinal->Eval(fRange.first)};
    auto max {fCalibFinal->Eval(fRange.second)};
    fHistFinal->GetXaxis()->SetRangeUser(min, max);
}

void Calibration::Runner::DoFinalPlots()
{
    // Refit the fully-calibrated histogram to get nice plots for thesis
    auto [energies, sigmas, brs] {fSource->GetComponents()};
    for(int s = 0; s < fSource->GetLabels().size(); s++)
    {
        auto name {fSource->GetLabels()[s]};
        TString sumf {};
        // Add all peaks (main + satellite)
        for(int p = 0; p < energies[s].size(); p++)
            sumf += fSatelliteStr[name][p];
        fGaussFinal[name] = std::make_shared<TF1>(("final" + name).c_str(), sumf);
        fGaussFinal[name]->SetNpx(1000);
        // Set range form source
        auto [min, max] {fSource->GetLimits(name)};
        fGaussFinal[name]->SetRange(min, max);
        // Set some limits on parameters to avoid wrong fits
        fGaussFinal[name]->SetParLimits(0, 0, 1e8);   // amp
        fGaussFinal[name]->SetParLimits(1, min, max); // mean
        fGaussFinal[name]->SetParLimits(2, 0.01, 5);  // sigma c [0, 5] MeV
        // And now initial parameters!
        auto [mean, amp] {GetAmpMeanInRange(fHistFinal.get(), min, max)};
        fGaussFinal[name]->SetParameters(amp, mean, sigmas[s].back());
    }
    // Fit!
    for(auto& [name, f] : fGaussFinal)
        fHistFinal->Fit(f.get(), fFitOpts.c_str());
    // Add satellites
    AddSatellite(fFinalSat, fGaussFinal);
}

std::pair<double, double> Calibration::Runner::GetParameters()
{
    if(!fCalibFinal)
        return {};
    else
        return {fCalibFinal->GetParameter(0), fCalibFinal->GetParameter(1)};
}

void Calibration::Runner::PrintRes() const
{
    std::cout << "---- Calibration::Runner: Resolutions ----" << '\n';
    for(const auto& [name, f] : fGaussFinal)
        std::cout << "-> " << name << " sigma : " << f->GetParameter(2) * 1e3 << " +/- " << f->GetParError(2) * 1e3
                  << " keV" << '\n';
    std::cout << "------------------------------" << '\n';
}

void Calibration::Runner::Debug(TCanvas* c) const
{
    TCanvas* canv {};
    bool toClone {false};
    if(c)
    {
        canv = c;
        toClone = true;
    }
    else
    {
        auto* inner {(TCanvas*)gROOT->GetListOfCanvases()->FindObject("cRunner")};
        if(inner)
            canv = inner;
        else
            canv = new TCanvas("cRunner", "Calibration::Runner canvas");
    }
    // workaround to allow plotting of clones once WaitPrimitive has been called:
    // reset the selected pad to allow the usage of gPad in DrawClone
    gROOT->SetSelectedPad(nullptr);
    canv->Clear();
    canv->DivideSquare(6);
    canv->cd(1);
    fData->Draw();
    DrawAll(fData->GetListOfFunctions());

    // Precalibration
    if(fGraphPre)
    {
        canv->cd(2);
        if(toClone)
            fGraphPre->DrawClone("ap");
        else
            fGraphPre->Draw("ap");
        DrawAll(fGraphPre->GetListOfFunctions(), toClone);
    }
    if(fHistPre)
    {
        canv->cd(3);
        if(toClone)
            fHistPre->DrawClone();
        else
            fHistPre->Draw();
        DrawAll(fHistPre->GetListOfFunctions(), toClone);
        // Draw all satellites
        DrawSat(fFineSat);
    }
    // Fine calibration
    if(fGraphFine)
    {
        canv->cd(4);
        if(toClone)
            fGraphFine->DrawClone("apl");
        else
            fGraphFine->Draw("apl");
        DrawAll(fGraphFine->GetListOfFunctions(), toClone);
    }
    // Final
    if(fCalibFinal)
    {
        canv->cd(5);
        if(toClone)
            fHistFinal->DrawClone();
        else
            fHistFinal->Draw();
        DrawAll(fHistFinal->GetListOfFunctions(), toClone);
        DrawSat(fFinalSat);
    }
    canv->cd();
    canv->Update();
    // Wait
    if(!toClone)
    {
        canv->WaitPrimitive("dummy", "");
        canv->Update();
    }
    // Just for safety, reset again
    gROOT->SetSelectedPad(nullptr);
}

#endif // !CalibrationRunner_cxx
