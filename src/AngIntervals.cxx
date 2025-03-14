#include "AngIntervals.h"

#include "ROOT/RDF/HistoModels.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TMath.h"
#include "TString.h"
#include "TVirtualPad.h"

#include "AngGlobals.h"
#include "FitUtils.h"

#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>

Angular::Intervals::Intervals(double xmin, double xmax, const ROOT::RDF::TH1DModel& model, double step, int nps)
{
    // Avoid adding to gROOT list of objects
    TH1::AddDirectory(false);
    if(step < 0)
        fRanges.push_back({xmin, xmax});
    else
    {
        for(double theta = xmin; theta < xmax; theta += step)
            fRanges.push_back({theta, theta + step});
    }
    // Init vector for phase spaces
    fHsPS.resize(nps);
    // Init histograms and solid angle values
    TString label {gIsLab ? "#theta_{Lab}" : "#theta_{CM}"};
    int idx {};
    for(const auto& [min, max] : fRanges)
    {
        int nbins {model.fNbinsX};
        double hxmax {model.fXUp};
        double hxmin {model.fXLow};
        fHs.push_back(new TH1D {TString::Format("hIvs%d", idx),
                                TString::Format("%s #in [%.2f, %.2f)#circ;E_{x} [MeV];Counts per %.0f keV",
                                                label.Data(), min, max, (hxmax - hxmin) / nbins * 1e3),
                                nbins, hxmin, hxmax});
        fHs.back()->SetDirectory(
            nullptr); // not added to gROOT list of histograms (avoid warning if creating Intevals in for loop)
        fOmegas.push_back(ComputeSolidAngle(min, max));
        // Add also histograms for PS
        for(int ps = 0; ps < nps; ps++)
        {
            fHsPS[ps].push_back(
                new TH1D {TString::Format("hPS%dIvs%d", ps, idx),
                          TString::Format("PS %d %s #in [%.2f, %.2f)#circ;E_{x} [MeV]", ps, label.Data(), min, max),
                          model.fNbinsX, model.fXLow, model.fXUp});
        }
        idx++;
    }
}

double Angular::Intervals::ComputeSolidAngle(double min, double max)
{
    return TMath::TwoPi() * (TMath::Cos(min * TMath::DegToRad()) - TMath::Cos(max * TMath::DegToRad()));
}

void Angular::Intervals::Fill(double thetaCM, double Ex)
{
    for(int i = 0; i < fRanges.size(); i++)
    {
        auto min {fRanges[i].first};
        auto max {fRanges[i].second};
        if(min <= thetaCM && thetaCM < max)
        {
            {
                std::lock_guard<std::mutex> lock {fMutex};
                fHs[i]->Fill(Ex);
            }
            break;
        }
    }
}

void Angular::Intervals::FillPS(int idx, double thetaCM, double Ex, double weight)
{
    for(int i = 0; i < fRanges.size(); i++)
    {
        auto min {fRanges[i].first};
        auto max {fRanges[i].second};
        if(min <= thetaCM && thetaCM < max)
        {
            {
                std::lock_guard<std::mutex> lock {fMutex};
                fHsPS[idx][i]->Fill(Ex, weight);
            }
            break;
        }
    }
}

void Angular::Intervals::TreatPS(int nsmooth, double scale, const std::set<int>& which)
{
    int idx {}; // idx of PS
    for(auto& hps : fHsPS)
    {
        // Disable treat ps for ps not in set
        // By default set is empty so run for all
        if(which.size())
            if(which.find(idx) == which.end())
                continue;
        for(int i = 0; i < hps.size(); i++) // i = interval
        {
            Fitters::TreatPS(fHs[i], hps[i], nsmooth, scale);
        }
        idx++;
    }
}

void Angular::Intervals::FitPS(int ps, int iv, const std::string& pol)
{
    auto& h {fHsPS[ps][iv]};
    // Find fitting range
    auto bmax {h->GetMaximumBin()};
    auto max {h->GetBinContent(bmax)};
    auto thresh {0.04};
    auto blow {h->FindFirstBinAbove(max * thresh)};
    auto bup {h->FindLastBinAbove(max * thresh)};
    fHsPS[ps][iv]->Fit(pol.c_str(), "0QM", "", h->GetBinCenter(blow), h->GetBinCenter(bup));
    // To view fit better
    fHsPS[ps][iv]->SetLineWidth(2);
}

void Angular::Intervals::FitPS(const std::string& pol)
{
    for(int p = 0; p < fHsPS.size(); p++)
        for(int i = 0; i < fHsPS[p].size(); i++)
            FitPS(p, i, pol);
}

void Angular::Intervals::ReplacePSWithFit()
{
    for(int p = 0; p < fHsPS.size(); p++)
    {
        for(int i = 0; i < fHsPS[p].size(); i++)
        {
            // Get function
            auto& h {fHsPS[p][i]};
            TF1* func {};
            if(h->GetListOfFunctions()->GetSize() == 1)
            {
                func = (TF1*)h->GetListOfFunctions()->First()->Clone();
                h->Reset();
                h->Add(func);
                delete func;
            }
        }
    }
}

TCanvas* Angular::Intervals::Draw(const TString& title) const
{
    static int cIvsIdx {};
    auto* c {new TCanvas {TString::Format("cIvs%d", cIvsIdx), (title.Length()) ? title : "Angular::Intervals"}};
    cIvsIdx++;
    c->DivideSquare(fHs.size());
    for(int i = 0; i < fHs.size(); i++)
    {
        c->cd(i + 1);
        bool withFuncs {};
        if(fHsPS.size() == 0)
            fHs[i]->Draw();
        else
        {
            auto* hs {new THStack};
            hs->SetTitle(TString::Format("%s;E_{x} [MeV]", fHs[i]->GetTitle()));
            hs->Add(fHs[i]);
            for(auto& hps : fHsPS)
            {
                hps[i]->SetLineStyle(1);
                hs->Add(hps[i], "hist");
            }
            hs->Draw("nostack");
        }
        // Draw fitted functions if any
        int color {1};
        for(auto* f : *(fHs[i]->GetListOfFunctions()))
        {
            if(color == 10)
                color++;
            if(f)
            {
                ((TF1*)f)->SetLineColor(color);
                f->Draw("same");
                withFuncs = true;
                color++;
            }
        }
        for(int p = 0; p < fHsPS.size(); p++)
        {
            for(auto* f : *(fHsPS[p][i]->GetListOfFunctions()))
            {
                if(color == 10)
                    color++;
                if(f)
                {
                    ((TF1*)f)->SetLineColor(color);
                    f->Draw("same");
                    withFuncs = true;
                    color++;
                }
            }
        }
        if(withFuncs)
            gPad->BuildLegend();
    }
    return c;
}

void Angular::Intervals::Read(const std::string& file)
{
    auto f {std::make_unique<TFile>(file.c_str())};
    auto* intervals {f->Get<Intervals>("Intervals")};
    if(!intervals)
        throw std::runtime_error("Intervals::Read(): cannot locate key Intervals in file");
    fRanges = intervals->GetRanges();
    fOmegas = intervals->GetOmegas();
    fHs = intervals->GetHistos();
    fHsPS = intervals->GetHistosPS();
    delete intervals;
}

void Angular::Intervals::Write(const std::string& file)
{
    auto f {std::make_unique<TFile>(file.c_str(), "recreate")};
    f->WriteObject(this, "Intervals");
}
