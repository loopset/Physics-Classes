#include "AngIntervals.h"

#include "ROOT/RDF/HistoModels.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "THStack.h"
#include "TMath.h"
#include "TString.h"
#include "TVirtualPad.h"

#include "FitUtils.h"

#include <mutex>

Angular::Intervals::Intervals(double xmin, double xmax, const ROOT::RDF::TH1DModel& model, double step, int nps)
{
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
    int idx {};
    for(const auto& [min, max] : fRanges)
    {
        int nbins {model.fNbinsX};
        double hxmax {model.fXUp};
        double hxmin {model.fXLow};
        fHs.push_back(new TH1D {TString::Format("hCM%d", idx),
                                TString::Format("#theta_{CM} #in [%.2f, %.2f)#circ;E_{x} [MeV];Counts per %.0f keV",
                                                min, max, (hxmax - hxmin) / nbins * 1e3),
                                nbins, hxmin, hxmax});
        fHs.back()->SetDirectory(
            nullptr); // not added to gROOT list of histograms (avoid warning if creating Intevals in for loop)
        fOmegas.push_back(ComputeSolidAngle(min, max));
        // Add also histograms for PS
        for(int ps = 0; ps < nps; ps++)
        {
            fHsPS[ps].push_back(
                new TH1D {TString::Format("hPS%dCM%d", ps, idx),
                          TString::Format("PS % d #theta_{CM} #in [%.2f, %.2f)#circ;E_{x} [MeV]", ps, min, max),
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

void Angular::Intervals::TreatPS(int nsmooth, double scale)
{
    for(auto& hps : fHsPS)
    {
        for(int i = 0; i < hps.size(); i++)
        {
            Fitters::TreatPS(fHs[i], hps[i], nsmooth, scale);
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
                hps[i]->SetLineStyle(2);
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
        if(withFuncs)
            gPad->BuildLegend();
    }
    return c;
}
