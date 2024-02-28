#include "AngIntervals.h"

#include "ROOT/RDF/HistoModels.hxx"

#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"

#include <mutex>

Angular::Intervals::Intervals(double xmin, double xmax, const ROOT::RDF::TH1DModel& model, double step)
{
    if(step < 0)
        fRanges.push_back({xmin, xmax});
    else
    {
        for(double theta = xmin; theta < xmax; theta += step)
            fRanges.push_back({theta, theta + step});
    }
    // Init histograms and solid angle values
    int idx {};
    for(const auto& [min, max] : fRanges)
    {
        fHs.push_back(new TH1D {TString::Format("hCM%d", idx),
                                TString::Format("#theta_{CM} #in [%.2f, %.2f)#circ;E_{x} [MeV]", min, max),
                                model.fNbinsX, model.fXLow, model.fXUp});
        fHs.back()->SetDirectory(
            nullptr); // not added to gROOT list of histograms (avoid warning if creating Intevals in for loop)
        fOmegas.push_back(ComputeSolidAngle(min, max));
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

TCanvas* Angular::Intervals::Draw(const TString& title) const
{
    static int cIvsIdx {};
    auto* c {new TCanvas {TString::Format("cIvs%d", cIvsIdx), (title.Length()) ? title : "Angular::Intervals"}};
    cIvsIdx++;
    c->DivideSquare(fHs.size());
    for(int i = 0; i < fHs.size(); i++)
    {
        c->cd(i + 1);
        fHs[i]->Draw();
    }
    return c;
}
