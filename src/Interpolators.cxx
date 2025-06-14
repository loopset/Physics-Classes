#include "Interpolators.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TVirtualPad.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

void Interpolators::Efficiency::AddGraph(const std::string& peak, TEfficiency* eff)
{
    if(fGraph[peak])
    {
        delete fGraph[peak];
        fGraph[peak] = nullptr;
    }
    fGraph[peak] = eff->CreateGraph();
    fGraph[peak]->SetName(peak.c_str());
    fGraph[peak]->SetBit(TGraph::kIsSortedX); // mark X to be already sorted
}

void Interpolators::Efficiency::Add(const std::string& peak, const std::string& file, const std::string& name)
{
    // Read file
    auto f {std::make_unique<TFile>(file.c_str())};
    // Get TEff obj
    auto eff {f->Get<TEfficiency>(name.c_str())};
    if(!eff)
        throw std::runtime_error("Efficiency::Add(): could not read " + name + " key in file " + file);
    // Determine whether is lab from name
    TString str {eff->GetName()};
    str.ToLower();
    if(str.Contains("lab"))
        fIsLab = true;
    // Push to map
    eff->SetTitle(TString::Format("%s;#theta_{%s} [#circ];#epsilon", peak.c_str(), fIsLab ? "Lab" : "CM"));
    fEff[peak] = eff;
    // Create graph
    AddGraph(peak, eff);
}

void Interpolators::Efficiency::Add(const std::string& peak, TEfficiency* eff)
{
    // Determine whether is lab from name
    TString str {eff->GetName()};
    str.ToLower();
    if(str.Contains("lab"))
        fIsLab = true;
    if(fEff[peak])
    {
        delete fEff[peak];
        fEff[peak] = nullptr;
    }
    fEff[peak] = eff;
    fEff[peak]->SetNameTitle(peak.c_str(), peak.c_str());
    AddGraph(peak, eff);
}

double Interpolators::Efficiency::GetMeanEff(const std::string& peak, double min, double max) const
{
    auto* g {fGraph.at(peak)};
    std::vector<double> vals;
    double step {TMath::Abs(max - min) / fMeanDiv};
    for(auto t = min; t <= max; t += step)
        vals.push_back(fGraph.at(peak)->Eval(t, nullptr, "S"));
    return TMath::Mean(vals.begin(), vals.end());
}

double Interpolators::Efficiency::GetPointUEff(const std::string& peak, double thetaCM) const
{
    auto* eff {fEff.at(peak)};
    // Get bin
    auto bin {eff->FindFixBin(thetaCM)};
    // Error as mean of low and up (must be equal in principle)
    return (eff->GetEfficiencyErrorLow(bin) + eff->GetEfficiencyErrorUp(bin)) / 2;
}

void Interpolators::Efficiency::Scale(double factor)
{
    for(auto& [key, eff] : fEff)
    {
        auto* hAll {eff->GetTotalHistogram()};
        auto* hPassed {eff->GetPassedHistogram()};
        auto* hScaled {(TH1D*)hPassed->Clone()};
        hScaled->Reset();
        for(int b = 1; b <= hScaled->GetNbinsX(); b++)
        {
            auto counts {hPassed->GetBinContent(b)};
            counts *= factor;
            hScaled->SetBinContent(b, int(counts));
        }
        auto* neweff {new TEfficiency {*hScaled, *hAll}};
        // Delete auxiliaries
        delete hScaled;
        // Add (this automatically overrides TGraph)
        Add(key, neweff);
    }
}

std::pair<double, double> Interpolators::Efficiency::GetXaxisRange(TMultiGraph* mg)
{
    double xmax {-1111};
    double xmin {1111};
    for(auto* o : *mg->GetListOfGraphs())
    {
        auto* g {(TGraphAsymmErrors*)o};
        auto npoints {g->GetN()};
        for(int i = 0; i < npoints; i++)
        {
            auto x {g->GetPointX(i)};
            auto y {g->GetPointY(i)};
            if(y > 0)
            {
                if(x < xmin)
                    xmin = x;
                break;
            }
        }
        for(int i = npoints - 1; i >= 0; i--)
        {
            auto x {g->GetPointX(i)};
            auto y {g->GetPointY(i)};
            if(y > 0)
            {
                if(x > xmax)
                    xmax = x;
                break;
            }
        }
    }
    return {0.8 * xmin, xmax * 1.2};
}

double Interpolators::Efficiency::GetYaxisRange(TMultiGraph* mg)
{
    double ymax {-1111};
    for(auto* o : *mg->GetListOfGraphs())
    {
        auto* g {(TGraphAsymmErrors*)o};
        auto npoints {g->GetN()};
        for(int i = 0; i < npoints; i++)
        {
            auto y {g->GetPointY(i)};
            if(y > ymax)
                ymax = y;
        }
    }
    return 1.2 * ymax;
}

TVirtualPad* Interpolators::Efficiency::Draw(bool multigraph, const TString& title, TVirtualPad* pad)
{
    TVirtualPad* c {};
    if(pad)
        c = pad;
    else
    {
        static int cEffIdx {};
        auto canvas =
            new TCanvas {TString::Format("cEff%d", cEffIdx), (title.Length()) ? title : "Interpolators::Efficiency"};
        cEffIdx++;
        if(!multigraph)
            canvas->DivideSquare(fEff.size());
        c = canvas;
    }
    if(!multigraph)
    {
        int idx {1};
        for(const auto& [_, g] : fEff)
        {
            c->cd(idx);
            g->SetLineWidth(2);
            g->Draw("apl");
            idx++;
        }
    }
    else
    {
        auto* mg {new TMultiGraph};
        mg->SetTitle(TString::Format("%s;#theta_{%s} [#circ];#epsilon", ((pad && title.Length()) ? title : "").Data(),
                                     fIsLab ? "Lab" : "CM"));
        for(const auto& [_, eff] : fEff)
        {
            // Get new TGraph to keep this classes independent of TMultiGraph
            auto* g {eff->CreateGraph()};
            g->SetLineWidth(2);
            g->SetFillStyle(0);
            mg->Add(g, "lp");
        }
        mg->Draw("apl plc pmc");
        // Set ranges
        // X
        auto [xmin, xmax] {GetXaxisRange(mg)};
        mg->GetXaxis()->SetLimits(xmin, xmax);
        // Y
        auto ymax {GetYaxisRange(mg)};
        mg->SetMaximum(ymax);
        gPad->Update();
        auto leg {c->BuildLegend()};
        c->GetListOfPrimitives()->RemoveLast();
        leg->SetNColumns(fEff.size() / 5 + 1); // 5 items per legend column
        leg->Draw();
    }
    return c;
}

void Interpolators::Efficiency::SaveAs(const std::string& file)
{
    auto f {std::make_unique<TFile>(file.c_str(), "recreate")};
    for(const auto& [name, eff] : fEff)
        eff->Write(name.c_str());
}

void Interpolators::Sigmas::SaveAs(const std::string& file, const std::string& name)
{
    auto fout {std::make_unique<TFile>(file.c_str(), "recreate")};
    fGraph->Write(name.c_str());
}

void Interpolators::Sigmas::Read(const std::string& file, const std::string& name)
{
    auto fin {std::make_unique<TFile>(file.c_str())};
    fGraph = fin->Get<TGraphErrors>(name.c_str());
    if(!fGraph)
        throw std::runtime_error("Interpolators::Sigmas::Read(): cannot load " + name + " in file " + file);
}

void Interpolators::Sigmas::Draw()
{
    static int cSigmas {};
    auto* c {new TCanvas {TString::Format("cSigmas%d", cSigmas), "Interpolators::Sigmas canvas"}};
    cSigmas++;
    fGraph->SetTitle(";E_{x} [MeV];#sigma_{E} [MeV]");
    fGraph->SetLineWidth(2);
    fGraph->SetLineColor(8);
    fGraph->SetMarkerStyle(24);
    fGraph->Draw("apl");
}
