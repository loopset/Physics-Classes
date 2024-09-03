#include "Interpolators.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TString.h"

#include <memory>
#include <stdexcept>
#include <string>

void Interpolators::Efficiency::Add(const std::string& peak, const std::string& file, const std::string& name)
{
    // Read file
    auto f {std::make_unique<TFile>(file.c_str())};
    // Get TEff obj
    auto eff {f->Get<TEfficiency>(name.c_str())};
    if(!eff)
        throw std::runtime_error("Efficiency::Add(): could not read " + name + " key in file " + file);
    // Push to map
    eff->SetTitle((peak + ";#theta_{CM} [#circ];#epsilon").c_str());
    fEff[peak] = eff;
    // Create graph
    fGraph[peak] = eff->CreateGraph();
    fGraph[peak]->SetName(name.c_str());
    // Mark X as being already sorted
    fGraph[peak]->SetBit(TGraph::kIsSortedX);
}

double Interpolators::Efficiency::GetMeanEff(const std::string& peak, double min, double max) const
{
    auto* g {fGraph.at(peak)};
    std::vector<double> vals;
    for(auto t = min; t < max; t += fMeanStep)
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

TCanvas* Interpolators::Efficiency::Draw(bool multigraph, const TString& title)
{
    static int cEffIdx {};
    auto* c {new TCanvas {TString::Format("cEff%d", cEffIdx), (title.Length()) ? title : "Interpolators::Efficiency"}};
    cEffIdx++;
    if(!multigraph)
    {
        c->DivideSquare(fEff.size());
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
        mg->SetTitle(";#theta_{CM} [#circ];#epsilon");
        for(const auto& [_, eff] : fEff)
        {
            // Get new TGraph to keep this classes independent of TMultiGraph
            auto* g {eff->CreateGraph()};
            g->SetLineWidth(2);
            g->SetFillStyle(0);
            mg->Add(g, "lp");
        }
        mg->Draw("apl plc pmc");
        c->BuildLegend();
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
