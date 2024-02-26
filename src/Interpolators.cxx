#include "Interpolators.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TGraph.h"
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
    // Mark X as being already sorted
    fGraph[peak]->SetBit(TGraph::kIsSortedX);
}

double Interpolators::Efficiency::GetMeanEff(const std::string& peak, double min, double max) const
{
    auto* g {fGraph.at(peak)};
    std::vector<double> vals;
    for(int p = 0; p < g->GetN(); p++)
    {
        if(min <= g->GetPointX(p) && g->GetPointX(p) < max)
            vals.push_back(g->GetPointY(p));
    }
    return TMath::Mean(vals.begin(), vals.end());
}

TCanvas* Interpolators::Efficiency::Draw(bool multigraph)
{
    static int cEffIdx {};
    auto* c {new TCanvas {TString::Format("cEff%d", cEffIdx), "Efficiency canvas"}};
    cEffIdx++;
    if(!multigraph)
    {
        c->DivideSquare(fEff.size());
        int idx {1};
        for(const auto& [_, g] : fEff)
        {
            c->cd(idx);
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
