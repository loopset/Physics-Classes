#include "Interpolators.h"

#include "TEfficiency.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"

#include <memory>
#include <stdexcept>
#include <string>

Interpolators::Efficiency::Efficiency(const std::string& file, const std::string& key)
{
    // Read file
    auto f {std::make_unique<TFile>(file.c_str())};
    // Get TEff obj
    fEff = f->Get<TEfficiency>(key.c_str());
    if(!fEff)
        throw std::runtime_error("Efficiency(): could not read " + key + " key in file");
    // Create graph
    fGraph = fEff->CreateGraph();
    // Mark X as being already sorted
    fGraph->SetBit(TGraph::kIsSortedX);
}

double Interpolators::Efficiency::GetMeanEff(double min, double max) const
{
    std::vector<double> vals;
    for(int p = 0; p < fGraph->GetN(); p++)
    {
        if(min <= fGraph->GetPointX(p) && fGraph->GetPointX(p) < max)
            vals.push_back(fGraph->GetPointY(p));
    }
    return TMath::Mean(vals.begin(), vals.end());
}
