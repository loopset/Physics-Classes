#ifndef ModelPlotter_cxx
#define ModelPlotter_cxx
#include "ModelPlotter.h"
#include "TLatex.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void PlotUtils::ModelToPlot::SetUniqueColor(int col)
{
    //Resize vector of colors
    fColors.clear();
    fColors.resize(fEx.size());
    for(int i = 0; i < fEx.size(); i++)
        fColors[i] = col;
}

void PlotUtils::ModelToPlot::SetColorsFromPalette()
{
    //Resize vector of colors
    auto nex {fEx.size()};
    fColors.clear();
    fColors.resize(nex);
    //Get number of colors in current palette
    auto ncols {gStyle->GetNumberOfColors()};
    for(int i = 0; i < nex; i++)
    {
        auto color {static_cast<EColor>((float)ncols / nex * i)};
        if(color == 0)//skip white just in case
            color = kAzure;
        fColors[i] = color;
    }
}

bool PlotUtils::ModelToPlot::CheckLabelExists(int i, const std::string& side)
{
    auto lambda = [](const std::vector<std::string>& vec, const int index) -> bool
    {
        bool isOk {false};
        if(index < vec.size())
            isOk = true;
        return isOk;
    };
    if(side == "left")
        return lambda(fJp, i);
    else if(side == "right")
        return lambda(fSF, i);
    else
        throw std::runtime_error("Checking for existance of label: only left or side");
}

/////////////////////////////////////////////////////////////////////

void PlotUtils::ModelPointers::Draw()
{
    //Lines
    for(auto& line : fLines)
        line->Draw();
    //Right labels
    for(auto& label : fRLabels)
        label->Draw();
    //Left labels
    for(auto& label : fLLabels)
        label->Draw();
    //X label
    fXLabel->Draw();
}

//////////////////////////////////////////////////////////////////////

PlotUtils::ModelPlotter::ModelPlotter(double ymin, double ymax, int nmodels)
    : fYRange({ymin, ymax}), fNModels(nmodels)
{
    //Init histogram
    fHist = new TH2I("hmodel", "", fNModels, 0,
                     (fInitialGap + fWidth) + (fNModels - 1) * (fGap + fWidth),
                     100, fYRange.first, fYRange.second);
    //Hide stats
    fHist->SetStats(false);
    //Hide X axis
    fHist->GetXaxis()->SetLabelOffset(999);
    fHist->GetXaxis()->SetAxisColor(0);
    fHist->GetYaxis()->SetLabelOffset(999);
    fHist->GetYaxis()->SetAxisColor(0);

    //Init new Y axis
    fYAxis = new TGaxis(0, fYRange.first, 0, fYRange.second,
                        fYRange.first, fYRange.second, 510, "-R");
    fYAxis->CenterTitle();
}

void PlotUtils::ModelPlotter::SetYaxisLabel(const std::string &label)
{
    fYAxis->SetTitle(label.c_str());
}

void PlotUtils::ModelPlotter::AddModel(const ModelToPlot &model)
{
    fModels.push_back(model);
    fPointers.push_back({});
    InitModel(fPointers.size() - 1);
}

void PlotUtils::ModelPlotter::InitModel(int i)
{
    //Beginning of line
    double begin {fInitialGap + (fGap + fWidth) * i};
    double end   {begin + fWidth};
    double center {begin + fWidth / 2};
    std::cout<<"Model "<<fModels[i].GetName()<<'\n';
    std::cout<<"Begin = "<<begin<<'\n';
    std::cout<<"End   = "<<end<<'\n';
    std::cout<<"Center = "<<center<<'\n';
    std::cout<<"==============="<<'\n';
    //Label on X axis
    fPointers[i].SetXLabel(new TLatex(center, fXaxisYpos, fModels[i].GetName().c_str()));
    fPointers[i].GetXLabel()->SetTextAlign(21);
    //Lines and labels
    int idx {};
    for(const auto& ex : fModels[i].GetExs())
    {
        //lines
        auto line {new TLine(begin, ex, end, ex)};
        line->SetLineWidth(2);
        line->SetLineColor(fModels[i].GetColor(idx));
        fPointers[i].AddLine(line);
        //Right labels
        if(fModels[i].CheckLabelExists(idx, "right"))
        {
            auto label {new TLatex(end + fLabelOffset, ex, fModels[i].GetSF(idx).c_str())};
            label->SetTextSize(0.035);
            label->SetTextAlign(12);
            fPointers[i].AddRightLabel(label);
        }
        //Left labels
        if(fModels[i].CheckLabelExists(idx, "left"))
        {
            auto label {new TLatex(begin - fLabelOffset, ex, fModels[i].GetJp(idx).c_str())};
            label->SetTextSize(0.035);
            label->SetTextAlign(32);
            fPointers[i].AddLeftLabel(label);
        }
        idx++;
    }
}

void PlotUtils::ModelPlotter::DrawModel(int i)
{
    fPointers[i].Draw();
}

TCanvas* PlotUtils::ModelPlotter::Draw()
{
    auto* cret {new TCanvas("cmodel")};
    //Delete frame
    cret->SetFrameLineColor(0);
    //Draw base
    fHist->Draw();
    //Draw axis
    fYAxis->Draw();
    //Draw models!
    for(auto& pointers : fPointers)
        pointers.Draw();
    
    return cret;
}
#endif
