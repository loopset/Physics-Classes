#include "TBox.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"

#include "ModelPlotter.h"
#include "PhysSM.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void PlotUtils::ModelToPlot::SetDefaults()
{
    auto size {fEx.size()};
    // Gammas
    fGammas = std::vector<double>(size, 0);
    // SF
    fSF = std::vector<std::string>(size, "");
    // Jpi
    fJp = std::vector<std::string>(size, "");
    // Colors
    fColors = std::vector<int>(size, 1);
}

void PlotUtils::ModelToPlot::SetUniqueColor(int col)
{
    // Resize vector of colors
    fColors.clear();
    fColors.resize(fEx.size());
    for(int i = 0; i < fEx.size(); i++)
        fColors[i] = col;
}

void PlotUtils::ModelToPlot::SetColorsFromPalette()
{
    // Resize vector of colors
    auto nex {fEx.size()};
    fColors.clear();
    fColors.resize(nex);
    // Get number of colors in current palette
    auto ncols {gStyle->GetNumberOfColors()};
    for(int i = 0; i < nex; i++)
    {
        auto color {static_cast<EColor>((float)ncols / nex * i)};
        if(color == 0) // skip white just in case
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

std::string PlotUtils::ModelToPlot::FormatSF(double value, double unc)
{
    // Round uncertainty to 2 significant digits int exponent = static_cast<int>(std::floor(std::log10(unc)));
    int exponent = static_cast<int>(std::floor(std::log10(unc)));
    int decimals = -exponent + 1;
    double rounded_unc = std::round(unc * std::pow(10, decimals)) / std::pow(10, decimals);
    int unc_digits = std::round(rounded_unc * std::pow(10, decimals)); // Extract digits (e.g., 0.0026 → 26)

    // Format value to match decimals
    std::ostringstream oss;
    oss.precision(decimals);
    oss << std::fixed << value;
    std::string str_val = oss.str();

    return str_val + "(" + std::to_string(unc_digits) + ")";
}

void PlotUtils::ModelToPlot::SetFromParser(PhysUtils::SMParser* parser)
{
    std::vector<double> exs;
    std::vector<double> gammas;
    std::vector<std::string> sfs;
    std::vector<std::string> jpis;
    for(const auto& [_, vals] : parser->GetMap())
    {
        for(const auto& val : vals)
        {
            // Ex
            exs.push_back(val.fEx);
            // Gamma
            gammas.push_back(val.fGamma);
            // SF
            sfs.push_back(TString::Format("%.2f", val.fSF).Data());
            // Jpi
            jpis.push_back(val.fQ.Format());
        }
    }
    SetEx(exs);
    SetGammas(gammas);
    SetSF(sfs);
    SetJp(jpis);
}

/////////////////////////////////////////////////////////////////////

void PlotUtils::ModelPointers::Draw()
{
    // Lines
    for(auto& o : fObjs)
        o->Draw();
    // Right labels
    for(auto& label : fRLabels)
        label->Draw();
    // Left labels
    for(auto& label : fLLabels)
        label->Draw();
    // X label
    fXLabel->Draw();
}

//////////////////////////////////////////////////////////////////////

PlotUtils::ModelPlotter::ModelPlotter(double ymin, double ymax, int nmodels) : fYRange({ymin, ymax}), fNModels(nmodels)
{
    // Init histogram
    fHist = new TH2I("hmodel", "", fNModels, 0, (fInitialGap + fWidth) + (fNModels - 1) * (fGap + fWidth), 100,
                     fYRange.first, fYRange.second);
    std::cout << "Xfinal: " << (fInitialGap + fWidth) + (fNModels - 1) * (fGap + fWidth) << '\n';
    // Hide stats
    fHist->SetStats(false);
    // Hide X axis
    fHist->GetXaxis()->SetLabelOffset(999);
    fHist->GetXaxis()->SetAxisColor(0);
    fHist->GetYaxis()->SetLabelOffset(999);
    fHist->GetYaxis()->SetAxisColor(0);
    fHist->GetYaxis()->SetNdivisions(520);
    gStyle->SetGridColor(1); // in case one wants to draw a grid.
    // otherwise the grid takes the color of the axis, and in this case it is white

    // Init new Y axis
    fYAxis = new TGaxis(0, fYRange.first, 0, fYRange.second, fYRange.first, fYRange.second, 510, "-R");
    fYAxis->CenterTitle();
    fYAxis->SetTitleFont(gStyle->GetTitleFont("Y"));
    fYAxis->SetTitleSize(gStyle->GetTitleSize("Y"));
    fYAxis->SetLabelFont(gStyle->GetLabelFont("Y"));
    fYAxis->SetLabelSize(gStyle->GetLabelSize("Y"));
}

void PlotUtils::ModelPlotter::SetYaxisLabel(const std::string& label)
{
    fYAxis->SetTitle(label.c_str());
}

void PlotUtils::ModelPlotter::AddModel(const ModelToPlot& model)
{
    fModels.push_back(model);
    fPointers.push_back({});
    InitModel(fPointers.size() - 1);
}

void PlotUtils::ModelPlotter::InitModel(int i)
{
    // Beginning of line
    double begin {fInitialGap + (fGap + fWidth) * i};
    double end {begin + fWidth};
    double center {begin + fWidth / 2};
    std::cout << "Model " << fModels[i].GetName() << '\n';
    std::cout << "Begin = " << begin << '\n';
    std::cout << "End   = " << end << '\n';
    std::cout << "Center = " << center << '\n';
    std::cout << "===============" << '\n';
    // Label on X axis
    fPointers[i].SetXLabel(new TLatex(center, fXaxisYpos, fModels[i].GetName().c_str()));
    fPointers[i].GetXLabel()->SetTextAlign(23);
    fPointers[i].GetXLabel()->SetTextFont(gStyle->GetTitleFont("X"));
    fPointers[i].GetXLabel()->SetTextSize(gStyle->GetTitleSize("X"));
    // Lines and labels
    int idx {};
    for(const auto& ex : fModels[i].GetExs())
    {
        auto gamma {fModels[i].GetGamma(idx)};
        if(gamma == 0)
        {
            // Line
            auto* line {new TLine(begin, ex, end, ex)};
            line->SetLineWidth(2);
            line->SetLineColor(fModels[i].GetColor(idx));
            fPointers[i].AddObj(line);
        }
        else
        {
            auto w {gamma / 2};
            // Box
            auto* box {new TBox {begin, ex - w, end, ex + w}};
            box->SetLineWidth(2);
            box->SetLineColor(fModels[i].GetColor(idx));
            box->SetFillColor(fModels[i].GetColor(idx));
            fPointers[i].AddObj(box);
        }
        // Right labels
        if(fModels[i].CheckLabelExists(idx, "right"))
        {
            auto label {new TLatex(end + fLabelOffset, ex, fModels[i].GetSF(idx).c_str())};
            label->SetTextFont(gStyle->GetTextFont());
            label->SetTextSize(gStyle->GetTextSize());
            label->SetTextAlign(12);
            fPointers[i].AddRightLabel(label);
        }
        // Left labels
        if(fModels[i].CheckLabelExists(idx, "left"))
        {
            auto label {new TLatex(begin - fLabelOffset, ex, fModels[i].GetJp(idx).c_str())};
            label->SetTextFont(gStyle->GetTextFont());
            label->SetTextSize(gStyle->GetTextSize());
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
    auto* cret {new TCanvas {"cmodel", "Comparison of models"}};
    // Delete frame
    cret->SetFrameLineColor(0);
    // Draw base
    fHist->Draw();
    // Draw axis
    fYAxis->Draw();
    // Draw models!
    for(auto& pointers : fPointers)
        pointers.Draw();
    cret->cd();
    cret->Modified();
    cret->Update();

    return cret;
}
