#ifndef ModelPlotter_h
#define ModelPlotter_h

#include "TCanvas.h"
#include "TGaxis.h"
#include "TH2.h"
#include "TLatex.h"
#include "TObject.h"
#include "TROOT.h"

#include "PhysSM.h"

#include <string>
#include <utility>
#include <vector>

namespace PlotUtils
{
class ModelToPlot
{
private:
    std::string fName;
    std::vector<double> fEx;
    std::vector<double> fGammas;
    std::vector<std::string> fSF;
    std::vector<std::string> fJp;
    std::vector<int> fColors;

public:
    ModelToPlot() = default;
    ModelToPlot(const std::string& name) : fName(name) {}

    // Setters
    void SetEx(const std::vector<double>& ex)
    {
        fEx = ex;
        SetDefaults();
    }
    void SetGammas(const std::vector<double>& gammas) { Fill(fGammas, gammas); }
    void SetSF(const std::vector<std::string>& sf) { Fill(fSF, sf); }
    void SetJp(const std::vector<std::string>& jp) { Fill(fJp, jp); }
    void SetColors(const std::vector<int>& colors) { Fill(fColors, colors); };
    void SetFromParser(PhysUtils::SMParser* parser);
    void SetUniqueColor(int col);
    void SetColorsFromPalette();

    // Getters
    std::string GetName() const { return fName; }
    std::vector<double>& GetExs() { return fEx; }
    double GetEx(int i) const { return fEx.at(i); }
    const std::vector<double>& GetGammas() const { return fGammas; }
    double GetGamma(int i) const { return fGammas.at(i); }
    const std::vector<std::string>& GetSFs() const { return fSF; }
    std::string GetSF(int i) const { return fSF.at(i); }
    const std::vector<std::string>& GetJps() const { return fJp; }
    std::string GetJp(int i) const { return fJp.at(i); }
    const std::vector<int>& GetColors() const { return fColors; }
    int GetColor(int i) const { return fColors.at(i); }

    // Miscellanea
    bool CheckLabelExists(int i, const std::string& side);

    static std::string FormatSF(double sf, double usf);

private:
    void SetDefaults();
    template <typename T>
    inline void Fill(std::vector<T>& to, const std::vector<T>& from)
    {
        // if(!(from.size() <= to.size()))
        //     return;
        for(int i = 0; i < to.size(); i++)
            to[i] = from[i];
    }
};

class ModelPointers
{
private:
    // Lines
    std::vector<TObject*> fObjs;
    // Right labels on lines (SFs)
    std::vector<TLatex*> fRLabels;
    // Left labels = Jpi
    std::vector<TLatex*> fLLabels;
    // X axis label
    TLatex* fXLabel;

public:
    // Setters
    void AddObj(TObject* line) { fObjs.push_back(line); }
    void AddRightLabel(TLatex* latex) { fRLabels.push_back(latex); }
    void AddLeftLabel(TLatex* latex) { fLLabels.push_back(latex); }
    void SetXLabel(TLatex* label) { fXLabel = label; }
    // Getters
    const std::vector<TObject*> GetObjs() const { return fObjs; }
    TObject* GetObj(int i) const { return fObjs.at(i); }
    std::vector<TLatex*> GetRightLabels() const { return fRLabels; }
    TLatex* GetRightLabel(int i) const { return fRLabels.at(i); }
    std::vector<TLatex*> GetLeftLabels() const { return fLLabels; }
    TLatex* GetLeftLabel(int i) const { return fLLabels.at(i); }
    TLatex* GetXLabel() const { return fXLabel; }
    // Draw
    void Draw();
};

class ModelPlotter
{
private:
    std::pair<double, double> fYRange;
    int fNModels;
    double fWidth {5};
    double fGap {8};
    double fXaxisYpos {-0.75};
    double fLabelOffset {+0.2};
    double fInitialGap {5};
    TH2I* fHist;
    TGaxis* fYAxis;
    // Experimental info
    std::vector<ModelToPlot> fModels;
    // Pointer to plot things
    std::vector<ModelPointers> fPointers;

public:
    ModelPlotter(double ymin, double ymax, int nmodels);

    // Setters
    void AddModel(const ModelToPlot& model);
    void SetYaxisLabel(const std::string& label);

    // Others
    TCanvas* Draw();

private:
    void InitModel(int i);
    void DrawModel(int i);
};

} // namespace PlotUtils

#endif
