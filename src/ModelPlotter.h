#ifndef ModelPlotter_h
#define ModelPlotter_h
#include "Rtypes.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TH2I.h"
#include "TFrame.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"

#include <iostream>
#include <stdexcept>
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
        std::vector<std::string> fSF;
        std::vector<std::string> fJp;
        std::vector<int> fColors;
    public:
        ModelToPlot() = default;
        ModelToPlot(const std::string& name) : fName(name) {}
        //Setters
        void SetEx(const std::vector<double>& ex){fEx = ex;}
        void SetSF(const std::vector<std::string>& sf){fSF = sf;}
        void SetJp(const std::vector<std::string>& jp){fJp = jp;}
        void SetColors(const std::vector<int>& colors){fColors = colors;};
        void SetUniqueColor(int col);
        void SetColorsFromPalette();
        //Getters
        std::string GetName() const {return fName;}
        std::vector<double> GetExs() const {return fEx;}
        double GetEx(int i) const {return fEx.at(i);}
        std::vector<std::string> GetSFs() const {return fSF;}
        std::string GetSF(int i) const {return fSF.at(i);}
        std::vector<std::string> GetJps() const {return fJp;}
        std::string GetJp(int i) const {return fJp.at(i);}
        std::vector<int> GetColors() const {return fColors;}
        int GetColor(int i) const {return fColors.at(i);}
        //Miscellanea
        bool CheckLabelExists(int i, const std::string& side);
    };

    class ModelPointers
    {
    private:
        //Lines
        std::vector<TLine*> fLines;
        //Right labels on lines (SFs)
        std::vector<TLatex*> fRLabels;
        //Left labels = Jpi
        std::vector<TLatex*> fLLabels;
        //X axis label
        TLatex* fXLabel;
    public:
        //Setters
        void AddLine(TLine* line){fLines.push_back(line);}
        void AddRightLabel(TLatex* latex){fRLabels.push_back(latex);}
        void AddLeftLabel(TLatex* latex){fLLabels.push_back(latex);}
        void SetXLabel(TLatex* label){fXLabel = label;}
        //Getters
        std::vector<TLine*> GetLines() const {return fLines;}
        TLine* GetLine(int i) const {return fLines.at(i);}
        std::vector<TLatex*> GetRightLabels() const {return fRLabels;}
        TLatex* GetRightLabel(int i) const {return fRLabels.at(i);}
        std::vector<TLatex*> GetLeftLabels() const {return fLLabels;}
        TLatex* GetLeftLabel(int i) const {return fLLabels.at(i);}
        TLatex* GetXLabel() const {return fXLabel;}
        //Draw
        void Draw();
    };

    class ModelPlotter
    {
    private:
        std::pair<double, double> fYRange;
        int fNModels;
        double fWidth {5};
        double fGap {4.5};
        double fXaxisYpos {-0.5};
        double fLabelOffset {+0.1};
        double fInitialGap {0.75};
        TH2I* fHist;
        TGaxis* fYAxis;
        //Experimental info
        std::vector<ModelToPlot> fModels;
        //Pointer to plot things
        std::vector<ModelPointers> fPointers;
    
    public:
        ModelPlotter(double ymin, double ymax, int nmodels);

        void AddModel(const ModelToPlot& model);
        void SetYaxisLabel(const std::string& label);
        TCanvas* Draw();
    private:
        void InitModel(int i);
        void DrawModel(int i);
    };
    
}

#endif
