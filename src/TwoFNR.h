#ifndef TwoFNR_h
#define TwoFNR_h

#include "RtypesCore.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
namespace TheoreticalUtils
{
    class TwoFNR
    {
    private:
        std::vector<std::string> fKeys;
        std::unordered_map<std::string, TGraphErrors*> fTheo;
        std::unordered_map<std::string, TGraph*> fAsym;
        std::unordered_map<std::string, TGraphErrors*> fFits;
        //Save colors and markers in plc and pmc to share between theo and fits
        std::unordered_map<std::string, std::pair<Color_t, Marker_t>> fStyle;
        //Store fit parameters
        std::unordered_map<std::string, std::pair<double, double>> fPar;
        
    public:
        TwoFNR() = default;
        TwoFNR(const std::string& name, const std::string& file);

        //Setters
        void Add(const std::string& name, const std::string& file);
        //Getters
        TGraphErrors* GetTheoretical(const std::string& name) const {return fTheo.at(name);}
        TGraph* GetAsymmetry(const std::string& name) const {return fAsym.at(name);}
        TGraphErrors* GetFitted(const std::string& name) const {return fFits.at(name);}
        void FitToExperimental(TGraphErrors* gexp, double xmin = 0, double xmax = 0);
        //Draw
        void DrawTheoretical();
        void DrawFitted();
        TLegend* DrawLegend(bool fancy = false);
        TCanvas* GetCanvas(TGraphErrors* gexp, bool asym = false);
        TCanvas* GetCanvasPublication(TGraphErrors* gexp, double xmin = -1, double xmax = -1);
        
    private:
        TGraphErrors* ReadDiffXS(const std::string& file);
        TGraph* ReadAsym(const std::string& file);
    };
}

#endif