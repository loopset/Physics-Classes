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
        std::string fName {};
        std::vector<std::string> fKeys;
        std::unordered_map<std::string, TGraphErrors*> fTheo;
        std::unordered_map<std::string, TGraph*> fAsym;
        std::unordered_map<std::string, TGraphErrors*> fFits;
        //Save colors and markers in plc and pmc to share between theo and fits
        std::unordered_map<std::string, std::pair<Color_t, Marker_t>> fStyle;
        //Store fit parameters
        std::unordered_map<std::string, std::pair<double, double>> fSF;
        
    public:
        TwoFNR() = default;
        TwoFNR(const std::string& name) : fName(name) {}

        //Setters
        void Add(const std::string& name, const std::string& file);
        //Getters
        TGraphErrors* GetTheoretical(const std::string& name) const {return fTheo.at(name);}
        TGraph* GetAsymmetry(const std::string& name) const {return fAsym.at(name);}
        TGraphErrors* GetFitted(const std::string& name) const {return fFits.at(name);}
        void FitToExperimental(TGraphErrors* gexp, double xmin = 0, double xmax = 0);
        double Integral(const std::string& key, double thetamin, double thetamax);
        void IntegralAll(double thetamin, double thetamax, double exp = -1, double uexp = -1);
        //Draw
        void DrawTheoretical(const std::string& opts = "");
        void DrawFitted();
        TLegend* DrawLegend(bool fancy = false);
        TCanvas* GetCanvas(TGraphErrors* gexp, double xmin = -1, double xmax = -1, bool asym = false);
        TCanvas* GetCanvasPublication(TGraphErrors* gexp, double xmin = -1, double xmax = -1, const std::vector<int>& ls = {}, const std::vector<int>& colors = {}, bool sf = false);
        
    private:
        TGraphErrors* ReadDiffXS(const std::string& file);
        TGraph* ReadAsym(const std::string& file);
    };
}

#endif
