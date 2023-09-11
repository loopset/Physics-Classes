#ifndef AngularDistribution_h
#define AngularDistribution_h

#include "Fit/FitResult.h"
#include "Fitter.h"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TSpline.h"
#include "TVirtualPad.h"
#include "TEfficiency.h"
#include "TSpline.h"
#include "TMultiGraph.h"
#include "ROOT/RDataFrame.hxx"

#include "PhysicsUtils.h"

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace AngularDistribution
{
    class AngularFitter
    {
    private:
        std::vector<TH1F*> fHistos;
        std::vector<Fitters::SpectrumData> fData;
        std::vector<Fitters::SpectrumFunction> fFuncs;
        std::vector<ROOT::Fit::FitResult> fRes;
        Fitters::SpectrumFitter::InitPars fInitPars;
        std::vector<std::unordered_map<std::string, double>> fIntegrals;
    public:
        //Setters
        void ReadFromFile(const std::string& file, const std::string& treename = "InitPars");
        void SetData(double xmin, double xmax, std::vector<TH1F*> hEx);
        void AddPhaseSpace(TH1* hPS);
        void SetFuncModel(int ngaus, int nvoigt, bool cte);
        void SetInitPars(const Fitters::SpectrumFitter::InitPars& initPars){fInitPars = initPars;}

        //Exec calculations
        void Run();
        void ComputeYield(const std::string& mode = "integral", int nsgima = 1);
        //Draw
        TCanvas* Print(double xmin = 0, double xmax = 0);
        //Get results of integral
        const std::vector<std::unordered_map<std::string, double>>& GetAllIntegrals() const {return fIntegrals;};
        std::vector<double> GetIntegralsForPeak(const std::string& key);

    private:
        Fitters::SpectrumFitter::InitBounds GetInitBounds();
        Fitters::SpectrumFitter::FixedPars GetFixedPars();
        void ComputeFuncIntegral(int idx);
        void ComputeCountSum(int idx, int nsigma);
    };

    class Efficiency
    {
    private:
        TEfficiency* fEff;
        TGraphAsymmErrors* fGraph;
        TSpline3* fSpe;
        TSpline3* fUSpe;

    public:
        Efficiency(const std::string& file, const std::string& name);

        double GetPointEff(double thetaCM);
        double GetAveragedEff(double thetaCMmin, double thetaCMmax);
        double GetPointUncertainty(double thetaCM);
        TEfficiency* GetEfficiency() const {return fEff;}
        TGraphAsymmErrors* GetGraph() const {return fGraph;}
        TSpline3* GetSpline() const {return fSpe;}
        TCanvas* GetCanvas(const std::string& opt = "") const;
    private:
        void InitUncertaintySpline();
    };

    class ThetaCMIntervals
    {
    private:
        std::vector<std::pair<double, double>> fVals;
        std::vector<TH1F*> fHistos;
        std::vector<double> fOmega;
        
    public:
        //For differential xs computation
        ThetaCMIntervals(double min, double max, double step, TH1F* hmodel);
        //For absolute xs computation
        ThetaCMIntervals(double min, double max, TH1F* hmodel);
        //Setters and fillers
        void Fill(double thetaCM, double Ex);
        void FillHisto(int idx, double val);
        void FillFromDF(std::vector<ROOT::RDF::RNode> dfs);
        //Getters
        std::pair<double, double> GetInterval(int idx) const {return fVals.at(idx);}
        TH1F* GetHisto(int idx) const {return fHistos.at(idx);}
        double GetOmega(int idx) const {return fOmega.at(idx);}
        std::vector<TH1F*> GetHistos() const {return fHistos;}
        std::vector<std::pair<double, double>> GetIntervals() const {return fVals;}
        std::vector<double> GetOmegas() const {return fOmega;}
        double GetIntervalCentre(int idx);
        int GetSize() const {return fVals.size();}
        TCanvas* Draw();
    private:
        double ComputeAngleSolidElement(double min, double max);
                
    };

    class DiffCrossSection
    {
    private:
        TGraphErrors* fxs;
        
    public:
        DiffCrossSection(const std::string& peak, AngularFitter& ang,
                         ThetaCMIntervals& ivs,
                         Efficiency& eff,
                         PhysicsUtils::ExperimentInfo& exp);

        //Getters
        TGraphErrors* GetExperimental() const {return fxs;}
    private:
        void PerformCalculation(const std::string& peak, AngularFitter& ang,
                                ThetaCMIntervals& ivs,
                                Efficiency& eff,
                                PhysicsUtils::ExperimentInfo& exp);
        double UncertaintyXS(double N, double Nt, double Nb, double Omega, double epsilon,
                             double thetaCMcentre, Efficiency& eff, PhysicsUtils::ExperimentInfo& exp);
    };

    class AbsCrossSection
    {
    private:
        std::pair<double, double> fRange;
        ThetaCMIntervals fIvs;
        AngularFitter fFitter;
        double fxs;
        double fuxs;
    public:
        AbsCrossSection(double thetamin, double thetamax,
                        std::vector<ROOT::RDF::RNode> dfs, TH1F* hmodel);
        void SetRunFitter(double xmin, double xmax, const std::string& reffile, TH1* hPS = nullptr);
        double PerformCalculation(const std::string& peak, Efficiency& eff, PhysicsUtils::ExperimentInfo& exp);
        TCanvas* Draw(double xmin, double xmax);

        //Getters
        double Get() const {return fxs;}
        double GetU() const {return fuxs;}
    private:
        double UncertaintyAbsXS(double N, double epsilon, double uepsilon, PhysicsUtils::ExperimentInfo& exp);
    };
    TMultiGraph* CompareMethods(AngularFitter& ang, const std::string& peak,
                                ThetaCMIntervals& ivs, Efficiency& eff, PhysicsUtils::ExperimentInfo& exp);
}
#endif
