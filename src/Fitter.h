#ifndef SpectraFitter_h
#define SpectraFitter_h

#include "TF1.h"
#include "TFitResultPtr.h"
#include "TH1.h"
#include "Fit/Fitter.h"
#include "Fit/FitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TString.h"

#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <map>
#include <vector>

namespace Fitters
{
    class SpectrumData
    {
    private:
        std::vector<double> fX;
        std::vector<double> fY;
        std::vector<std::vector<double>> fPS;
        std::pair<double, double> fRange;
        double fBinWidth;
    public:
        SpectrumData() = default;
        SpectrumData(double min, double max, TH1* h);
        ~SpectrumData() = default;

        //Getters
        const std::vector<double> GetX() const {return fX;}
        const std::vector<double> GetY() const {return fY;}
        const std::vector<double> GetPS(int idx) const {return fPS.at(idx);}
        std::pair<double, double> GetRange() const {return fRange;}
        int GetNBinsX() const {return fX.size();}
        int GetNPS() const {return fPS.size();}
        double GetBinWidth() const {return fBinWidth;}
        std::pair<double, double> GetDataPair(int bin) const {return {fX[bin], fY[bin]};}
        double GetX(int bin) const {return fX[bin];}
        double GetY(int bin) const {return fY[bin];}
        double GetYPS(int idx, int bin) const {return fPS[idx][bin];}
        //Setters
        void SetRange(int min, int max) {fRange = {min, max};}
        void AddPhaseSpace(TH1* h){FillPS(h);}
        int FindBin(double x) const;
        double Integral(double xmin, double xmax) const;
        
    private:
        void SetRangeFromHisto(double min, double max, TH1* h);
        void FillVectors(TH1* h);
        void FillPS(TH1* h);
    };

    class SpectrumFunction
    {
    public:
        typedef std::unordered_map<int, std::vector<double>> ParGroup;
        typedef std::vector<ParGroup> ParPack;
    private:
        const SpectrumData* fData;
        //Functions
        int fNGauss;
        int fNVoigt;
        int fNPS;
        bool fCte {false};
        //Number of pars per func
        int fNParGauss {3};
        int fNParVoigt {4};
        int fNParPS {1};//just amplitude
        //Configs
        int fNDivBin {20};//number of divisions in bin
        //Interal map
        std::vector<std::pair<std::string, int>> fChart;

    public:
        SpectrumFunction() = default;
        SpectrumFunction(int ngaus, int nvoigt, const SpectrumData* data);
        ~SpectrumFunction() = default;

        //Overload operator eval
        double operator() (const double* pars) const;
        //Setters
        void EnableCteBackground() {fCte = true; CreateChart();}
        //Getters
        bool GetEnableCteBackground() const {return fCte;}
        int GetNParFunc(const std::string& key);
        int GetNPars() const;
        int GetNFuncs() const;
        int GetNGaus() const {return fNGauss;}
        int GetNVoigt() const {return fNVoigt;}
        std::vector<std::pair<std::string, int>> GetChart() const {return fChart;}
        int GetDataSize() const {return fData->GetNBinsX();}
        void PrintChart() const;
        inline std::pair<double, double> EvalAfterFit(int bin, const double* pars) const
        {
            return {fData->GetX(bin), EvalInBin(bin, pars)};
        }
        inline double EvalAfterFitByX(double x, const double* pars) const
        {
            return EvalInBin(x, pars, true);
        }
        ParPack UnpackCParameters(const double* pars) const;
    private:
        void CreateChart();
        double EvalInBin(double bin, const double* pars, bool customx = false) const;
        double EvalSigma(double nexp, double nfit) const;
    };

    class SpectrumFitter
    {
    public:
        typedef std::vector<double> DoubleVec;
        typedef std::vector<bool> BoolVec;
        typedef std::vector<std::pair<double, double>> PairVec;
        typedef std::unordered_map<std::string, DoubleVec> InitPars;
        typedef std::unordered_map<std::string, PairVec> InitBounds;
        typedef std::unordered_map<std::string, BoolVec> FixedPars;
        typedef std::unordered_map<std::string, DoubleVec> StepPars;

    private:
        ROOT::Fit::Fitter fFitter;
        ROOT::Fit::FitResult fFitResult;
        SpectrumFunction* fFunc;
        //Fitter setting
        InitPars fInitPars;
        InitBounds fInitBounds;
        FixedPars fFixedPars;
        StepPars fStepPars;

    public:
        SpectrumFitter() = default;
        SpectrumFitter(SpectrumFunction* func, const std::string& minimizer = "Minuit2",
                       const std::string& algorithm = "Migrad", int strategy = 1);
        ~SpectrumFitter() = default;

        bool Fit(bool minos = false);

        //Setters
        void SetInitPars(const InitPars& pars){fInitPars = pars;}
        void SetInitBounds(const InitBounds& bounds){fInitBounds = bounds;}
        void SetFixedPars(const FixedPars& fixed){fFixedPars = fixed;}
        void SetSteps(const StepPars& steps){fStepPars = steps;}
        //Getters
        ROOT::Fit::FitResult GetFitResult() const {return fFitResult;}
        const double* GetFitParameters() const {return fFitResult.GetParams();}
        //Write results to file
        void WriteToFile(const std::string& file);
        void PrintShiftedMeans();
    private:
        void InitCParameters(double* pars);
        std::pair<int, int> LocatePars(const std::string& key);
        template<typename T>
        void AssertDimensions(const std::string& key, const std::vector<T>& vec);
        void ConfigFit();
        void ConfigLabels();
        void ConfigBounds();
        void ConfigFixed();
        void ConfigSteps();
        void PrintParametersAtLimit();
        bool CompareDoubles(double a, double b, double tol = 0.0001);
    };

    class SpectrumPlotter
    {
    private:
        SpectrumData* fData;
        SpectrumFunction* fFunc;
        ROOT::Fit::FitResult fRes;

    public:
        SpectrumPlotter() = default;
        SpectrumPlotter(SpectrumData* data, SpectrumFunction* func, ROOT::Fit::FitResult res);
        ~SpectrumPlotter() = default;

        TGraph* GetGlobalFitGraph();
        TGraphErrors* GetFitResiduals();
        std::unordered_map<std::string, TF1*> GetIndividualFuncs();
        std::unordered_map<std::string, TH1F*> GetIndividualHists();
        TH1F* GetIndividualPS(int idx, const TString& hname = "");

    private:
        void FillHistoFromFunc(TH1F* h, TF1* f);
    };
}

#endif
