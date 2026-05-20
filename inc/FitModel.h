#ifndef FitModel_h
#define FitModel_h

#include "TF1.h"
#include "TF1Convolution.h"
#include "TH1.h"
#include "TSpline.h"

#include "Math/IParamFunction.h"

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Fitters
{
class Model : public ROOT::Math::IParametricFunctionMultiDimTempl<double>
{
public:
    typedef std::vector<std::vector<double>> ParPack;
    typedef std::vector<ParPack> ParVec;
    using GammaFunc = std::function<double(double, double, double)>;

private:
    // PS data (by copied histograms)
    std::vector<TH1D> fPS {};
    // Splines of PS to plot global fit
    std::vector<std::shared_ptr<TSpline3>> fSpePS {};
    // Number of inner functions to use: gauss, voigt, ps and cte
    int fNGauss {};
    int fNVoigt {};
    int fNPS {};
    bool fCte {false};
    // Number of pars per func
    int fNParGauss {3};
    int fNParVoigt {4};
    int fNParPS {1}; // just amplitude
    // Names of parameters per function
    std::unordered_map<std::string, std::vector<std::string>> fFuncParNames {};
    // Parameters
    std::vector<double> fPars {};
    std::vector<std::string> fParNames {};
    std::vector<std::pair<std::string, unsigned int>> fChart;
    // Store penetrability functions for each voigt
    std::map<int, GammaFunc> fGammaFuncs {};
    std::map<int, std::shared_ptr<TF1>> fGaussians {};
    std::map<int, std::shared_ptr<TF1>> fBWs {};
    std::map<int, std::shared_ptr<TF1Convolution>> fConvObjs;
    // Number of points to compute convolution (integral) for each value of Ex
    int fNConvolutionPoints {};
    // Configuration options
    bool fUseSpline {false};

public:
    // Constructor
    Model() = default;
    Model(int ngaus, int nvoigt, const std::vector<TH1D>& ps = {}, bool withCte = false);

    // Derived functions from IBasePars
    unsigned int NPar() const override;
    std::string ParameterName(unsigned int i) const override { return fParNames[i]; }
    const double* Parameters() const override { return &fPars.front(); }
    void SetParameters(const double* p) override
    {
        for(unsigned int i = 0, size = NPar(); i < size; i++)
        {
            fPars[i] = p[i];
        }
    }

    // Derived functions from IBaseFunction
    Model* Clone() const override
    {
        auto ret {new Model {fNGauss, fNVoigt, fPS, fCte}};
        ret->SetParameters(Parameters());
        ret->SetUseSpline(fUseSpline);
        ret->SetGammaFuncs(fGammaFuncs);
        ret->SetBWs(fBWs);
        ret->SetGaussians(fGaussians);
        ret->SetConvObjs(fConvObjs);
        return ret;
    };
    bool HasGradient() const override { return false; }
    unsigned int NDim() const override { return 1; }

    // Custom getters and setters
    void SetUseSpline(bool use) { fUseSpline = use; }
    bool GetUseSpline() const { return fUseSpline; }
    double EvalPS(unsigned int i, double x) const;
    double EvalWithPacks(double x, ParPack& gaus, ParPack& voigt, ParPack& phase, ParPack& cte) const;

    // Other custom functions to get type func and idx from par name and viceversa
    unsigned int GetIdxFromLabel(const std::string& typeIdx, unsigned int par) const;
    ParVec UnpackParameters(const double* pars) const;

    void Print() const;

    // Add option to have penetrabilities
    void AddBWL(int vIdx, int l, double s, double mu, double R);
    void SetGammaFuncs(const std::map<int, GammaFunc>& funcs) { fGammaFuncs = funcs; }
    const std::map<int, GammaFunc>& GetGammaFuncs() const { return fGammaFuncs; }
    void SetBWs(const std::map<int, std::shared_ptr<TF1>>& bws) { fBWs = bws; }
    const std::map<int, std::shared_ptr<TF1>>& GetBWs() const { return fBWs; }
    void SetGaussians(const std::map<int, std::shared_ptr<TF1>>& gauss) { fGaussians = gauss; }
    const std::map<int, std::shared_ptr<TF1>>& GetGaussians() const { return fGaussians; }
    void SetConvObjs(const std::map<int, std::shared_ptr<TF1Convolution>>& convObjs) { fConvObjs = convObjs; }
    const std::map<int, std::shared_ptr<TF1Convolution>>& GetConvObjs() const { return fConvObjs; }
    // Function to evaluate manual convolution with Gaussian
    void TriggerConvolution(const double* p, double xMin, double xMax);

private:
    void InitSplines();
    void InitFuncParNames();
    void InitParNames();
    std::pair<std::string, int> GetTypeIdx(const std::string& name) const;
    // Override of IBaseFunction
    double DoEvalPar(const double* x, const double* p) const override;
    // Function to initialize penetrability function for given l, s, mu and R
    GammaFunc InitLambda(int l, double s, double mu, double R);
    // Helper methods
    void InitOrUpdateBW(int vIdx, double mean, double Gamma0, double xMin, double xMax);
    void InitOrUpdateGaussian(int vIdx, double sigma, double xMin, double xMax);
    void InitOrUpdateConvolution(int vIdx, double mean, double sigma, double Gamma0, double xMin, double xMax);
};
} // namespace Fitters

#endif
