#ifndef FitModel_h
#define FitModel_h

#include "TH1.h"
#include "TSpline.h"

#include "Math/IParamFunction.h"

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
    // Configuration options
    bool fUseSpline {false};

public:
    // Constructor
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

private:
    void InitSplines();
    void InitFuncParNames();
    void InitParNames();
    std::pair<std::string, int> GetTypeIdx(const std::string& name) const;
    // Override of IBaseFunction
    double DoEvalPar(const double* x, const double* p) const override;
};
} // namespace Fitters

#endif
