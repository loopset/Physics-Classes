#ifndef FitModel_h
#define FitModel_h

#include "TH1.h"

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
#include "Math/WrappedMultiTF1.h"

#include <string>
#include <vector>
namespace Fitters
{
class Model : public ROOT::Math::IParametricFunctionMultiDimTempl<double>
{
public:
    typedef std::unordered_map<int, std::vector<double>> ParGroup;
    typedef std::vector<ParGroup> ParPack;

private:
    // PS data
    std::vector<TH1*> fPS {};
    // Inner function settings
    int fNGauss {};
    int fNVoigt {};
    int fNPS {};
    bool fCte {false};
    // Number of pars per func
    int fNParGauss {3};
    int fNParVoigt {4};
    int fNParPS {1}; // just amplitude
    // Configs
    int fNDivBin {20}; // number of divisions in bin
    bool fIsBinned {true};
    // Parameters
    std::vector<double> fPars {};
    std::vector<std::string> fParNames {};

public:
    // Constructor
    Model(int ngaus, int nvoigt, const std::vector<TH1*>& ps = {}, bool withCte = false);

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
    ROOT::Math::IBaseFunctionMultiDimTempl<double>* Clone() const override
    {
        return new Model {fNGauss, fNVoigt, fPS, fCte};
    };
    bool HasGradient() const override { return false; }
    unsigned int NDim() const override { return 1; }

    // Other custom functions
    std::pair<std::string, int> GetTypeIdx(const std::string& name) const;

    // Workaround: wrap this into a TF1
    std::pair<TF1*, ROOT::Math::WrappedMultiTF1> Wrap(double xmin, double xmax);

private:
    void InitParNames();
    ParPack UnpackParameters(const double* pars) const;
    double EvalPS(unsigned int i, double x) const;
    // Override of IBaseFunction
    double DoEvalPar(const double* x, const double* p) const override;
};
} // namespace Fitters

#endif
