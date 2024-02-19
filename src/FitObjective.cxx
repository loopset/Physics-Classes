#include "FitObjective.h"

#include "TMath.h"

#include "Fit/FitUtil.h"

#include "PhysColors.h"

#include <ios>
#include <iostream>

double Fitters::Objective::DoEvalSigma(double nexp, double nfit) const
{
    double sigma {};
    if(nexp == 0)
        sigma = 1.84;
    if(nexp == 1 && nfit <= 1)
        sigma = 0.827;
    if(nexp == 1 && nfit > 1)
        sigma = 2.3;
    if(nexp == 2 && nfit <= 2)
        sigma = 1.292;
    if(nexp == 2 && nfit > 2)
        sigma = 2.64;
    if(nexp > 2 && (nfit > nexp))
        sigma = TMath::Sqrt(nexp) + 1; // upper limit. see eq.(13) in Some remarks on the error analysis in the case of
                                       // poor statistics, by K.h. Schmidt
    if(nexp > 2 && (nfit <= nexp))
        sigma = TMath::Sqrt(nexp); // lower limit
    return sigma;
}

double Fitters::Objective::DoEval(const double* p) const
{
    // Do a chi2 fit
    double res {};
    for(int i = 0, size = fData->GetSize(); i < size; i++)
    {
        auto x {fData->GetX(i)};
        auto yexp {fData->GetY(i)};
        double yfit {};
        if(fUseIntegral)
            yfit = DoEvalWithIntegral(i, p);
        else if(fUseDivisions)
            yfit = DoEvalWithDivisions(x, p);
        else
            yfit = (*fModel)(&x, p);
        // Numerator of Chi2 func
        auto diff {yexp - yfit};
        // Compute sigma
        auto sigma {DoEvalSigma(yexp, yfit)};
        // Chi2 value!
        res += TMath::Power(diff / sigma, 2);
    }
    return res;
}

double Fitters::Objective::DoEvalWithIntegral(int i, const double* p) const
{
    ROOT::Fit::FitUtil::IntegralEvaluator<> ig {*fModel, p};
    double low {fData->GetXLowEdge(i)};
    double up {fData->GetXUpEdge(i)};
    return ig(&low, &up);
}

double Fitters::Objective::DoEvalWithDivisions(double x, const double* p) const
{
    // Unpack parameters (once per call)
    auto pack = fModel->UnpackParameters(p);
    auto gaus = pack[0];
    auto voigt = pack[1];
    auto phase = pack[2];
    auto cte = pack[3];
    // Define steps
    double start {x - 0.5 * fData->GetBinWidth()};
    double step {fData->GetBinWidth() / fNdiv};
    double sum {};
    for(int i = 0; i < fNdiv; i++)
    {
        // Center of division (index + 0.5)
        double xi {start + (i + 0.5) * step};
        sum += fModel->EvalWithPacks(xi, gaus, voigt, phase, cte);
    }
    // Get mean
    sum /= fNdiv;
    return sum;
}

void Fitters::Objective::Print() const
{
    std::cout << BOLDGREEN << "---- Objective func settings ----" << '\n';
    std::cout << "-> UseDivisions ? " << std::boolalpha << fUseDivisions << '\n';
    std::cout << "-> Ndiv         : " << fNdiv << '\n';
    std::cout << "-> UseIntegral ? " << std::boolalpha << fUseIntegral << '\n';
    std::cout << "------------------------------" << RESET << '\n';
}
