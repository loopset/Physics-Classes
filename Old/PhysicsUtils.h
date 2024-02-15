#ifndef PhysicsUtils_h
#define PhysicsUtils_h

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpline.h"

#include <string>
namespace PhysicsUtils
{
    class ExperimentInfo
    {
    private:
        double fNt {}; //!< Number of target particles
        double fNtriggers {}; //!< Number of beam triggers
        double fNb {}; //!< Number of beam particles (after applying division factor)
        double fDiv {};//!< Division factor of CATS
        
    public:
        ExperimentInfo(double ntargets, double ntriggers, double div = 1)
            : fNt(ntargets), fNtriggers(ntriggers), fDiv(div)
        {
            fNb = fNtriggers * fDiv * fDiv;
        }
        //Getters
        double GetNt() const {return fNt;}
        double GetNtriggers() const {return fNtriggers;}
        double GetNb() const {return fNb;}
        double GetUNb() const;
    };

    class Uncertainty
    {
    private:
        double fVal;
        double fUVal;
    public:
        Uncertainty(double val, double uval)
            : fVal(val), fUVal(uval)
        {}

        std::string GetStr() const;
    };

    class SigmaInterpolator
    {
    private:
        TGraphErrors* fGraph;
        TGraphErrors* fScaled;//!< Just if we need to scale simulated sigmas against experimental
        TSpline* fSpe;
        double fFactor;
    public:
        SigmaInterpolator(const std::string& file, const std::string& name = "gsigma");

        void ComputeAndSetScalingFactorInGS(double expsimgags);
        void SetScalingFactor(double scaling);
        double EvalSpline(double Ex){return fSpe->Eval(Ex);}
    };
}

#endif
