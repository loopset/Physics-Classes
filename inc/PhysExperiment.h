#ifndef PhysExperiment_h
#define PhysExperiment_h

#include "TMath.h"
namespace PhysUtils
{
class Experiment
{
private:
    double fNt {};
    double fNtriggers {};
    double fNb {};
    double fNdiv {};

public:
    Experiment(double ntargets, double ntriggers, double div, int ndiv = 1) : fNt(ntargets), fNtriggers(ntriggers)
    {
        fNdiv = TMath::Power(div, ndiv);
        fNb = fNtriggers * fNdiv;
    }
    // Getters
    double GetNt() const { return fNt; }
    double GetUNt() const { return 0; }
    double GetNb() const { return fNb; }
    double GetUNb() const { return TMath::Sqrt(fNtriggers) * fNdiv; }
    double GetDivFactor() const { return fNdiv; }
};
} // namespace PhysUtils
#endif // !PhysExperiment_h
