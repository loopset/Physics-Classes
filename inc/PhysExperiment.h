#ifndef PhysExperiment_h
#define PhysExperiment_h

#include "TMath.h"

#include <string>
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
    Experiment(const std::string& file) { Read(file); }
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
    void Write(const std::string& file) const;
    void Read(const std::string& file);
    void Print() const;
};
} // namespace PhysUtils
#endif // !PhysExperiment_h
