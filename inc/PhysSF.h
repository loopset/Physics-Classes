#ifndef PhySF_h
#define PhySF_h

namespace PhysUtils
{
class SpectroscopicFactor
{
private:
    double fSF {};      //!< SF from fit to experimental angular distribution
    double fUSF {};     //!< uncertainty from fit
    double fChi2Red {}; //!< Chi2 reduced from fit
    int fNdf {};        //!< N degrees of freedom in fit

public:
    SpectroscopicFactor() = default;
    SpectroscopicFactor(double sf, double usf, double chi2red, int dof)
        : fSF(sf),
          fUSF(usf),
          fChi2Red(chi2red),
          fNdf(dof)
    {
    }

    double GetSF() const { return fSF; }
    double GetUSF() const { return fUSF; }
    double GetChi2Red() const { return fChi2Red; }
    int GetNDF() const { return fNdf; }

    void Print() const;
};
} // namespace PhysUtils
#endif
