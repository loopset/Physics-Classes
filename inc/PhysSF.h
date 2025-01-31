#ifndef PhySF_h
#define PhySF_h

#include "Rtypes.h"

#include "TObject.h"

#include <string>
#include <vector>
namespace PhysUtils
{
//! A class representing a SF and its associated minimization
class SpectroscopicFactor : public TObject
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

    ClassDef(SpectroscopicFactor, 1);
};

//! A class representing a collection of SFs
class SFCollection : public TObject
{
private:
    std::vector<std::string> fModels {};
    std::vector<SpectroscopicFactor> fSFs {};

public:
    SFCollection() = default;

    void Add(const std::string& model, const SpectroscopicFactor& sf);
    const std::vector<std::string>& GetModels() const { return fModels; }
    const std::vector<SpectroscopicFactor>& GetSFs() const { return fSFs; }
    SpectroscopicFactor* Get(const std::string& model);
    SpectroscopicFactor* GetApprox(const std::string& model);
    SpectroscopicFactor* GetBestChi2();
    void Print() const;

    ClassDef(SFCollection, 1);
};
} // namespace PhysUtils
#endif
