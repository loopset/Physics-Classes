#ifndef PhysOMP_h
#define PhysOMP_h

#include <string>
namespace PhysOMP
{
class OMP
{
public:
    // All parameters of a general OM potential
    double frc {};   //!< Coulomb radius
    double fVr {};   //!< Real volume depth
    double frv {};   //!< Real volume radius
    double fav {};   //!< Real volume diff
    double fWv {};   //!< Imaginary volume depth
    double frw {};   //!< Imaginary volume radiu
    double faw {};   //!< Imaginary volume diff
    double fWs {};   //!< Imaginary surface depth
    double frs {};   //!< Imaginary surface radius
    double fas {};   //!< Imaginary surface diff
    double fVso {};  //!< Real spin-orbit depth
    double frvso {}; //!< Real spin-orbit radius
    double favso {}; //!< Real spin-orbit diff
    double fWso {};  //!< Imaginary spin-orbit depth
    double frwso {}; //!< Imaginary spin-orbit radius
    double fawso {}; //!< Imaginary spin-orbit diff

    // Particle information
    double fZ {};
    double fN {};
    double fA {};
    double fEnergy {};
    // Potential characteristics
    double fTargetZ {};
    double fTargetA {};
    std::string fName {};

    // OMP
    OMP(int Z, int A, double energy, int targetZ, int targetA, const std::string& name)
        : fZ(Z),
          fA(A),
          fN(A - Z),
          fEnergy(energy),
          fTargetZ(targetZ),
          fTargetA(targetA),
          fName(name)
    {
    }

    void Print() const;
};

class Daehnick : public OMP
{
public:
    Daehnick(int Z, int A, double energy);
};

class Haixia : public OMP
{
public:
    Haixia(int Z, int A, double energy);
};

class Pang : public OMP
{
private:
    bool fIsTriton {}; //!< Bool to set whether is triton or 3He

public:
    Pang(int Z, int A, double energy, bool isTriton = false);
};

class HT1p : public OMP
{
private:
    bool fIsTriton {}; //!< Bool to set whether is triton or 3He

public:
    HT1p(int z, int A, double energy, bool isTrition = false);
};

class KoningDelaroche : public OMP
{
public:
    KoningDelaroche(int Z, int A, double energy);
};
} // namespace PhysOMP

#endif
