#ifndef PhysSM_h
#define PhysSM_h

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace PhysUtils
{

class QuantumNumbers
{
public:
    int fN {};
    int fL {};
    double fJ {};
    // Default ctor
public:
    QuantumNumbers() = default;
    QuantumNumbers(int n, int l, double j) : fN(n), fL(l), fJ(j) {}
    // Allows usage as map key
    friend bool operator<(const QuantumNumbers& a, const QuantumNumbers& b)
    {
        if(a.fN != b.fN)
            return a.fN < b.fN;
        if(a.fL != b.fL)
            return a.fL < b.fL;
        return a.fJ < b.fJ;
    }
    inline void Print() const
    {
        std::cout << "-- QuantumNumers --" << '\n';
        std::cout << " (n, l, j) : (" << fN << ", " << fL << ", " << fJ << ")" << '\n';
    }
    std::string Format() const;
};

class SMData
{
public:
    QuantumNumbers fQ {};
    double fEx {};
    double fSF {};
    double fuSF {};
    double fGamma {};

public:
    SMData() = default;
    SMData(double ex, double sf, double gamma = 0) : fEx(ex), fSF(sf), fGamma(gamma) {}
    void SetuSF(double usf) { fuSF = usf; }
    void AddQuantumNumbers(const QuantumNumbers& q) { fQ = q; }
    friend bool operator<(const SMData& a, const SMData& b)
    {
        if(a.fEx == b.fEx)
            return a.fSF > b.fSF;
        return a.fEx < b.fEx;
    }
    inline void Print() const
    {
        std::cout << "---- ModelData ----" << '\n';
        std::cout << " (Ex, SF) : (" << fEx << ", " << fSF << " +- " << fuSF << ")" << '\n';
    }
};

class SMParser
{
public:
    using QuantumMap = std::map<QuantumNumbers, std::set<SMData>>;

private:
    QuantumMap fMap {};
    double fSFgs {}; //!< Ground-state SF
    double fExgs {}; //!< Ground-state Ex

public:
    SMParser() = default;
    SMParser(const std::vector<std::string>& files);

    // Setters
    void AddFile(const std::string& file) { ParseFile(file); }
    // Getters
    QuantumMap& GetMap() { return fMap; }
    std::set<SMData>& GetDataFor(const QuantumNumbers& q) { return fMap.at(q); }
    double GetGroundStateEx() const { return fExgs; }
    double GetGroundStateSF() const { return fSFgs; }

    // Actions
    void ShiftEx();
    void SFRelativeToGS();
    void MaskSFBelow(double thresh);
    void MaskExAbove(double thres);


private:
    void ParseFile(const std::string& file);
    std::string StripSpaces(std::string str);
};
} // namespace PhysUtils
#endif // !PhysSM_h
