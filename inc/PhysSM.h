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

class ModelData
{
public:
    QuantumNumbers fQ {};
    double fEx {};
    double fSF {};
    double fGamma {};

public:
    ModelData() = default;
    ModelData(double ex, double sf, double gamma = 0) : fEx(ex), fSF(sf), fGamma(gamma) {}
    void AddQuantumNumbers(const QuantumNumbers& q) { fQ = q; }
    friend bool operator<(const ModelData& a, const ModelData& b)
    {
        if(a.fEx == b.fEx)
            return a.fSF > b.fSF;
        return a.fEx < b.fEx;
    }
    inline void Print() const
    {
        std::cout << "---- ModelData ----" << '\n';
        std::cout << " (Ex, SF) : (" << fEx << ", " << fSF << ")" << '\n';
    }
};

class ModelParser
{
public:
    using QuantumMap = std::map<QuantumNumbers, std::set<ModelData>>;

private:
    QuantumMap fMap {};

public:
    ModelParser() = default;
    ModelParser(const std::vector<std::string>& files);

    // Setters
    void AddFile(const std::string& file) { ParseFile(file); }
    // Getters
    QuantumMap& GetMap() { return fMap; }
    std::set<ModelData>& GetDataFor(const QuantumNumbers& q) { return fMap.at(q); }

    // Actions
    void ShiftEx();
    void MaskSFBelow(double thresh);
    void MaskExAbove(double thres);


private:
    void ParseFile(const std::string& file);
    std::string StripSpaces(std::string str);
};
} // namespace PhysUtils
#endif // !PhysSM_h
