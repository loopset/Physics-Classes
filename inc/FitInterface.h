#ifndef FitInterface_h
#define FitInterface_h

#include "TGraphErrors.h"

#include "AngComparator.h"

#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace Fitters
{
class Interface
{
public:
    using Key = std::string;
    using Info = std::string;
    using Pair = std::pair<double, double>;
    using DoubleVec = std::vector<double>;
    using BoolVec = std::vector<bool>;
    using PairVec = std::vector<Pair>;
    typedef std::unordered_map<std::string, DoubleVec> Initial;
    typedef std::unordered_map<std::string, PairVec> Bounds;
    typedef std::unordered_map<std::string, BoolVec> Fixed;
    typedef std::unordered_map<std::string, DoubleVec> Step;

private:
    // Commom settings
    std::vector<Key> fKeys {};                     //!< Keys used for this fit
    std::map<std::string, double> fGuess {};       //!< Ex initial guess, which tags the state
    std::map<std::string, std::string> fLabels {}; //!< Labels for state
    bool fIsRead {};                               //! !< Whether this object is constructed from file or not
    // Settings for Fit stage
    Initial fInitial {}; //!
    Bounds fBounds {};   //!
    Fixed fFix {};       //!
    Step fStep {};       //!
    std::unordered_map<std::string, unsigned int> fNParsType {
        {"g", 3}, {"v", 4}, {"ps", 1}, {"cte", 1}};                       //! !< Number of pars pers type
    std::vector<std::string> fParNames {"Amp", "Mean", "Sigma", "Gamma"}; //! !< Names of parameters
    unsigned int fNGaus {};
    unsigned int fNVoigt {};
    unsigned int fNPS {};
    bool fCte {};
    // Settings for Angular distribution stage
    std::map<std::string, Angular::Comparator> fComparators {}; //!

public:
    Interface() = default;

    // Setters
    // Fit
    void AddState(const Key& key, const DoubleVec& init, const Info& info = {});
    void SetInitial(const Key& key, unsigned int idx, double val);
    void SetBounds(const Key& key, unsigned int idx, const Pair& pair);
    void SetBoundsAll(unsigned int idx, const Pair& pair);
    void DisableBounds(unsigned int idx);
    void SetFix(const Key& key, unsigned int idx, bool fix);
    void SetFixAll(unsigned int idx, bool fix);
    void EndAddingStates();
    // Angular
    void AddAngularDistribution(const Key& key, TGraphErrors* gexp);
    void Do(std::function<void(Angular::Comparator& comp)> func);

    // Getters
    const Initial& GetInitial() const { return fInitial; }
    const Bounds& GetBounds() const { return fBounds; }
    const Fixed& GetFixed() const { return fFix; }
    int GetNGauss() const { return fNGaus; }
    int GetNVoigt() const { return fNVoigt; }
    int GetNPS() const { return fNPS; }
    bool GetCte() const { return fCte; }
    const std::vector<Key>& GetKeys() const { return fKeys; }
    double GetGuess(const Key& key) const;

    // Other methods
    void Print() const;
    void Write(const std::string& file) const;
    void Read(const std::string& file);


private:
    std::string GetType(const Key& key);
    int GetIndex(const Key& key);
    bool CheckInit(const Key& key, const DoubleVec& init);
    void SetSorting();
    unsigned int GetIdxKey(const Key& key) const;
    void SetDefaults();
    void CleanNotStates();
    std::string FormatLabel(const std::string& label);
};
} // namespace Fitters

#endif // !FitInterface_h
