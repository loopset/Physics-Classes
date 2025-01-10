#ifndef FitInterface_h
#define FitInterface_h

#include <string>
#include <unordered_map>
#include <vector>

namespace Fitters
{
class Interface
{
public:
    using Key = std::string;
    using DoubleVec = std::vector<double>;
    using BoolVec = std::vector<bool>;
    using Pair = std::pair<double, double>;
    using PairVec = std::vector<Pair>;
    typedef std::unordered_map<std::string, DoubleVec> Initial;
    typedef std::unordered_map<std::string, PairVec> Bounds;
    typedef std::unordered_map<std::string, BoolVec> Fixed;
    typedef std::unordered_map<std::string, DoubleVec> Step;

private:
    Initial fInitial {};           //!
    Bounds fBounds {};             //!
    Fixed fFix {};                 //!
    Step fStep {};                 //!
    std::vector<Key> fKeys {};     //!< Keys used for this fit
    std::vector<double> fGuess {}; //!< Ex initial guess, which tags the state
    std::unordered_map<std::string, unsigned int> fNParsType {
        {"g", 3}, {"v", 4}, {"ps", 1}, {"cte", 1}};                       //! Number of pars pers type
    std::vector<std::string> fParNames {"Amp", "Mean", "Sigma", "Gamma"}; //!< Names of parameters
    unsigned int fNGaus {};
    unsigned int fNVoigt {};
    unsigned int fNPS {};
    bool fCte {};

public:
    Interface() = default;

    // Setters
    void AddState(const Key& key, const DoubleVec& init);
    void SetInitial(const Key& key, unsigned int idx, double val);
    void SetBounds(const Key& key, unsigned int idx, const Pair& pair);
    void SetBoundsAll(unsigned int idx, const Pair& pair);
    void DisableBounds(unsigned int idx);
    void SetFix(const Key& key, unsigned int idx, bool fix);
    void SetFixAll(unsigned int idx, bool fix);
    void EndAddingStates();

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


private:
    std::string GetType(const Key& key);
    int GetIndex(const Key& key);
    bool CheckInit(const Key& key, const DoubleVec& init);
    void SetSorting();
    unsigned int GetIdxKey(const Key& key) const;
    void SetDefaults();
};
} // namespace Fitters

#endif // !FitInterface_h
