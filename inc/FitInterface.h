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
    Initial fInitial {};                  //!
    Bounds fBounds {};                    //!
    Fixed fFix {};                        //!
    Step fStep {};                        //!
    std::vector<Key> fKeys {};            //!< Keys used for this fit
    std::vector<double> fInitialGuess {}; //!< Ex initial guess, which tags the state
    std::unordered_map<std::string, unsigned int> fNParsType {
        {"g", 3}, {"v", 4}, {"ps", 1}, {"cte", 1}}; //! Number of pars pers type
    unsigned int fNGaus {};
    unsigned int fNVoigt {};
    unsigned int fNPS {};
    bool fCte {};

public:
    Interface() = default;

    // Setters
    void Init(const Key& key, const DoubleVec& init);
    void SetBounds(const Key& key, unsigned int idx, const Pair& pair);
    void SetFix(const Key&, unsigned int idx, bool fix);

    // Getters
    const Initial& GetInitial() const { return fInitial; }

private:
    std::string GetType(const Key& key);
    bool CheckInit(const Key& key, const DoubleVec& init);
    void SetDefaults();
};
} // namespace Fitters

#endif // !FitInterface_h
