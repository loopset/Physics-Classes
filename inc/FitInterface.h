#ifndef FitInterface_h
#define FitInterface_h

#include "AngComparator.h"

#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

// forward declaration
class TGraphErrors;

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
    std::vector<Key> fKeys {};       //!< Keys used for this fit
    std::map<Key, double> fGuess {}; //!< Ex initial guess, which tags the state
    std::map<Key, Info> fLabels {};  //!< Labels for state
    bool fIsRead {};                 //! !< Whether this object is constructed from file or not
    // Settings for Fit stage
    Initial fPars {};  //! //!< Initial or (once fitted) parameter values
    Initial fUncs {};  //! //!< Uncertainties of parameters
    Bounds fBounds {}; //!
    Fixed fFix {};     //!
    Step fStep {};     //!
    std::unordered_map<std::string, unsigned int> fNParsType {
        {"g", 3}, {"v", 4}, {"ps", 1}, {"cte", 1}};                       //! !< Number of pars pers type
    std::vector<std::string> fParNames {"Amp", "Mean", "Sigma", "Gamma"}; //! !< Names of parameters
    unsigned int fNGaus {};
    unsigned int fNVoigt {};
    unsigned int fNPS {};
    bool fCte {};
    // Settings for Angular distribution stage
    std::map<std::string, Angular::Comparator> fComparators {};                           //!
    std::map<std::string, std::vector<std::pair<std::string, std::string>>> fCompConf {}; //!

public:
    Interface() = default;

    // Setters
    // Fit
    void AddState(const Key& key, const DoubleVec& init, const Info& info = {});
    void SetInitial(const Key& key, unsigned int idx, double val);
    void SetBounds(const Key& key, unsigned int idx, const Pair& pair);
    void SetOffsetMeanBounds(double offset);
    void SetBoundsAll(unsigned int idx, const Pair& pair);
    void DisableBounds(unsigned int idx);
    void SetFix(const Key& key, unsigned int idx, bool fix);
    void SetFixAll(unsigned int idx, bool fix);
    void EndAddingStates();
    void ReadPreviousFit(const std::string& file);
    void EvalSigma(TGraphErrors* gsigma);
    // Angular
    void AddAngularDistribution(const Key& key, TGraphErrors* gexp);
    void ReadCompConfig(const std::string& file);
    void SetCompConfig(const std::string& key, const std::string& value);
    void FillComp();
    void FitComp();
    void DoComp();
    void WriteComp(const std::string& file);

    // Getters
    const Initial& GetInitial() const { return fPars; }
    const Bounds& GetBounds() const { return fBounds; }
    const Fixed& GetFixed() const { return fFix; }
    int GetNGauss() const { return fNGaus; }
    int GetNVoigt() const { return fNVoigt; }
    int GetNPS() const { return fNPS; }
    bool GetCte() const { return fCte; }
    const std::vector<Key>& GetKeys() const { return fKeys; }
    const std::vector<Key>& GetPeaks() const { return GetKeys(); }
    double GetParameter(const Key& key, unsigned int idx);
    std::vector<double> GetParameterAll(unsigned int idx);
    double GetUnc(const Key& key, unsigned int idx);
    std::vector<double> GetUncAll(unsigned int idx);
    double GetGuess(const Key& key) const;
    Angular::Comparator* GetComp(const Key& key) { return &fComparators.at(key); }
    Key GetKeyOfGuess(double guess, double width = 0.15);
    std::string GetTheoCrossSection(const Key& key);

    // Other methods
    void Print() const;
    void Write(const std::string& file) const;
    void Read(const std::string& file, const std::string& fitfile = "");


private:
    std::string GetType(const Key& key);
    int GetIndex(const Key& key);
    bool CheckInit(const Key& key, const DoubleVec& init);
    void SetSorting();
    unsigned int GetIdxKey(const Key& key) const;
    void SetDefaults();
    void CleanNotStates();
    std::string FormatLabel(const std::string& label);
    template <typename T>
    std::optional<T> GetCompOpt(const std::string& opt, const std::string& header = "Draw");
    std::vector<int> SplitString(const std::string& str);
};
} // namespace Fitters

#endif // !FitInterface_h
