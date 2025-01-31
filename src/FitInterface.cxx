#include "FitInterface.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include "AngComparator.h"
#include "FitUtils.h"
#include "PhysColors.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <ios>
#include <iostream>
#include <iterator>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>


void Fitters::Interface::AddState(const Key& key, const DoubleVec& init, const Info& info)
{
    if(CheckInit(key, init))
    {
        fKeys.push_back(key);
        fPars.insert({key, init});
        fLabels.insert({key, info});
    }
}

std::string Fitters::Interface::GetType(const Key& key)
{
    auto it {key.find_first_of("0123456789")};
    if(it == std::string::npos)
        throw std::invalid_argument("Interface::GetType(): received not numbered parameter " + key);
    auto type {key.substr(0, it)};
    if(!(type == "g" || type == "v" || type == "ps" || type == "cte"))
        throw std::invalid_argument("Interface::GetType(): received wrong key " + key);
    return type;
}

int Fitters::Interface::GetIndex(const Key& key)
{
    auto it {key.find_first_of("0123456789")};
    if(it == std::string::npos)
        throw std::invalid_argument("Interface::GetIndex(): received not numbered parameter " + key);
    return std::stoi(key.substr(it));
}

bool Fitters::Interface::CheckInit(const Key& key, const DoubleVec& init)
{
    auto type {GetType(key)};
    // Check size
    auto size {init.size()};
    if(size != fNParsType[type])
        throw std::runtime_error("Interface::CheckInit(): " + type + "parameter must have only" +
                                 std::to_string(fNParsType[type]) + " initial values");
    // Check if it is already in keys
    auto isRepeated {std::find(fKeys.begin(), fKeys.end(), key) != fKeys.end()};
    if(isRepeated)
        throw std::runtime_error("Interface::CheckInit(): parameter " + key + " is repeated!");
    if(type == "g")
        fNGaus++;
    else if(type == "v")
        fNVoigt++;
    else if(type == "ps")
        fNPS++;
    else if(type == "cte")
        fCte = true;
    return true;
}

void Fitters::Interface::SetSorting()
{
    // This function hardcodes the sorting of the parameters
    std::sort(fKeys.begin(), fKeys.end(),
              [this](const std::string& a, const std::string& b)
              {
                  auto sortType = [this](const std::string& s)
                  {
                      auto type {GetType(s)};
                      if(type == "g")
                          return 0;
                      if(type == "v")
                          return 1;
                      if(type == "ps")
                          return 2;
                      if(type == "cte")
                          return 3;
                      return 4; // default in case unknown type (never happens, guaranteed by CheckInit)
                  };
                  int rankA = sortType(a), rankB = sortType(b);
                  if(rankA != rankB)
                      return rankA < rankB;
                  return GetIndex(a) < GetIndex(b);
              });
    // And also the initial guess
    fGuess.clear();
    for(const auto& key : fKeys)
    {
        auto& initial {fPars[key]};
        if(initial.size() < 2)
            fGuess[key] = initial.front();
        else
            fGuess[key] = initial[1];
    }
}

void Fitters::Interface::SetDefaults()
{
    int idx {};
    for(const auto& key : fKeys)
    {
        // Bounds
        PairVec bounds;
        auto type {GetType(key)};
        for(int i = 0; i < fNParsType[type]; i++)
        {
            // 0 -> Amplitude
            if(i == 0)
                bounds.push_back({0, 100000});
            // 1 -> Mean
            if(i == 1)
            {
                double offset {0.2};
                bounds.push_back({fGuess[key] - offset, fGuess[key] + offset});
            }
            // 2 -> Sigma
            if(i == 2)
                bounds.push_back({0, 1});
            // 3 -> Gamma
            if(i == 3)
                bounds.push_back({0, 10});
        }
        fBounds.insert({key, bounds});

        // Fixed
        BoolVec fixed(fNParsType[type], false); // released by default
        fFix.insert({key, fixed});

        idx++;
    }
}

void Fitters::Interface::EndAddingStates()
{
    SetSorting();
    SetDefaults();
}

void Fitters::Interface::SetInitial(const Key& key, unsigned int idx, double val)
{
    auto& initial {fPars.at(key)};
    if(idx < initial.size())
    {
        initial.at(idx) = val;
        // Update mean guess only if idx == 1
        if(idx == 1)
            fGuess[key] = val;
    }
    // else idx is out of bounds for key. could happen if we set a sigma value for a ps, for instance
}

void Fitters::Interface::SetBounds(const Key& key, unsigned int idx, const Pair& pair)
{
    auto& bounds {fBounds.at(key)};
    if(idx < bounds.size())
        bounds[idx] = pair;
}

void Fitters::Interface::SetOffsetMeanBounds(double offset)
{
    for(auto& [key, vec] : fBounds)
    {
        if(vec.size() > 1)
        {
            auto mean {fPars[key][1]};
            vec[1] = {mean - offset, mean + offset};
        }
    }
}

void Fitters::Interface::SetBoundsAll(unsigned int idx, const Pair& pair)
{
    for(auto& [key, _] : fBounds)
        SetBounds(key, idx, pair);
}

void Fitters::Interface::DisableBounds(unsigned int idx)
{
    for(auto& [key, bounds] : fBounds)
        bounds.at(idx) = {-11, -11};
}

void Fitters::Interface::SetFix(const Key& key, unsigned int idx, bool fix)
{
    auto& fixed {fFix.at(key)};
    if(idx < fixed.size())
        fixed[idx] = fix;
}

void Fitters::Interface::SetFixAll(unsigned int idx, bool fix)
{
    for(auto& [key, _] : fFix)
        SetFix(key, idx, fix);
}

unsigned int Fitters::Interface::GetIdxKey(const Key& key) const
{
    auto it {std::find(fKeys.begin(), fKeys.end(), key)};
    if(it == fKeys.end())
        throw std::runtime_error("Interface::GetIdxKey() : cannot locate idx for " + key);
    auto idx {std::distance(fKeys.begin(), it)};
    return idx;
}

double Fitters::Interface::GetGuess(const Key& key) const
{
    return fGuess.at(key);
}

double Fitters::Interface::GetParameter(const Key& key, unsigned int idx)
{
    if(fPars.count(key))
    {
        auto& vec {fPars[key]};
        if(idx < vec.size())
            return vec[idx];
    }
    return -11; // default value to avoid raising an exception
}

std::vector<double> Fitters::Interface::GetParameterAll(unsigned int idx)
{
    std::vector<double> ret;
    for(const auto& key : fKeys)
        ret.push_back(GetParameter(key, idx));
    return ret;
}

double Fitters::Interface::GetUnc(const Key& key, unsigned int idx)
{
    if(fUncs.count(key))
    {
        auto& vec {fUncs[key]};
        if(idx < vec.size())
            return vec[idx];
    }
    return -11;
}

void Fitters::Interface::Print() const
{
    std::cout << BOLDCYAN << "----- Fitters::Interface -----" << '\n';
    if(!fIsRead)
    {
        for(const auto& key : fKeys)
        {
            const auto& initial {fPars.at(key)};
            const auto& bounds {fBounds.at(key)};
            const auto& fixed {fFix.at(key)};
            std::cout << "-> " << key << " : " << fLabels.at(key) << '\n';
            for(int i = 0; i < initial.size(); i++)
            {
                std::cout << "    " << fParNames[i] << '\n';
                std::cout << "      initial : " << initial[i] << '\n';
                std::cout << "      bounds  : [" << bounds[i].first << ", " << bounds[i].second << "]" << '\n';
                std::cout << "      fixed   ? " << std::boolalpha << fixed[i] << '\n';
            }
        }
    }
    else
    {
        for(const auto& key : fKeys)
        {
            std::cout << "-> " << key << '\n';
            std::cout << "    label : " << fLabels.at(key) << '\n';
            std::cout << "    guess : " << fGuess.at(key) << '\n';
            if(fPars.count(key))
            {
                const auto& pars {fPars.at(key)};
                const auto& uncs {fUncs.at(key)};
                for(int i = 0; i < pars.size(); i++)
                {
                    std::cout << "    " << fParNames[i] << '\n';
                    std::cout << "      fitted : " << pars[i] << " +/- " << uncs[i] << '\n';
                }
            }
        }
    }
    std::cout << RESET;
}

void Fitters::Interface::Write(const std::string& file) const
{
    auto f {std::make_unique<TFile>(file.c_str(), "recreate")};
    f->WriteObject(this, "inter");
}

void Fitters::Interface::CleanNotStates()
{
    auto toClean {[this](const std::string& key)
                  {
                      auto type {GetType(key)};
                      if(type == "ps" || type == "cte")
                          return true;
                      else
                          return false;
                  }};
    // Vector of keys
    for(auto it = fKeys.begin(); it != fKeys.end();)
    {
        if(toClean(*it))
            it = fKeys.erase(it);
        else
            it++;
    }
    // Labels
    for(auto it = fLabels.begin(); it != fLabels.end();)
    {
        if(toClean(it->first))
            it = fLabels.erase(it);
        else
            it++;
    }
    // Guesses
    for(auto it = fGuess.begin(); it != fGuess.end();)
    {
        if(toClean(it->first))
            it = fGuess.erase(it);
        else
            it++;
    }
}

void Fitters::Interface::ReadPreviousFit(const std::string& file)
{
    auto [vals, uncs] {Fitters::ReadInit(file)};
    for(auto& [key, vec] : vals)
    {
        fPars[key] = vec;
        fUncs[key] = uncs[key];
    }
    std::cout << BOLDYELLOW << "Interface::ReadPreviousFit(): read file " << file << RESET << '\n';
}

void Fitters::Interface::EvalSigma(TGraphErrors* gsigma)
{
    for(auto& [key, vec] : fPars)
    {
        if(vec.size() > 2) // sigma is 3rd parameter
            vec[2] = gsigma->Eval(vec[1]);
    }
}

void Fitters::Interface::Read(const std::string& file, const std::string& fitfile)
{
    auto f {std::make_unique<TFile>(file.c_str())};
    auto* inter {f->Get<Fitters::Interface>("inter")};
    if(!inter)
        throw std::runtime_error("Interface::Read(): cannot read written object");
    *this = *inter;
    fIsRead = true;
    delete inter;
    // Clean states that are not gaussian or voigt
    CleanNotStates();
    // Read previous fit if specified
    if(fitfile.length())
        ReadPreviousFit(fitfile);
}

std::string Fitters::Interface::FormatLabel(const std::string& label)
{
    auto it {label.find_first_of("+-")};
    if(it == std::string::npos)
        return label;
    auto j {label.substr(0, it)};
    auto pi {label.substr(it, 1)};
    auto end {label.substr(it + 1)};
    return j + "^{" + pi + "}" + (end.length() == 1 ? "_{" + end + "}" : " " + end);
}

void Fitters::Interface::AddAngularDistribution(const Key& key, TGraphErrors* gexp)
{
    auto str {key + " = " + TString::Format("%.1f MeV ", fGuess[key]).Data() + FormatLabel(fLabels[key])};
    fComparators[key] = Angular::Comparator {str, gexp};
}

void Fitters::Interface::ReadCompConfig(const std::string& file)
{
    // Code partially generated with chatgpt
    // We could have used ActRoot::InputParser but I dont want
    // to introduce library dependences
    std::ifstream streamer(file);
    if(!streamer)
        throw std::runtime_error("Interface::ReadComparator: cannot open config file");

    std::string currentHeader {"global"};
    std::string line;
    while(std::getline(streamer, line))
    {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        // Skip comments and empty lines
        if(line.empty() || line[0] == '#' || line[0] == '%')
            continue;

        // Check for headers
        if(line[0] == '[' && line.back() == ']')
        {
            currentHeader = line.substr(1, line.size() - 2); // Extract header name
            continue;
        }

        // Treat everything else as a key-value pair
        auto sep {line.find(':')};
        if(sep != std::string::npos)
        {
            auto key {line.substr(0, sep)};
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);

            auto value {line.substr(sep + 1)};
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);

            fCompConf[currentHeader].push_back({key, value});
        }
    }
}

std::vector<int> Fitters::Interface::SplitString(const std::string& str)
{
    std::vector<int> ret;
    std::stringstream ss {str};
    std::string token;
    while(std::getline(ss, token, ','))
    {
        token.erase(0, token.find_first_not_of(" \t\n\r"));
        token.erase(token.find_last_not_of(" \t\n\r") + 1);
        ret.push_back(std::stoi(token));
    }
    return ret;
}

template <typename T>
std::optional<T> Fitters::Interface::GetCompOpt(const std::string& opt, const std::string& header)
{
    if(!fCompConf.count(header))
        return std::nullopt;
    auto& vec {fCompConf[header]};
    auto it {std::find_if(vec.begin(), vec.end(),
                          [&](const std::pair<std::string, std::string>& pair) { return pair.first == opt; })};
    if(it != vec.end())
    {
        if constexpr(std::is_arithmetic_v<T>)                   // int, double, bool
            return std::optional<T> {(T)std::stod(it->second)}; // str to double and then to T type
        else if constexpr(std::is_same_v<T, std::string>)       // only string
            return std::optional<std::string>(it->second);
        else if constexpr(!std::is_scalar_v<T>) // vector
            return std::optional<std::vector<int>>(SplitString(it->second));
        else
            throw std::runtime_error("Interface::GetCompOpt(): if else constexpr failed");
    }
    else
        return std::nullopt;
}

void Fitters::Interface::SetCompConfig(const std::string& key, const std::string& value)
{
    if(fCompConf.count("Draw"))
    {
        // Locate key
        auto& vec {fCompConf["Draw"]};
        auto it {std::find_if(vec.begin(), vec.end(),
                              [&](const std::pair<std::string, std::string>& pair) { return pair.first == key; })};
        if(it != vec.end())
            it->second = value;
    }
}

void Fitters::Interface::FillComp()
{
    // Styling options
    auto lc {GetCompOpt<std::vector<int>>("lc")};
    auto ls {GetCompOpt<std::vector<int>>("ls")};
    auto useStyle {lc.has_value()};

    for(auto& [state, comp] : fComparators)
    {
        if(fCompConf.count(state))
        {
            int idx {};
            for(const auto& [key, file] : fCompConf[state])
            {
                // Check whether is file path; if not, it is a setting
                if(file.find("/") != std::string::npos)
                {
                    // if available lc and ls
                    if(useStyle)
                    {
                        int c {1};
                        if(idx < lc.value().size())
                            c = lc.value()[idx];
                        int s {1};
                        if(ls.has_value())
                            if(idx < ls.value().size())
                                s = ls.value()[idx];
                        comp.Add(key, file, c, s);
                    }
                    else
                        comp.Add(key, file);
                    idx++;
                }
            }
        }
    }
}

void Fitters::Interface::FitComp()
{
    // Get general configuration
    auto logy {GetCompOpt<bool>("logy").value()};
    auto withSF {GetCompOpt<bool>("withSF").value()};
    auto offset {GetCompOpt<double>("offset").value()};
    auto save {GetCompOpt<bool>("save").value()};

    // Canvas layout
    // Get number of canvas
    auto size {(int)fKeys.size()};
    int npads {size <= 4 ? size : 4}; // pads per canvas
    auto ncanv {static_cast<int>(std::ceil((double)fKeys.size() / npads))};
    std::vector<TCanvas*> cs;
    for(int c = 0; c < ncanv; c++)
    {
        cs.push_back(new TCanvas {TString::Format("cComp%d", c), TString::Format("Comp canvas %d", c)});
        cs.back()->DivideSquare(npads);
    }

    // Iterate!
    int ic {};  // index of canvas
    int ip {1}; // index of pad
    for(auto& [state, comp] : fComparators)
    {
        comp.Fit();
        if(ip > npads)
        {
            ip = 1;
            ic++;
        }
        // Draw!
        // But first check for local options!
        auto localLogy {GetCompOpt<bool>("logy", state)};
        auto localOffset {GetCompOpt<int>("offset", state)};
        comp.Draw(state, (localLogy ? localLogy.value() : logy), withSF, (localOffset ? localOffset.value() : offset),
                  cs[ic]->cd(ip));
        ip++;
    }
    // And save
    if(save)
    {
        gSystem->mkdir("./Outputs");
        for(int i = 0; i < ncanv; i++)
        {
            cs[i]->cd();
            cs[i]->Modified();
            cs[i]->Update();
            cs[i]->SaveAs(TString::Format("./Outputs/ang_%d.png", i));
        }
        gROOT->SetSelectedPad(nullptr); // to avoid issues later with DrawClone
    }
}

void Fitters::Interface::DoComp()
{
    FillComp();
    FitComp();
}

void Fitters::Interface::WriteComp(const std::string& file)
{
    // Open file
    auto f {std::make_unique<TFile>(file.c_str(), "recreate")};
    // For each comparator, write!
    for(const auto& key : fKeys)
        fComparators[key].Write(key, ""); // file not specified: save in current one
    // Save also list of keys
    f->WriteObject(&fKeys, "Keys");
}

Fitters::Interface::Key Fitters::Interface::GetKeyOfGuess(double guess, double w)
{
    for(int i = 0; i < fKeys.size(); i++)
    {
        auto& k {fKeys[i]};
        auto& g {fGuess[k]};
        if(std::abs(guess - g) < w)
            return k;
    }
    std::cout << BOLDRED << "Interface::GetKeyOfGuess(): cannot locate state for " << guess << " MeV" << RESET << '\n';
    return "";
}

std::string Fitters::Interface::GetTheoCrossSection(const Key& key)
{
    if(!fCompConf.count(key))
        std::cout << BOLDCYAN << "Interface::GetTheoXS(): no state " << key << " in conf file" << RESET << '\n';
    auto simu {GetCompOpt<std::string>("simu", key)};
    if(simu)
    {
        for(const auto& [model, file] : fCompConf.at(key))
            if(model == simu)
                return file;
        std::cout << BOLDRED << "Interface::GetTheoXS(): cannot find model for simulation named " << simu.value()
                  << RESET << '\n';
        return "";
    }
    else
    {
        std::cout << BOLDRED << "Interface::GetTheoXS(): no simu keyword in header " << key << RESET << '\n';
        return "";
    }
}
