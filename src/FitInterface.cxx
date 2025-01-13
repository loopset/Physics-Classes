#include "FitInterface.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include "AngComparator.h"
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
#include <vector>


void Fitters::Interface::AddState(const Key& key, const DoubleVec& init, const Info& info)
{
    if(CheckInit(key, init))
    {
        fKeys.push_back(key);
        fInitial.insert({key, init});
        fLabels.insert({key, info});
    }
}

std::string Fitters::Interface::GetType(const Key& key)
{
    auto it {key.find_first_of("0123456789")};
    if(it == std::string::npos)
        throw std::invalid_argument("Fitters::Interface::GetType(): received not numbered parameter " + key);
    auto type {key.substr(0, it)};
    if(!(type == "g" || type == "v" || type == "ps" || type == "cte"))
        throw std::invalid_argument("Fitters::Interface::GetType(): received wrong key " + key);
    return type;
}

int Fitters::Interface::GetIndex(const Key& key)
{
    auto it {key.find_first_of("0123456789")};
    if(it == std::string::npos)
        throw std::invalid_argument("Fitters::Interface::GetIndex(): received not numbered parameter " + key);
    return std::stoi(key.substr(it));
}

bool Fitters::Interface::CheckInit(const Key& key, const DoubleVec& init)
{
    auto type {GetType(key)};
    // Check size
    auto size {init.size()};
    if(size != fNParsType[type])
        throw std::runtime_error("Fitters::Interface::CheckInit(): " + type + "parameter must have only" +
                                 std::to_string(fNParsType[type]) + " initial values");
    // Check if it is already in keys
    auto isRepeated {std::find(fKeys.begin(), fKeys.end(), key) != fKeys.end()};
    if(isRepeated)
        throw std::runtime_error("Fitters::Interface::CheckInit(): parameter " + key + " is repeated!");
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
        auto& initial {fInitial[key]};
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
                bounds.push_back({0, 0.5});
            // 3 -> Gamma
            if(i == 3)
                bounds.push_back({0, 2});
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
    auto& initial {fInitial.at(key)};
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
        throw std::runtime_error("Fitters::Interface::GetIdxKey() : cannot locate idx for " + key);
    auto idx {std::distance(fKeys.begin(), it)};
    return idx;
}

double Fitters::Interface::GetGuess(const Key& key) const
{
    return fGuess.at(key);
}

void Fitters::Interface::Print() const
{
    std::cout << BOLDCYAN << "----- Fitters::Interface -----" << '\n';
    if(!fIsRead)
    {
        for(const auto& key : fKeys)
        {
            const auto& initial {fInitial.at(key)};
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

void Fitters::Interface::Read(const std::string& file)
{
    auto f {std::make_unique<TFile>(file.c_str())};
    auto* inter {f->Get<Fitters::Interface>("inter")};
    if(!inter)
        throw std::runtime_error("Fitters::Interface::Read(): cannot read written object");
    *this = *inter;
    fIsRead = true;
    delete inter;
    // Clean states that are not gaussian or voigt
    CleanNotStates();
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
    auto str {key + " = " + FormatLabel(fLabels[key])};
    fComparators[key] = Angular::Comparator {str, gexp};
}

void Fitters::Interface::ReadComparatorConfig(const std::string& file)
{
    // Code partially generated with chatgpt
    // We could have used ActRoot::InputParser but I dont want
    // to introduce library dependences
    std::ifstream streamer(file);
    if(!streamer)
        throw std::runtime_error("Fitters::Interface::ReadComparator: cannot open config file");

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

template <typename T>
std::optional<T> Fitters::Interface::GetCompDrawOpt(const std::string& opt, const std::string& header)
{
    if(!fCompConf.count(header))
        return std::nullopt;
    auto& vec {fCompConf[header]};
    auto it {std::find_if(vec.begin(), vec.end(),
                          [&](const std::pair<std::string, std::string>& pair) { return pair.first == opt; })};
    if(it != vec.end())
        return std::optional<T> {(T)std::stod(it->second)};
    else
        return std::nullopt;
}

void Fitters::Interface::DoComp()
{
    // Get general
    auto logy {GetCompDrawOpt<bool>("logy").value()};
    auto withSF {GetCompDrawOpt<bool>("withSF").value()};
    auto offset {GetCompDrawOpt<double>("offset").value()};
    auto save {GetCompDrawOpt<bool>("save").value()};
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
        // If there are theoretical predictions
        if(fCompConf.count(state))
        {
            for(const auto& [key, file] : fCompConf[state])
            {
                // Check whether is file
                if(file.find("/") != std::string::npos)
                    comp.Add(key, file);
            }
            comp.Fit();
        }
        if(ip > npads)
        {
            ip = 1;
            ic++;
        }
        // cd into current pad (TVirtualPad)
        auto* base {cs[ic]->cd(ip)};
        // cast into TPad
        auto* pad {dynamic_cast<TPad*>(base)};
        // Draw!
        // But first check for local options!
        auto localLogy {GetCompDrawOpt<bool>("logy", state)};
        comp.Draw(state, (localLogy ? localLogy.value() : logy), withSF, offset, pad);
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
