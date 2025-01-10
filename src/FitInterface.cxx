#include "FitInterface.h"

#include "PhysColors.h"

#include <algorithm>
#include <ios>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>


void Fitters::Interface::AddState(const Key& key, const DoubleVec& init)
{
    if(CheckInit(key, init))
    {
        fKeys.push_back(key);
        fInitial.insert({key, init});
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
        fGuess.push_back(fInitial[key][1]);
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
                bounds.push_back({fGuess[idx] - offset, fGuess[idx] + offset});
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
            fGuess[GetIdxKey(key)] = val;
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
    return fGuess.at(GetIdxKey(key));
}

void Fitters::Interface::Print() const
{
    std::cout << BOLDCYAN << "----- Fitters::Interface -----" << '\n';
    for(const auto& key : fKeys)
    {
        const auto& initial {fInitial.at(key)};
        const auto& bounds {fBounds.at(key)};
        const auto& fixed {fFix.at(key)};
        std::cout << "-> " << key << '\n';
        for(int i = 0; i < initial.size(); i++)
        {
            std::cout << "    " << fParNames[i] << '\n';
            std::cout << "      initial : " << initial[i] << '\n';
            std::cout << "      bounds  : [" << bounds[i].first << ", " << bounds[i].second << "]" << '\n';
            std::cout << "      fixed   ? " << std::boolalpha << fixed[i] << '\n';
        }
    }
    std::cout << RESET;
}
