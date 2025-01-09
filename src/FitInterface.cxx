#include "FitInterface.h"

#include <stdexcept>
#include <string>


void Fitters::Interface::Init(const Key& key, const DoubleVec& init)
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
    if(type != "g" || type != "v" || type != "ps" || type != "cte")
        throw std::invalid_argument("Fitters::Interface::GetType(): received wrong key " + key);
    return type;
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
                bounds.push_back({0, 10000});
            // 1 -> Mean
            if(i == 1)
            {
                double offset {0.2};
                bounds.push_back({fInitialGuess[idx] - offset, fInitialGuess[idx] + offset});
            }
            // 2 -> Sigma
            if(i == 3)
                bounds.push_back({0, 0.5});
            // 3 -> Gamma
            if(i == 4)
                bounds.push_back({0, 2});
        }
        fBounds.insert({key, bounds});

        // Fixed
        BoolVec fixed(fNParsType[type], false); // released by default
        fFix.insert({key, fixed});

        idx++;
    }
}
