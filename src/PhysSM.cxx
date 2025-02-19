#include "PhysSM.h"

#include "TString.h"

#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

std::string PhysUtils::QuantumNumbers::Format() const
{
    auto str {TString::Format("%d/2%s", static_cast<int>(2 * fJ), (fL % 2 == 0) ? "+" : "-")};
    return str.Data();
}

PhysUtils::SMParser::SMParser(const std::vector<std::string>& files)
{
    for(const auto& file : files)
        ParseFile(file);
}

std::string PhysUtils::SMParser::StripSpaces(std::string line)
{
    // Remove preceding spaces
    while(*line.begin() == ' ')
        line = line.substr(1, line.length());
    // Remove trailing spaces
    if(line.length() > 0)
        while(*line.rbegin() == ' ')
            line = line.substr(0, line.length() - 1);
    // Remove preceding tabs
    while(*line.begin() == '\t')
        line = line.substr(1, line.length());
    // Remove trailing tabs
    if(line.length() > 0)
        while(*line.rbegin() == '\t')
            line = line.substr(0, line.length() - 1);
    return line;
}

void PhysUtils::SMParser::ParseFile(const std::string& file)
{
    std::ifstream streamer {file};
    if(!streamer)
        throw std::invalid_argument("ModelParser::ParseFile(): cannot load file " + file);
    std::string line {};
    QuantumNumbers currentQuantum {};
    while(std::getline(streamer, line))
    {
        line = StripSpaces(line);
        if(!line.length())
            continue;
        // Use TString to check contains
        TString tstr {line};
        if(tstr.Contains("orbit"))
        {
            auto n {line.substr(9, 2)};
            auto l {line.substr(12, 2)};
            auto j {line.substr(15, 2)};
            for(auto str : {&n, &l, &j})
                *str = StripSpaces(*str);
            currentQuantum = QuantumNumbers(std::stoi(n), std::stoi(l), std::stod(j) / 2);
        }
        else if(tstr.BeginsWith("0("))
        {
            auto ex {line.substr(34, 7)};
            auto sf {line.substr(44)};
            for(auto str : {&ex, &sf})
                *str = StripSpaces(*str);
            SMData data {std::stod(ex), std::stod(sf)};
            data.AddQuantumNumbers(currentQuantum);
            fMap[currentQuantum].insert(data);
        }
    }
}

void PhysUtils::SMParser::ShiftEx()
{
    // Find states with largest SF == ground state
    double max {-1111};
    double shift {};
    for(const auto& [key, vals] : fMap)
        for(const auto& val : vals)
            if(val.fSF > max)
            {
                max = val.fSF;
                shift = val.fEx;
            }
    // Print
    // Shift
    for(auto& [key, vals] : fMap)
    {
        std::set<SMData> newSet {};
        for(auto& val : vals)
        {
            auto newData {val};
            newData.fEx -= shift;
            newSet.insert(newData);
        }
        fMap[key] = newSet;
    }
}

void PhysUtils::SMParser::MaskSFBelow(double thres)
{
    for(auto& [key, vals] : fMap)
    {
        for(auto it = vals.begin(); it != vals.end();)
        {
            if(it->fSF < thres)
                it = vals.erase(it);
            else
                it++;
        }
    }
}

void PhysUtils::SMParser::MaskExAbove(double thres)
{
    for(auto& [key, vals] : fMap)
    {
        for(auto it = vals.begin(); it != vals.end();)
        {
            if(it->fEx > thres)
                it = vals.erase(it);
            else
                it++;
        }
    }
}
