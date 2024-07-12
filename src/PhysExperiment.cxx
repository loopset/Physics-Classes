#include "PhysExperiment.h"

#include "PhysColors.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

void PhysUtils::Experiment::Write(const std::string& file) const
{
    std::ofstream streamer {file};
    streamer << "Ntargets(/cm2): " << fNt << '\n';
    streamer << "BeamTriggers: " << fNtriggers << '\n';
    streamer << "DivFactor: " << fNdiv << '\n';
    streamer << "Nb: " << fNb << '\n';
    streamer.close();
}

void PhysUtils::Experiment::Read(const std::string& file)
{
    std::ifstream streamer {file};
    std::string line {};
    int row {};
    while(std::getline(streamer, line))
    {
        std::istringstream it {line};
        std::string val {};
        int col {};
        while(std::getline(it, val, ' '))
        {
            if(col == 1)
            {
                if(row == 0)
                    fNt = std::stod(val);
                else if(row == 1)
                    fNtriggers = std::stod(val);
                else if(row == 2)
                    fNdiv = std::stod(val);
                else if(row == 3)
                    fNb = std::stod(val);
                else
                    continue;
            }
            col++;
        }
        row++;
    }
}

void PhysUtils::Experiment::Print() const
{
    std::cout << BOLDYELLOW;
    std::cout << "---- PhysUtils::Experiment ----" << '\n';
    std::cout << "Ntargets(/cm2): " << fNt << '\n';
    std::cout << "BeamTriggers: " << fNtriggers << '\n';
    std::cout << "DivFactor: " << fNdiv << '\n';
    std::cout << "Nb: " << fNb << '\n';
    std::cout << RESET << '\n';
}
