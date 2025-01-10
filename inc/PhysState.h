#ifndef PhysState_h
#define PhysState_h

#include <cmath>
#include <string>
namespace PhysUtils
{
class State
{
private:
    float fJ {};          //!< Spin J of the state
    int fPi {};           //!< Parity
    unsigned int fIdx {}; //!< Repetition in case of rotational bands

public:
    State() = default;
    State(double j, int p, unsigned int idx = 0) : fJ(j), fPi(p), fIdx(idx) {}

    // Getters
    double J() const { return fJ; }
    int Pi() const { return fPi; }
    int Idx() const { return fIdx; }
    std::string Format() const
    {
        // Format spin
        std::string jstr {};
        if(std::floor(fJ) == fJ) // then is integer
            jstr = std::to_string((int)fJ);
        else
        {
            auto num {(int)(fJ * 2)}; // spin is always integer / 2
            jstr = std::to_string(num) + "/2";
        }
        // Format parity
        std::string pstr {fPi < 0 ? "^{-}" : "^{+}"};
        // Format idx
        std::string istr {fIdx > 0 ? (std::string("_{") + std::to_string(fIdx) + "}") : ""};
        return jstr + pstr + istr;
    }
};
} // namespace PhysUtils

#endif
