#ifndef PhysColors_h
#define PhysColors_h

// This file just defines the preprocessor variables allowing
// colorized output in std::cout
#include "TColor.h"

#include <string>
#include <unordered_map>
#include <vector>
#define RESET "\033[0m"
#define BLACK "\033[30m"              /* Black */
#define RED "\033[31m"                /* Red */
#define GREEN "\033[32m"              /* Green */
#define YELLOW "\033[33m"             /* Yellow */
#define BLUE "\033[34m"               /* Blue */
#define MAGENTA "\033[35m"            /* Magenta */
#define CYAN "\033[36m"               /* Cyan */
#define WHITE "\033[37m"              /* White */
#define BOLDBLACK "\033[1m\033[30m"   /* Bold Black */
#define BOLDRED "\033[1m\033[31m"     /* Bold Red */
#define BOLDGREEN "\033[1m\033[32m"   /* Bold Green */
#define BOLDYELLOW "\033[1m\033[33m"  /* Bold Yellow */
#define BOLDBLUE "\033[1m\033[34m"    /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m" /* Bold Magenta */
#define BOLDCYAN "\033[1m\033[36m"    /* Bold Cyan */
#define BOLDWHITE "\033[1m\033[37m"   /* Bold White */

// Class to hold custom colors
namespace PhysUtils
{
class Colors
{
private:
    std::vector<TColor*> fColors {};                                    //!< General colors
    std::unordered_map<std::string, std::vector<TColor*>> fPalettes {}; //!< Palette colors
    static Colors* fInstance;

    Colors() { Init(); }

    void Init();

    void AddPalette(const std::string& name, const std::vector<std::vector<int>>& vrgb, double norm = 255);
    void AddMatplotlib();

    void DrawFrom(const std::vector<TColor*>& colors, const std::string& title) const;

public:
    static Colors* GetInstance();

    Colors(Colors&) = delete;
    void operator=(const Colors&) = delete;

    int Get(int i, const std::string& pal = "");
    void Draw() const;
};
} // namespace PhysUtils

// Define a global variable
#define gPhysColors (PhysUtils::Colors::GetInstance())

#endif
