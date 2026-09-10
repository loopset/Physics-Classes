#ifndef FitGlobals_h
#define FitGlobals_h

namespace Fitters
{

extern bool gIsLogL; //!< Global extern variable to represent whether ALL (including in Angular::) count fits are done
                     //!< with Chi2 (default) or LogL

void ToggleIsLogL();

bool GetIsLogL();

} // namespace Fitters
#endif // !FitGlobals_h
