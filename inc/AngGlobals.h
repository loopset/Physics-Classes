#ifndef AngGlobals_h
#define AngGlobals_h

namespace Angular
{

extern bool gIsLab; //!< Global extern variable to represent whether xs calculations are in CM or Lab

void ToggleIsLab();

bool GetIsLab();

} // namespace Angular
#endif // !AngGlobals_h
