#include "TROOT.h"
//#include "readMVAtree.C"

void readMVATreeMacroFull()
{
   gROOT->ProcessLine(".L readMVAtreeFull.C++");
   gROOT->ProcessLine("readMVAtreeFull(0)");
//   gROOT->ProcessLine("readMVAtree(1)");
//   gROOT->ProcessLine("readMVAtree(2)");
}
