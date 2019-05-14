#include "TROOT.h"
//#include "readMVAtree.C"

void readMVATreeMacro()
{
   gROOT->ProcessLine(".L readMVAtree.C++");
   gROOT->ProcessLine("readMVAtree(0)");
   gROOT->ProcessLine("readMVAtree(1)");
   gROOT->ProcessLine("readMVAtree(2)");
}
