#include <cstdio>
#include <cstdlib>

#include "RbGTKGui.h"


int main( int argc, char *argv[] )
{
    
    RevBayesGTK::RbGTKGui& gui = RevBayesGTK::RbGTKGui::globalInstanceGUI();
    gui.start(argc, argv);
    
    return 0;
}
