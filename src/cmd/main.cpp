#include <stdio.h>
#include <stdlib.h>

#include "RbGTKGui.h"


int main( int argc, char *argv[] )
{
    
    RevBayesGTK::RbGTKGui& gui = RevBayesGTK::RbGTKGui::globalInstanceGUI();
    gui.start(argc, argv);
    
    return 0;
}
