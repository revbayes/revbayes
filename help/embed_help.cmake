# Writes each help page in HELP_DIR as a { name, text } entry of C++ raw string literals.
# Usage: cmake -DHELP_DIR=<dir> -DOUT=<file> -P embed_help.cmake

file(GLOB pages "${HELP_DIR}/*.md")

set(text "")
foreach(page ${pages})
    get_filename_component(name "${page}" NAME)
    string(REGEX REPLACE "\\.md$" "" name "${name}")
    file(READ "${page}" content)

    string(FIND "${content}" ")rbhelp\"" found)
    if(NOT found EQUAL -1)
        message(FATAL_ERROR "${page} contains the string delimiter )rbhelp\"")
    endif()

    string(APPEND text "{ \"${name}\", R\"rbhelp(${content})rbhelp\" },\n")
endforeach()

file(WRITE "${OUT}" "${text}")
