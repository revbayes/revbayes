#!/usr/bin/env python3
# Writes each help page as a { name, text } entry of C++ raw string literals.
# Usage: embed_help.py <out> <page.md>...

import pathlib
import sys

out = pathlib.Path(sys.argv[1])
pages = sorted(pathlib.Path(p) for p in sys.argv[2:])

text = ""
for page in pages:
    content = page.read_text()
    if ')rbhelp"' in content:
        sys.exit(f'{page} contains the string delimiter )rbhelp"')
    text += f'{{ "{page.stem}", R"rbhelp({content})rbhelp" }},\n'

out.write_text(text)
