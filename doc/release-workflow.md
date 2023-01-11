# Release workflow

0. Update NEWS.md file

1. Build packages

    1. Checkout the development branch
        * `git checkout development`
        * `git pull`
    2. Create a tag for the release
        * `git tag vD.D.D`
        * tag names
            * v2.0 for a major release
    3. Push the tag to github
        * `git push --tags`
        * This will create a draft release.
        * GitHub actions will upload packages for Linux, Linux/Singularity, Mac, Windows
    4. Wait for packages to build.
    5. Edit draft-release page on github
        * See https://github.com/revbayes/revbayes/releases
        * Add some basic release notes.
        * Make the draft release public.

2. Update website

    1. Update version number on download page.
    2. Update links on download page.
