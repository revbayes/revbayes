# Release workflow

Neither the [master](https://github.com/revbayes/revbayes/tree/master) branch nor the [development](https://github.com/revbayes/revbayes/tree/development) branch is a direct ancestor of the other: each has at least one unique commit. Therefore we cannot merge either one into the other using a PR, because the PR will always say that the branch to be merged is out-of-date with respect to the base branch. What we can do is create our own commit that is an ancestor of both, instead of relying on GitHub PRs to merge things for us. But that has to be on a non-protected branch, since protected branches can only be modified by PRs.

1. Create a (non-protected) release branch

    1. Checkout the development branch
        * `git checkout development`
        * `git pull`
    2. Create the release branch off of development
        * `git branch release-vX.Y.Z`
        * `git checkout release-vX.Y.Z`
    3. Update the version number
        * modify `meson.build`
        * modify `src/revlanguage/utils/RbVersion.cpp`
        * modify `flake.nix`
    4. Update `NEWS.md` file
    5. Make sure the release workflow (`release.yml`) is in sync with the build workflow (`build.yml`)

2. Prepare the release

    1. Merge master into the release branch using git
        * `git merge master`
    2. Create a tag for the release
        * `git tag vX.Y.Z`
        * tag names: no strict semantic versioning; use v2.0 for a major release
    3. Push the tag to github
        * `git push --tags`
        * This will create a draft release.
        * GitHub actions will upload packages for Linux, Linux/Singularity, macOS, Windows
    4. Wait for packages to build.
    5. Edit the draft release page on GitHub
        * See https://github.com/revbayes/revbayes/releases
        * Add release notes from `NEWS.md`
        * Make the draft release public
        
3. Open a pull request to merge the release branch into development, making development a strict descendant of master

4. On development, update the URLs for the following badges in `README.md`:
    * New commits in last release
    * Downloads

5. Merge development into master
    * On doing this, GitHub creates a merge commit, giving master 1 unique change that is not on development, and making this whole dance necessary next time
    
6. Check the website

    Publishing the release runs `notify-website.yml`, which starts the `update-help` workflow in [revbayes.github.io](https://github.com/revbayes/revbayes.github.io). That workflow downloads the release's linux64 package, regenerates `_data/help.json` with `rb --dump-help-json`, sets the version in `_config.yml`, and opens a pull request with both that merges by itself once the website's checks pass. If the workflow fails, rerun it from that repository's Actions tab with the release tag. If a check fails, the pull request stays open until it is fixed.
