# Release workflow

Neither the [master](https://github.com/revbayes/revbayes/tree/master) branch nor the [development](https://github.com/revbayes/revbayes/tree/development) branch is a direct ancestor of the other: each has at least one unique commit. Therefore we cannot merge either one into the other using a PR, because the PR will always say that the branch to be merged is out-of-date with respect to the base branch. What we can do is create our own commit that is an ancestor of both, instead of relying on GitHub PRs to merge things for us. But that has to be on a non-protected branch, since protected branches can only be modified by PRs.

1. Create a (non-protected) release branch

    1. Checkout the development branch
        * `git checkout development`
        * `git pull`
    2. Create the release branch off of development
        * `git branch release-vX.Y.Z`
        * `git checkout release-vX.Y.Z`
    3. Update `NEWS.md` file
    4. Make sure the release workflow (`release.yml`) is in sync with the build workflow (`build.yml`)

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

4. Merge development into master
    * On doing this, GitHub creates a merge commit, giving master 1 unique change that is not on development, and making this whole dance necessary next time

5. Update website from within the [revbayes.github.io](https://github.com/revbayes/revbayes.github.io) repo

    1. Checkout the source branch
        * `git checkout source`
        * `git pull`
    2. Create a release branch off of source
        * `git branch release-vX.Y.X`
        * `git checkout release-vX.Y.Z`
    3. Update download version by modifying `_config.yml`
    4. Re-build RevBayes from scratch. This automatically generates a `help.yml` file from the contents of [`help/md`](https://github.com/revbayes/revbayes/tree/master/help/md).
    5. Copy `help.yml` to the `_data` directory
    6. Commit the changes
        * `git add _config.yml`
        * `git add _data/help.yml`
        * `git commit -m "Update RB version to <X.Y.Z>"`
    7. Open a pull request to merge the release branch into source
