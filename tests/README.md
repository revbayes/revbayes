# RevBayes Integration Tests

This repository contains the integration test suite for RevBayes. Individual tests are stored in directories beginning with `test_` and are run automatically by the script `run_integration_tests.sh`. The return status of this script is used to determine whether the tests have passed.

## Website tutorials

The RevBayes website tutorials can also be tested using this script by placing a `test.sh` script in the appropriate tutorial directory on the [RevBayes website repository](https://github.com/revbayes/revbayes.github.io). The RevBayes website is used as a submodule of this repository and is automatically cloned and updated by the `run_integration_tests.sh` script. In order to clone and update it manually, you can use

```
git submodule update --init --recursive
```

This will clone the RevBayes website as a submodule directory. In order to test against a specific commit in the website repository, you should checkout that commit in the submodule.

## Updating the tests

Updating the test output or adding new tests requires both changes to this repository, as well as a change to the main RevBayes repository to update the version of the tests used by the main repository. Step-by-step instructions:

 * Make changes to the tests.
 * You need to add -DCONTINUOUS_INTEGRATION=TRUE to the arguments when compiling with cmake. This turns on using openlibm instead of the regular math library, and also adds the flag -mfpmath=sse . ./build.sh can also take -travis true which is a hold-over from when we used travis to do the tests.
 * Make sure you are on a branch of the `tests` repository that is separate from `development`. If necessary, checkout a new branch `git checkout -b <branch>`.
 * Add and commit changes to the tests repository using `git add` and `git commit` from within the `tests` folder.
 * Push the changes to the revbayes repository (using `git push`).
 * Update the tests version used by the main repository (using `git add tests`, `git commit` and `git push` from within the `revbayes` folder).
 * Note that the tests version is specific to each branch, so when merging with the main revbayes `development` branch, you should first merge the tests `development` into your test branch, then follow the above steps to update the tests submodule to point to your new merge commit.

