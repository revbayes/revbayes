{
  description = "RevBayes -- Bayesian phylogenetic inference using probabilistic graphical models";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in {
        packages.default = pkgs.stdenv.mkDerivation {
          pname = "revbayes";
          version = "1.4.2-preview";

          src = ./.;

          nativeBuildInputs = with pkgs; [ cmake ninja ];
          buildInputs = with pkgs; [ boost zlib ];

          # CMakeLists.txt lives in src/, not the repo root
          cmakeDir = "../src";

          preConfigure = ''
            # Generate GitVersion.cpp — git is not available in the Nix sandbox
            {
              echo '#include "GitVersion.h"'
              echo 'const char *build_git_sha = "1.4.0";'
              echo 'const char *build_date = "unknown";'
              echo 'const char *build_git_branch = "unknown";'
            } > src/revlanguage/utils/GitVersion.cpp

            # Generate generated_include_dirs.cmake and per-subdir CMakeLists.txt
            # (normally done by build.sh before invoking cmake)
            (cd projects/cmake && bash regenerate.sh)
          '';

          # CMake reads these when BOOST_INCLUDEDIR/BOOST_LIBRARYDIR are present
          BOOST_INCLUDEDIR = "${pkgs.lib.getDev pkgs.boost}/include";
          BOOST_LIBRARYDIR = "${pkgs.lib.getLib pkgs.boost}/lib";

          meta = with pkgs.lib; {
            description = "Bayesian phylogenetic inference using probabilistic graphical models";
            homepage = "https://revbayes.com";
            license = licenses.gpl2Only;
            mainProgram = "rb";
            platforms = platforms.unix;
          };
        };

        devShells.default = pkgs.mkShell {
          nativeBuildInputs = with pkgs; [ cmake ninja ];
          buildInputs = with pkgs; [ boost zlib ];

          BOOST_INCLUDEDIR = "${pkgs.lib.getDev pkgs.boost}/include";
          BOOST_LIBRARYDIR = "${pkgs.lib.getLib pkgs.boost}/lib";
        };
      });
}
