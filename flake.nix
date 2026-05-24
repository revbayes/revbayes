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

        mkRevBayes = { withOpenMP ? false }:
          pkgs.stdenv.mkDerivation {
            pname = if withOpenMP then "revbayes-omp" else "revbayes";
            version = "1.4.0";

            src = ./.;

            nativeBuildInputs = with pkgs; [ cmake ninja ];
            buildInputs = with pkgs; [ boost zlib ]
              ++ pkgs.lib.optionals withOpenMP [ pkgs.llvmPackages.openmp ];

            cmakeDir = "../src";
            cmakeFlags = pkgs.lib.optionals withOpenMP [ "-DOPENMP=ON" ];

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

      in {
        packages.default     = mkRevBayes {};
        packages.revbayes-omp = mkRevBayes { withOpenMP = true; };

        devShells.default = pkgs.mkShell {
          nativeBuildInputs = with pkgs; [ cmake ninja ];
          buildInputs = with pkgs; [ boost zlib ];
          BOOST_INCLUDEDIR = "${pkgs.lib.getDev pkgs.boost}/include";
          BOOST_LIBRARYDIR = "${pkgs.lib.getLib pkgs.boost}/lib";
        };

        devShells.omp = pkgs.mkShell {
          nativeBuildInputs = with pkgs; [ cmake ninja ];
          buildInputs = with pkgs; [ boost zlib pkgs.llvmPackages.openmp ];
          BOOST_INCLUDEDIR = "${pkgs.lib.getDev pkgs.boost}/include";
          BOOST_LIBRARYDIR = "${pkgs.lib.getLib pkgs.boost}/lib";
        };
      });
}
