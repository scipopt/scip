{
  description = "SCIP Optimization Suite";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    soplex-src = {
      url = "github:scipopt/soplex";
      flake = false;
    };
    papilo-src = {
      url = "github:scipopt/papilo";
      flake = false;
    };
    # just to get zimpl, would be awesome zimpl be just git repo in future:)
    suite-src = {
      url = "https://scipopt.org/download/release/scipoptsuite-8.1.0.tgz";
      flake = false;
    };
  };

  outputs = inputs@{ flake-parts, suite-src, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [ "x86_64-linux" "aarch64-linux" "aarch64-darwin" "x86_64-darwin" ];
      perSystem = { config, self', inputs', pkgs, system, ... }: {
        packages =
          let
            # need cmake support with binayr libs be same folder as dev texts of very specific version
            tbb-all-cmake = pkgs.tbb_2020_3.overrideAttrs (old: rec {
              doCheck = false;
              name = "tbb-all-cmake";
              outputs = [ "out" ];
              installPhase = old.installPhase + ''            
                mkdir --parents "$out"/lib/cmake/TBB
                ${pkgs.cmake}/bin/cmake \
                    -DINSTALL_DIR="$out"/lib/cmake/TBB \
                    -DSYSTEM_NAME=Linux -DTBB_VERSION_FILE="$out"/include/tbb/tbb_stddef.h \
                    -P cmake/tbb_config_installer.cmake
                cp --recursive --force ${pkgs.tbb_2020_3.dev}/include "$out"
                cp --recursive --force ${pkgs.tbb_2020_3.dev}/lib/pkgconfig "$out"/lib
                cp --recursive --force ${pkgs.tbb_2020_3.dev}/nix-support "$out"
              '';
            });
            scip-src = pkgs.nix-gitignore.gitignoreSource [ ] ./.;
            # all is here, really need to fine tune per package
            build-inputs = [
              pkgs.boost.dev
              pkgs.cmake
              pkgs.gmp.dev
              pkgs.pkg-config
              pkgs.gfortran
              pkgs.blas
              pkgs.bison
              pkgs.flex
              pkgs.libamplsolver
            ];

          in
          rec {
            soplex = pkgs.stdenv.mkDerivation {
              name = "soplex";
              src = inputs.soplex-src;
              buildInputs = build-inputs;
              nativeBuildInputs = build-inputs;
              doCheck = false;
            };
            papilo = pkgs.stdenv.mkDerivation {
              name = "papilo";
              src = inputs.papilo-src;
              nativeBuildInputs = build-inputs;
              buildInputs = build-inputs ++ [ tbb-all-cmake soplex ];
              doCheck = false;
            };
            zimpl = pkgs.stdenv.mkDerivation {
              name = "zimpl";
              src = "${suite-src}/zimpl";
              nativeBuildInputs = build-inputs;
              buildInputs = build-inputs;
              doCheck = false;
            };
            default = scip;
            inherit tbb-all-cmake;
            scip =
              pkgs.stdenv.mkDerivation {
                name = "scip";
                src = scip-src;
                cmakeFlags = [
                  "-D AMPL=on"
                  "-D AUTOBUILD=off" # by no means in nix can use internet and scan outer host system
                  "-D CMAKE_BUILD_TYPE=Release"
                  "-D CMAKE_CXX_COMPILER_ID=GNU"
                  "-D COVERAGE=off"
                  "-D DEBUGSOL=on"
                  "-D GMP=on"
                  "-D IPOPT=on" # really it has no NLP parts, which all behind auth/name/pay-walls even to download as of now
                  "-D LPS=spx"
                  "-D LPSCHECK=off"
                  "-D OPT=opt"
                  "-D PAPILO=on"
                  "-D READLINE=on"
                  "-D STATIC_GMP=on"
                  "-D TBB_DIR=${tbb-all-cmake}"
                  "-D THREADSAFE=off"
                  "-D WORHP=off"
                  "-D ZIMPL_DIR=${zimpl}"
                  "-D ZLIB=on"
                ];
                dontFixCmake = true;
                dontUseCmakeConfigure = false;
                nativeBuildInputs = build-inputs ++ [ papilo soplex tbb-all-cmake pkgs.criterion ];
                buildInputs = build-inputs ++ [
                  pkgs.bliss
                  pkgs.readline.dev
                  pkgs.zlib.dev
                  papilo
                  soplex
                  pkgs.ipopt
                  tbb-all-cmake
                  pkgs.criterion.dev
                ];

                doCheck = false;
              };
          };
      };
    };
}
