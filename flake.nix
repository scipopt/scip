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
    tbb-src = {
      url = github:oneapi-src/oneTBB/tbb_2020;
      flake = false;
    };
  };

  outputs = inputs@{ flake-parts, suite-src, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [ "x86_64-linux" "aarch64-linux" "aarch64-darwin" "x86_64-darwin" ];
      perSystem = { config, self', inputs', pkgs, system, ... }: {
        packages =
          let
            scip-src = pkgs.nix-gitignore.gitignoreSource [ ] ./.;
            build-inputs = [
              pkgs.boost.dev
              pkgs.cmake
              pkgs.gmp.dev
              pkgs.pkg-config
              pkgs.tbb.dev
            ];

          in
          rec {
            soplex = pkgs.stdenv.mkDerivation {
              name = "soplex";
              src = inputs.soplex-src;
              nativeBuildInputs = build-inputs;
              doCheck = false;
            };
            papilo = pkgs.stdenv.mkDerivation {
              name = "papilo";
              src = inputs.papilo-src;
              nativeBuildInputs = build-inputs;
              buildInputs = build-inputs;
              doCheck = false;
            };
            zimpl = pkgs.stdenv.mkDerivation {
              name = "zimpl";
              src = "${suite-src}/zimpl";
              nativeBuildInputs = build-inputs ++ [ pkgs.bison pkgs.flex ];
              buildInputs = build-inputs;
              doCheck = false;
            };
            default = scip;
            scip =
              let
                tbb-find = "${inputs.papilo-src}/cmake/Modules/FindTBB.cmake";
              in
              pkgs.stdenv.mkDerivation {
                name = "scip";
                src = scip-src;
                cmakeFlags = [
                  "-DZIMPL_DIR=${zimpl}"
                  "-DTBB_DIR=${inputs.tbb-src}/cmake"
                  "-DCMAKE_PREFIX_PATH=${tbb-find}"
                  "-DPAPILO=off"
                  "-DAUTOBUILD=off" # by no means in nix can use internet and scan outer host system
                ];
                dontFixCmake = true;
                dontUseCmakeConfigure = false;
                nativeBuildInputs = build-inputs ++ [ papilo soplex ];
                buildInputs = build-inputs ++ [
                  pkgs.bliss
                  pkgs.readline.dev
                  pkgs.zlib.dev
                  papilo
                  soplex
                  pkgs.ipopt
                ];
                doCheck = false;
              };
          };
      };
    };
}
