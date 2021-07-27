{
  description = ''
    This project provides several additional classes for the LHCb
    `Analysis/Phys/DecayTreeTuple` package for some semileptonic analyses.
  '';

  inputs = {
    root-curated.url = "github:umd-lhcb/root-curated";
    nixpkgs.follows = "root-curated/nixpkgs";
    flake-utils.follows = "root-curated/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils, root-curated }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; };
          overlays = [
            root-curated.overlay
          ];
        };
      in
      {
        devShell = pkgs.mkShell {
          name = "TupleToolSemiLeptonic";
          buildInputs = with pkgs; [
            clang-tools # For clang-format
          ];

          shellHook = ''
            export PATH=$(pwd)/test:$(pwd)/tools:$(pwd)/ganga:$(pwd)/scripts:$PATH
          '';
        };
      });
}
