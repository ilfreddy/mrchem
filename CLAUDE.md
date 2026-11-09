# CLAUDE.md

Guidance for Claude Code when working in the MRChem repository.

MRChem is a multiwavelet, real-space code for molecular electronic structure
(HF and DFT) built on top of the MRCPP library. It is a C++17 core (`src/`) with
a Python front end (`python/`) that parses the user input and launches the
binary.

## Configure and build

Use the `setup` script, which wraps CMake (run `./setup --help` for all options):

```bash
./setup --type=release --omp build          # OpenMP build in ./build
./setup --type=debug --mpi --omp build-mpi  # hybrid MPI/OpenMP debug build
./setup --show build                        # print the CMake command, don't run it
cmake --build build -j
```

Useful options: `--cxx=<compiler>`, `--type=debug|release|relwithdebinfo`,
`--arch-flags=false` (portable binaries, as used in CI), `--generator=Ninja`,
`--coverage`, `--prefix=<path>`, `--cmake-options="-D..."`.
These map to the CMake options `ENABLE_OPENMP` and `ENABLE_MPI` (both OFF by default).

### Dependencies

Fetched by `external/upstream/fetch_*.cmake` if no suitable installed version is
found: **MRCPP**, XCFun, Libxc, Eigen3 and nlohmann_json.

- MRCPP is pinned to a specific commit in `external/upstream/fetch_mrcpp.cmake`.
  Point CMake at a local install with `-DMRCPP_DIR=<prefix>/share/cmake/MRCPP`.
- An installed MRCPP must match MRChem's parallel configuration: MPI MRChem with
  non-MPI MRCPP is a configuration error; OpenMP MRChem with non-OpenMP MRCPP only
  gives a warning.
- Externals are always built in Release mode, whatever MRChem's build type is.

## Running the tests

From the build directory:

```bash
ctest --output-on-failure            # everything
ctest -L unit                        # Catch2 unit tests (mrchem-tests binary)
ctest -L integration                 # full mrchem runs compared to reference output
ctest -L short                       # quick subset; other labels include medium, scf, dft, xcfun, libxc
ctest -R H2_SCF_HF                   # a single test by name
ctest -j4                            # in parallel
```

- Unit tests live in `tests/qmfunctions`, `tests/qmoperators` and `tests/solventeffect`,
  linked into `mrchem-tests` (Catch2 v3).
- Each integration test is a directory under `tests/` containing an input, a
  `reference/` folder and a `test` script that uses `tests/tester.py` with
  per-quantity tolerances. They are registered in `tests/CMakeLists.txt` via the
  `add_integration_test` macro (`cmake/custom/testing_macros.cmake`).
- In MPI builds the tests run through `${MPIEXEC}` with a single process.

## Code conventions

- **C++17**, set in `cmake/compiler_flags/CXXFlags.cmake`.
- **Formatting:** clang-format using the repo `.clang-format` (4-space indent).
  Run it on changed C++ files before committing. The hooks in `.githooks/`
  (clang-format plus the license header maintainer) can be enabled with
  `git config core.hooksPath .githooks`.
- **License headers:** every source file has the LGPL header (`.githooks/LICENSE-C++`,
  `.githooks/LICENSE-Python`). Keep it on new files.
- **Naming:** code lives in `namespace mrchem`. Classes use `CamelCase` and have
  one class per file with a matching name (`Orbital.h`/`Orbital.cpp`). Free-function
  collections use `snake_case` files and namespaces (`orbital_utils`). Methods are
  `camelCase` (`getCharge`, `printGeometry`). Members are accessed via `this->`.
- **Documentation:** Doxygen comments (`/** @brief ... */`) in the C++ sources and Sphinx/RST in `doc/`.
- **Parallelism:** guard MPI/OpenMP-specific code with the `MRCHEM_HAS_MPI` /
  `MRCHEM_HAS_OMP` compile definitions.

<!-- TODO: add any further naming or style rules (e.g. member naming, header include order) -->

## How the Python parts fit in

- `python/mrchem/` is the `mrchem` driver script installed alongside the
  `mrchem.x` binary. It reads the user's `.inp` file, validates it, applies
  defaults and computes derived values (`validators.py`, `helpers.py`), and
  writes the JSON program input that the C++ binary consumes (`src/driver.cpp`).
- The input parser in `python/mrchem/input_parser/` is **generated** by
  `parselglossy` from `python/template.yml`, which is the single source of truth for input
  keywords, types, defaults and their documentation.
- Physical constants and element data come from `physical_constants.py` and
  `periodictable.py`.
- See `python/mrchem/input_parser/README.md` for how to regenerate the parser.

## Commits and pull requests

If you make a commit, follow conventional commits and add a trailer:
`Assisted-by: <harness>:<model>`, where `<harness>` is the current agent harness
(like ClaudeCode), and `<model>` is the AI model (Like claude-opus-4.8). You
don't need to add a coauthored-by Claude when you have this.

Prefix PR descriptions and comments on PRs with the line ":robot: _AI text
below_ :robot:" to indicate you are an agent speaking on a user's behalf.

## Things to avoid

- **Don't hand-edit generated files.** Change `python/template.yml` and regenerate
  with `python/mrchem/update_input_parser.py` instead of editing:
  - `python/mrchem/input_parser/` (`api.py`, `cli.py`, `plumbing/`, `docs/user_ref.rst`)
  - `doc/users/user_ref.rst`
- **Don't change numerical defaults without discussion.** This covers defaults in
  `template.yml`, thresholds and precisions, physical constants and integration
  test tolerances. These affect results and reference values.
- **Don't update test reference files** just to make a test pass. A changed
  reference needs an explanation of why the numbers moved.
- **Don't touch vendored or upstream-managed code:** `python/mrchem/input_parser/plumbing/pyparsing`,
  `cmake/autocmake` and `cmake/downloaded` (managed by `cmake/update.py`).
- **Don't change `VERSION` or tags.** Code administrators handle versioning (see `CONTRIBUTING.md`).
- **Don't bump the pinned MRCPP commit** in `fetch_mrcpp.cmake` without discussion.

<!-- TODO: add further items -->
