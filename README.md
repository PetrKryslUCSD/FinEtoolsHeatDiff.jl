[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build status](https://github.com/PetrKryslUCSD/FinEtoolsHeatDiff.jl/workflows/CI/badge.svg)](https://github.com/PetrKryslUCSD/FinEtoolsHeatDiff.jl/actions)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/FinEtoolsHeatDiff.jl/dev)

# FinEtoolsHeatDiff: Linear heat diffusion analysis

[`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl.git) is a package
for basic operations on finite element meshes. `FinEtoolsHeatDiff` is a package
using `FinEtools` to solve linear heat conduction (diffusion) problems.

## News

- 05/20/2023: Updated for Julia 1.9.
- 04/22/2023: Updated for generic FinEtools.
- 03/08/2022: Introduced incompatible change of the assemblers (FinEtools 6.0.1).
- 05/05/2022: Renamed the branch `main`. Updated for Julia 1.7.
- 05/23/2021: Updated for Julia 1.6.
- 10/11/2019: Corrected design flaw in matrix utilities.

[Past news](oldnews.md)

## How to test the package

Here is a record of a session to install this package and test it. You should
see something similar. The git bash running on Windows 10 was used in this
example.

Clone the repo:
```
PetrKrysl@Spectre MINGW64 /tmp/exp
$ git clone https://github.com/PetrKryslUCSD/FinEtoolsHeatDiff.jl.git
Cloning into 'FinEtoolsHeatDiff.jl'...
[Further output suppressed...]
```
Change your working directory, and run Julia:
```
$ cd FinEtoolsHeatDiff.jl/
PS julia-1.9.0\bin\julia.exe
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.0 (2023-05-07)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |
```
Activate and instantiate the environment:
```
(@v1.9) pkg> activate .; instantiate
  Activating project at `C:\Users\...\FinEtoolsHeatDiff.jl`
```
Test the package:
```
(FinEtoolsHeatDiff) pkg> test
     Testing FinEtoolsHeatDiff
[Output suppressed...]

     Testing Running tests...
Test Summary:  | Pass  Total     Time
Heat diffusion |   62     62  1m08.7s
 68.849887 seconds (85.11 M allocations: 6.588 GiB, 2.97% gc time, 85.82% compilation time)
     Testing FinEtoolsHeatDiff tests passed

(FinEtoolsHeatDiff) pkg>
```

## Examples

The examples have their own environment. Change the folder to `examples`.
Then activate and instantiate the `examples` environment.
```
(FinEtoolsHeatDiff) pkg>

shell> cd examples
C:\Users\...\FinEtoolsHeatDiff.jl\examples

julia> using Pkg

julia> Pkg.activate("."); Pkg.instantiate()
  Activating project at `C:\Users\...\FinEtoolsHeatDiff.jl\examples`
   [Output suppressed...]

julia>
```

There are a number of examples, which may
be executed as described in the  [conceptual guide to
`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest).
