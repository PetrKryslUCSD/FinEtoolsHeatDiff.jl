[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build status](https://github.com/PetrKryslUCSD/FinEtoolsHeatDiff.jl/workflows/CI/badge.svg)](https://github.com/PetrKryslUCSD/FinEtoolsHeatDiff.jl/actions)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/FinEtoolsHeatDiff.jl/latest)

# FinEtoolsHeatDiff: Linear heat diffusion analysis

[`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl.git) is a package
for basic operations on finite element meshes. `FinEtoolsHeatDiff` is a package
using `FinEtools` to solve linear heat conduction (diffusion) problems.

## News

- 02/25/2024: Enables tutorials.
- 02/23/2024: Updated documentation.
- 02/21/2024: Updated for FinEtools 8.0.
- 09/28/2023: Updated examples.

[Past news](#past-news)

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
For instance
```
julia> include("steady_state/3-d/actuator_examples.jl"); Main.actuator_examples.allrun()  
```

## <a name="past-news"></a>Past news

- 06/20/2023: Updated for FinEtools 7.0.
- 05/20/2023: Updated for Julia 1.9.
- 04/22/2023: Updated for generic FinEtools.
- 03/08/2022: Introduced incompatible change of the assemblers (FinEtools 6.0.1).
- 05/05/2022: Renamed the branch `main`. Updated for Julia 1.7.
- 05/23/2021: Updated for Julia 1.6.
- 10/11/2019: Corrected design flaw in matrix utilities.
- 08/30/2019: A transient heat-conduction example added.
- 06/11/2019: Applications are now separated  out from the `FinEtools` package.