# Test the package using Github Actions

We use the [Github Actions](https://docs.github.com/en/actions) to perform the CI test of the package.

We made the test codes of the `SeisMonitoring.jl` using some examples in [**SeisMonitoring_Example**](https://github.com/kura-okubo/SeisMonitoring_Example) and [**SeisMonitoring_Paper**](https://github.com/kura-okubo/SeisMonitoring_Paper).

We test the followings:

1. build the package in different os and latest version of Julia.
2. synthetic test of dv/v measurement
3. ambient noise processings associated with the downloading, removing transient signals, computing the cross-correlation, and stacking.

Test scripts can be found in **[SeisMonitoring.jl/test
](https://github.com/kura-okubo/SeisMonitoring.jl/tree/dev_beforeupdate%400.2.0/test)**.

To run the test, type the following in the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.test("SeisMonitoring")
```

## Github Actions test.yml

Here is the script of github actions located in [`.github/workflows/test.yml`](https://github.com/kura-okubo/SeisMonitoring.jl/blob/dev_beforeupdate%400.2.0/.github/workflows/test.yml). You can configure this from the tab of `Actions` in the [web page of repository](https://github.com/kura-okubo/SeisMonitoring.jl/actions).

```yml
name: Run tests

on:
  push:
    branches:
      - master
      - dev_beforeupdate@0.2.0
  pull_request:

jobs:
  test:
    runs-on: ${{ matrix.os }}
    strategy:
      # We test the latest Julia v1.x, nightly in 64bit PC on ubuntu-latest and macOS-latest
      matrix:
        julia-version: ['1', 'nightly']
        julia-arch: [x64]
        os: [ubuntu-latest, macOS-latest]

    timeout-minutes: 20
    steps:
      - uses: actions/checkout@v2
      - uses: julia-actions/setup-julia@v1
        with:
          version: ${{ matrix.julia-version }}
          arch: ${{ matrix.julia-arch }}
      - run: julia --color=yes --project="@." -e 'using Pkg; ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; Pkg.Registry.add("General"); Pkg.develop(url="https://github.com/kura-okubo/SeisDvv.jl"); Pkg.instantiate(); Pkg.resolve(); Pkg.build()'
      - uses: julia-actions/julia-runtest@v1
        # with:
        #   annotate: true
```


!!! tip
    - `SeisMonitoring.jl` contains the unregistered package of `SeisDvv.jl`. This caused the error in building the package when we used `julia-actions/julia-buildpkg@v1` such as `ERROR: expected package SeisDvv to be registered`. Therefore, we did the hard coding of this process in the `steps` above, which solved this issue.

    - `ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0` saved the computational time of the test.

    - The usage of Github Actions runtime is limited like 2000-3000 minutes per month. Use `timeout-minutes: 20` to avoid the unexpected long time of executing the test.

You can check the status of this test by the badge as below:

[![Run tests](https://github.com/kura-okubo/SeisMonitoring.jl/actions/workflows/test.yml/badge.svg)](https://github.com/kura-okubo/SeisMonitoring.jl/actions/workflows/test.yml)

## Build environment

| os | architecture | Julia version |
| --- | --- | --- |
| ubuntu-latest | x64 | latest (1.9.3 as of Sep. 2023) |
| ubuntu-latest | x64 | nightly |
| macOS-latest | x64 | latest (1.9.3 as of Sep. 2023) |
| macOS-latest | x64 | nightly |

!!! note
    The architecture of x86 (32-bit machines) is not supported as we cannot build the package there. The Windows os (windows-latest) also did not work due to the different syntax of the commands. However, we locally tested the `SeisMonitoring.jl` using the Docker Desktop in Windows os (Windows 10). Please use the docker-desktop to use the `SeisMonitoring.jl` for Windows users. The instruction can be found in the `README.md` of [**SeisMonitoring_Example**](https://github.com/kura-okubo/SeisMonitoring_Example#how-to-run-the-notebook).

## Reference
- [Julia Actions](https://github.com/julia-actions): setup julia in the virtual machine and run the test.
