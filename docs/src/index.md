# TransitionPathTheory.jl

## Introduction

This is a package for computing transition path theory[^1] (TPT) statistics a Markov chain on a discrete set $\mathbb{S}$.  All that is required to apply TPT is 

- A strongly connected, transition probability matrix $P$.
- A "source" set $\mathbb{A} \subset \mathbb{S}$ and a "target" set $\mathbb{B} \subset \mathbb{S}$. The set $\mathbb{A} \cap \mathbb{B}$ is automatically avoided by TPT.

The main function of TPT is to compute statistics of "reactive" trajectories, namely trajectories which travel directly from the source to the target with no intermediate visit to either.

### Features
- Automatic validation and fast computation with simple interface.
- Supports discrete, homogeneous (time-independent) Markov chains.
- Stationary and non-stationary TPT statistics.

## Documentation

[Documentation](https://70gage70.github.io/TransitionPathTheory.jl/)

## Installation

This package is in the Julia General Registry. In the Julia REPL, run the following code and follow the prompts:

```julia
import Pkg
Pkg.add("TransitionPathTheory")
```

Access the functionality of the package in your code by including the following line:

```julia
using TransitionPathTheory
```

## Quickstart

We will use a randomly generated stochastic transition matrix with $10$ states and source $\mathbb{A} = [1, 2, 3]$ and target $\mathbb{B} = [3, 4, 5]$.

```julia
P = TransitionMatrix(10)  
A = [1, 2, 3]                           
B = [3, 4, 5]                           
```              

Then, we set up `HomogeneousTPTProblem` and use `stationary_statistics` and `nonstationary_statistics` to compute all the relevant statistics in a `NamedTuple.`

```
tpt = HomogeneousTPTProblem(P, A, B)
stats = stationary_statistics(tpt)
stats_ns = nonstationary_statistics(tpt, 100)  # 100 time step horizon
```

## Citation

!!! note

    Please use the following citation if you use this package in your research.

    ```
    @article{bonner2023improving,
    title={Improving the stability of temporal statistics in transition path theory with sparse data},
    author={Bonner, Gage and Beron-Vera, FJ and Olascoaga, MJ},
    journal={Chaos: An Interdisciplinary Journal of Nonlinear Science},
    volume={33},
    number={6},
    year={2023},
    publisher={AIP Publishing}
    }
    ```

Initial development of this package was supported by the National Science Foundation.

[^1]: Vanden-Eijnden, Eric. "Transition path theory." Computer Simulations in Condensed Matter Systems: From Materials to Chemical Biology Volume 1. Springer, Berlin, Heidelberg, 2006. 453-493.
