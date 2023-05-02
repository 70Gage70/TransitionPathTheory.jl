# TransitionPathTheory.jl

## Introduction

This is a package for computing transition path theory[^fn1] (TPT) statistics of an ergodic Markov chain on a discrete set $\mathbb{S}$. All that is required to apply TPT is 

- A strongly connected, stochastic transition probability matrix $P$.
- A "source" set $\mathbb{A} \subset \mathbb{S}$ and a "target" set $\mathbb{B} \subset \mathbb{S}$. The set $\mathbb{A} \cap \mathbb{B}$ is automatically avoided by TPT.

The main function of TPT is to compute statistics of "reactive" trajectories, namely trajectories which travel directly from the source to the target with no intermediate visit to either.

## Installation

In the Julia REPL, run the following code and follow the prompts:

```julia
using Pkg
Pkg.add(TransitionPathTheory)
```

Make the package available to use in your code by including the following line:

```julia
using TransitionPathTheory
```

## Quickstart

```julia
P = P_stoc(10)      # a 10 x 10 stochastic matrix w/ random entries
A = [1, 2, 5]       # the source indices
B = [5, 9, 10]      # the target indices

# validate and setup a homogenous TPT problem
tpt_homog = TPTHomog(P, A, B)

# compute TPT statistics in the stationary case
tpt_stat = tpt_stationary_statistics(tpt_homog)
tpt_stat.q_plus # the forward committor; note that q_plus[5] = 0 
tpt_stat.normalized_reactive_density # normalized muAB

# compute TPT statistics in the nonstationary case
tpt_nonstat = tpt_nonstationary_statistics(tpt_homog)
tpt_nonstat.q_plus # the forward committor is the same for homogeneous problems
tpt_nonstat.normalized_reactive_density # normalized muAB; columns refer to increasing time

# write the results to an .h5 file
f_name = "tpt_results.h5"
tpt_write(f_name, tpt_stat, dir_name = "tpt_stat")
tpt_write(f_name, tpt_nonstat, dir_name = "tpt_nonstat")
```


[^fn1]: Vanden-Eijnden, Eric. "Transition path theory." Computer Simulations in Condensed Matter Systems: From Materials to Chemical Biology Volume 1. Springer, Berlin, Heidelberg, 2006. 453-493.
