# TransitionPathTheory.jl

## Introduction

This is a package for computing transition path theory[^fn1] (TPT) statistics of an ergodic Markov chain on a discrete set $\mathbb{S}$. All that is required to apply TPT is 

- A strongly connected, stochastic (or sub-stochastic) transition probability matrix $P$.
- A "source" set $\mathbb{A} \subset \mathbb{S}$ and a "target" set $\mathbb{B} \subset \mathbb{S}$.

The main function of TPT is to compute statistics of "reactive" trajectories, namely trajectories which travel directly from the source to the target with no intermediate visit to either.

## Installation

_Under construction._

## Quickstart

```julia
P = P_stoc(10)  # a 10 x 10 stochastic matrix w/ random entries
A = [1, 2]      # the source indices
B = [9, 10]     # the target indices

# validate and setup a homogenous TPT problem
tpt_homog = TPTHomog(P, A, B)

# compute all the statistics using infinite time TPT
tpt_stats = tpt_infinite(tpt_homog)

# write the results to an .h5 file
tpt_write("tpt_stats.h5", tpt_homog, tpt_stats)
```


[^fn1]: Vanden-Eijnden, Eric. "Transition path theory." Computer Simulations in Condensed Matter Systems: From Materials to Chemical Biology Volume 1. Springer, Berlin, Heidelberg, 2006. 453-493.
