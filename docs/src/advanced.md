# Advanced Usage

```@contents
Pages = ["advanced.md"]
Depth = 5
```

# High level overview

Use of this package proceeds as follows.

1. Define a `TransitionMatrix`, and provide source, target and avoid indices as desired.
2. Select a `TPTProblem`.
3. Compute the relevant statistics.

# Defining transition matrices

```@docs; canonical=false
TransitionMatrix
```

## Stochasticity

```@docs; canonical=false
Stochasticity
Stochastic
SuperStochastic
NonStochastic
```

## Connectivity

```@docs; canonical=false
Connectivity
StronglyConnected
WeaklyConnected
Disconnected
```

# Defining a TPTProblem

```@docs; canonical=false
TPTProblem
HomogeneousTPTProblem
```

## Basic sets in a TPTProblem

```@docs; canonical=false
𝒫
𝒜
ℬ
𝒮
Ω
𝒜_true
ℬ_true
𝒞
```

# TPT Statistics

## High level

```@docs; canonical=false
stationary_statistics
nonstationary_statistics
```

## Low level

```@docs; canonical=false
stationary_distribution
𝒫_backwards
forward_committor
backward_committor
𝒮_plus
𝒫_plus
remaining_time
hitting_location_distribution
```