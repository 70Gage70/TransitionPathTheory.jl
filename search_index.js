var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"These are the full docstrings for TransitionPathTheory.jl.","category":"page"},{"location":"api/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Pages = [\"api.md\"]","category":"page"},{"location":"api/#Functions","page":"API","title":"Functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [TransitionPathTheory]\nOrder = [:function]","category":"page"},{"location":"api/#TransitionPathTheory.backward_committor-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.backward_committor","text":"backward_committor(tpt)\n\nCompute the backward committor.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.forward_committor-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.forward_committor","text":"forward_committor(tpt)\n\nCompute the forward committor.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.hitting_location_distribution-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.hitting_location_distribution","text":"hitting_location_distribution(tpt)\n\nCompute the distribution of target-hitting locations. \n\nReturn a matrix R such that R[i, j] is the probability that - starting in state i - the first visit to B occurs in state j. Note that  R[i, j] == 0 for all j ∉ B.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.nonstationary_statistics-Tuple{HomogeneousTPTProblem, Integer}","page":"API","title":"TransitionPathTheory.nonstationary_statistics","text":"nonstationary_statistics(tpt, horizon)\n\nArguments\n\ntpt: The TPTProblem.\nhorizon: The time step at which to cut off the calculation. Note that the horizon value does NOT enforce that trajectories leaving A hit B by that time.\n\nCompute and return the following statistics in a NamedTuple:\n\ndensity\nreactive_density\nt_cdf\nt_cdf_AB\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.remaining_time-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.remaining_time","text":"remaining_time(tpt)\n\nCompute the remaining time, aka the lead time.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.stationary_distribution-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.stationary_distribution","text":"stationary_distribution(tpt)\n\nCompute the stationary distribution of tpt.transition_matrix.P.\n\nErrors when P is not strongly connected.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.stationary_statistics-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.stationary_statistics","text":"stationary_statistics(tpt)\n\nCompute and return the following statistics in a NamedTuple:\n\nstationary_distribution\nforward_committor\nbackward_committor\nreactive_density\nreactive_current\nforward_current\nreactive_rate\nreactive_time\nremaining_time\nhitting_location_distribution\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.Ω-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.Ω","text":"Ω(tpt)\n\nReturn the set of avoided states.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.ℬ-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.ℬ","text":"ℬ(tpt)\n\nReturn tpt.target.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.ℬ_true-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.ℬ_true","text":"ℬ_true(tpt)\n\nReturn the set of states that are members of the target, but are not avoided.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.𝒜-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.𝒜","text":"𝒜(tpt)\n\nReturn tpt.source.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.𝒜_true-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.𝒜_true","text":"𝒜_true(tpt)\n\nReturn the set of states that are members of the source, but are not avoided.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.𝒞-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.𝒞","text":"𝒞(tpt)\n\nReturn the set of states that are not in the source or target.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.𝒫-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.𝒫","text":"𝒫(tpt)\n\nReturn tpt.transition_matrix.P.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.𝒫_backwards-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.𝒫_backwards","text":"𝒫_backwards(tpt)\n\nCompute the \"backwards\" transition matrix; return a Matrix.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.𝒫_plus-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.𝒫_plus","text":"𝒫_plus(tpt)\n\nThe forward-reactive analogue of 𝒫.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.𝒮-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.𝒮","text":"𝒮(tpt)\n\nReturn all possible indices implied by the transition matrix, i.e. 1:size(P, 1).\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionPathTheory.𝒮_plus-Tuple{HomogeneousTPTProblem}","page":"API","title":"TransitionPathTheory.𝒮_plus","text":"𝒮_plus(tpt)\n\nThe set of indices i such that\n\nif i is in B_true, then i is in S_plus\nif i is outside B_true, then i is in S_plus if both \nq_minus[i] > 0.0\nsum(P[i, j]*qp[j] for j in sets.S) > 0.0\n\nIntuitively: i is either in B_true, or could have come from A_true and is connected to a state that is reactively connected to B.\n\n\n\n\n\n","category":"method"},{"location":"api/#Types","page":"API","title":"Types","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [TransitionPathTheory]\nOrder = [:type]","category":"page"},{"location":"api/#TransitionPathTheory.Connectivity","page":"API","title":"TransitionPathTheory.Connectivity","text":"abstract type Connectivity\n\nSupertype for all transition matrix connectivity types.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionPathTheory.Disconnected","page":"API","title":"TransitionPathTheory.Disconnected","text":"struct Disconnected\n\nA Connectivity such that the directed graph implied by the transition matrix is disconnected.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionPathTheory.HomogeneousTPTProblem","page":"API","title":"TransitionPathTheory.HomogeneousTPTProblem","text":"struct HomogeneousTPTProblem\n\nA TPTProblem where the transition matrix is homogeneous, i.e. time-independent.\n\nFields\n\ntransition_matrix: A TransitionMatrix.\nsource: A vector of indices defining the TPT source, aka A.\ntarget: A vector of indices defining the TPT target, aka B.\n\nConstructor\n\nHomogeneousTPTProblem(P, source, target; avoid = Int64[])\n\nBuild a HomogeneousTPTProblem where each argument of the constructor maps to the corresponding field.\n\nOptional Arguments \n\navoid: Indices provided here will be appended to both A and B and will hence be avoided by transition path theory.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionPathTheory.NonStochastic","page":"API","title":"TransitionPathTheory.NonStochastic","text":"struct NonStochastic\n\nA Stochasticity such that the row sums of the transition matrix are incomparable to 1.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionPathTheory.Stochastic","page":"API","title":"TransitionPathTheory.Stochastic","text":"struct Stochastic\n\nA Stochasticity such that the row sums of the transition matrix are all equal to 1.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionPathTheory.Stochasticity","page":"API","title":"TransitionPathTheory.Stochasticity","text":"abstract type Stochasticity\n\nSupertype for all transition matrix stochasticity types.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionPathTheory.StronglyConnected","page":"API","title":"TransitionPathTheory.StronglyConnected","text":"struct StronglyConnected\n\nA Connectivity such that the directed graph implied by the transition matrix is strongly connected.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionPathTheory.SubStochastic","page":"API","title":"TransitionPathTheory.SubStochastic","text":"struct SubStochastic\n\nA Stochasticity such that the row sums of the transition matrix are all strictly less than 1.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionPathTheory.SuperStochastic","page":"API","title":"TransitionPathTheory.SuperStochastic","text":"struct SuperStochastic\n\nA Stochasticity such that the row sums of the transition matrix are all strictly greater than 1.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionPathTheory.TPTProblem","page":"API","title":"TransitionPathTheory.TPTProblem","text":"abstract type TPTProblem\n\nSupertype for all TPT problems.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionPathTheory.TransitionMatrix","page":"API","title":"TransitionPathTheory.TransitionMatrix","text":"struct TransitionMatrix{S, C}\n\nA transition probability matrix with Stochasticity S and Connectivity C.\n\nFields\n\nP: The transition probability Matrix itself.\nstochasticity: S\nconnectivity: C\n\nConstructors\n\nTransitionMatrix(P_size; n_zeros = 0, normalize = true, seed = 1234)\n\nConstruct a randomly generated transition matrix.\n\nArguments \n\nP_size: The number of rows (== columns) of the matrix.\n\nOptional Arguments \n\nn_zeros: This many zeros will be placed in the matrix at random. If this results in a row of the matrix having all zeros, a 1 will be placed randomly in that row.\nnormalize: Whether to normalize the row sums of the resulting matrix.\nseed: A seed for reproducible randomness.\n\n\n\nTransitionMatrix(P)\n\nConstruct a TransitionMatrix from a matrix P. The properties of P are inferred automatically.\n\nArguments \n\nP: The transition matrix.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionPathTheory.WeaklyConnected","page":"API","title":"TransitionPathTheory.WeaklyConnected","text":"struct WeaklyConnected\n\nA Connectivity such that the directed graph implied by the transition matrix is weakly connected.\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"advanced/#Advanced-Usage","page":"Advanced Usage","title":"Advanced Usage","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"Pages = [\"advanced.md\"]\nDepth = 5","category":"page"},{"location":"advanced/#High-level-overview","page":"Advanced Usage","title":"High level overview","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"Use of this package proceeds as follows.","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"Define a TransitionMatrix, and provide source, target and avoid indices as desired.\nSelect a TPTProblem.\nCompute the relevant statistics.","category":"page"},{"location":"advanced/#Defining-transition-matrices","page":"Advanced Usage","title":"Defining transition matrices","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"TransitionMatrix","category":"page"},{"location":"advanced/#TransitionPathTheory.TransitionMatrix-advanced","page":"Advanced Usage","title":"TransitionPathTheory.TransitionMatrix","text":"struct TransitionMatrix{S, C}\n\nA transition probability matrix with Stochasticity S and Connectivity C.\n\nFields\n\nP: The transition probability Matrix itself.\nstochasticity: S\nconnectivity: C\n\nConstructors\n\nTransitionMatrix(P_size; n_zeros = 0, normalize = true, seed = 1234)\n\nConstruct a randomly generated transition matrix.\n\nArguments \n\nP_size: The number of rows (== columns) of the matrix.\n\nOptional Arguments \n\nn_zeros: This many zeros will be placed in the matrix at random. If this results in a row of the matrix having all zeros, a 1 will be placed randomly in that row.\nnormalize: Whether to normalize the row sums of the resulting matrix.\nseed: A seed for reproducible randomness.\n\n\n\nTransitionMatrix(P)\n\nConstruct a TransitionMatrix from a matrix P. The properties of P are inferred automatically.\n\nArguments \n\nP: The transition matrix.\n\n\n\n\n\n","category":"type"},{"location":"advanced/#Stochasticity","page":"Advanced Usage","title":"Stochasticity","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"Stochasticity\nStochastic\nSuperStochastic\nNonStochastic","category":"page"},{"location":"advanced/#TransitionPathTheory.Stochasticity-advanced","page":"Advanced Usage","title":"TransitionPathTheory.Stochasticity","text":"abstract type Stochasticity\n\nSupertype for all transition matrix stochasticity types.\n\n\n\n\n\n","category":"type"},{"location":"advanced/#TransitionPathTheory.Stochastic-advanced","page":"Advanced Usage","title":"TransitionPathTheory.Stochastic","text":"struct Stochastic\n\nA Stochasticity such that the row sums of the transition matrix are all equal to 1.\n\n\n\n\n\n","category":"type"},{"location":"advanced/#TransitionPathTheory.SuperStochastic-advanced","page":"Advanced Usage","title":"TransitionPathTheory.SuperStochastic","text":"struct SuperStochastic\n\nA Stochasticity such that the row sums of the transition matrix are all strictly greater than 1.\n\n\n\n\n\n","category":"type"},{"location":"advanced/#TransitionPathTheory.NonStochastic-advanced","page":"Advanced Usage","title":"TransitionPathTheory.NonStochastic","text":"struct NonStochastic\n\nA Stochasticity such that the row sums of the transition matrix are incomparable to 1.\n\n\n\n\n\n","category":"type"},{"location":"advanced/#Connectivity","page":"Advanced Usage","title":"Connectivity","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"Connectivity\nStronglyConnected\nWeaklyConnected\nDisconnected","category":"page"},{"location":"advanced/#TransitionPathTheory.Connectivity-advanced","page":"Advanced Usage","title":"TransitionPathTheory.Connectivity","text":"abstract type Connectivity\n\nSupertype for all transition matrix connectivity types.\n\n\n\n\n\n","category":"type"},{"location":"advanced/#TransitionPathTheory.StronglyConnected-advanced","page":"Advanced Usage","title":"TransitionPathTheory.StronglyConnected","text":"struct StronglyConnected\n\nA Connectivity such that the directed graph implied by the transition matrix is strongly connected.\n\n\n\n\n\n","category":"type"},{"location":"advanced/#TransitionPathTheory.WeaklyConnected-advanced","page":"Advanced Usage","title":"TransitionPathTheory.WeaklyConnected","text":"struct WeaklyConnected\n\nA Connectivity such that the directed graph implied by the transition matrix is weakly connected.\n\n\n\n\n\n","category":"type"},{"location":"advanced/#TransitionPathTheory.Disconnected-advanced","page":"Advanced Usage","title":"TransitionPathTheory.Disconnected","text":"struct Disconnected\n\nA Connectivity such that the directed graph implied by the transition matrix is disconnected.\n\n\n\n\n\n","category":"type"},{"location":"advanced/#Defining-a-TPTProblem","page":"Advanced Usage","title":"Defining a TPTProblem","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"TPTProblem\nHomogeneousTPTProblem","category":"page"},{"location":"advanced/#TransitionPathTheory.TPTProblem-advanced","page":"Advanced Usage","title":"TransitionPathTheory.TPTProblem","text":"abstract type TPTProblem\n\nSupertype for all TPT problems.\n\n\n\n\n\n","category":"type"},{"location":"advanced/#TransitionPathTheory.HomogeneousTPTProblem-advanced","page":"Advanced Usage","title":"TransitionPathTheory.HomogeneousTPTProblem","text":"struct HomogeneousTPTProblem\n\nA TPTProblem where the transition matrix is homogeneous, i.e. time-independent.\n\nFields\n\ntransition_matrix: A TransitionMatrix.\nsource: A vector of indices defining the TPT source, aka A.\ntarget: A vector of indices defining the TPT target, aka B.\n\nConstructor\n\nHomogeneousTPTProblem(P, source, target; avoid = Int64[])\n\nBuild a HomogeneousTPTProblem where each argument of the constructor maps to the corresponding field.\n\nOptional Arguments \n\navoid: Indices provided here will be appended to both A and B and will hence be avoided by transition path theory.\n\n\n\n\n\n","category":"type"},{"location":"advanced/#Basic-sets-in-a-TPTProblem","page":"Advanced Usage","title":"Basic sets in a TPTProblem","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"𝒫\n𝒜\nℬ\n𝒮\nΩ\n𝒜_true\nℬ_true\n𝒞","category":"page"},{"location":"advanced/#TransitionPathTheory.𝒫-advanced","page":"Advanced Usage","title":"TransitionPathTheory.𝒫","text":"𝒫(tpt)\n\nReturn tpt.transition_matrix.P.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.𝒜-advanced","page":"Advanced Usage","title":"TransitionPathTheory.𝒜","text":"𝒜(tpt)\n\nReturn tpt.source.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.ℬ-advanced","page":"Advanced Usage","title":"TransitionPathTheory.ℬ","text":"ℬ(tpt)\n\nReturn tpt.target.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.𝒮-advanced","page":"Advanced Usage","title":"TransitionPathTheory.𝒮","text":"𝒮(tpt)\n\nReturn all possible indices implied by the transition matrix, i.e. 1:size(P, 1).\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.Ω-advanced","page":"Advanced Usage","title":"TransitionPathTheory.Ω","text":"Ω(tpt)\n\nReturn the set of avoided states.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.𝒜_true-advanced","page":"Advanced Usage","title":"TransitionPathTheory.𝒜_true","text":"𝒜_true(tpt)\n\nReturn the set of states that are members of the source, but are not avoided.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.ℬ_true-advanced","page":"Advanced Usage","title":"TransitionPathTheory.ℬ_true","text":"ℬ_true(tpt)\n\nReturn the set of states that are members of the target, but are not avoided.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.𝒞-advanced","page":"Advanced Usage","title":"TransitionPathTheory.𝒞","text":"𝒞(tpt)\n\nReturn the set of states that are not in the source or target.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TPT-Statistics","page":"Advanced Usage","title":"TPT Statistics","text":"","category":"section"},{"location":"advanced/#High-level","page":"Advanced Usage","title":"High level","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"stationary_statistics\nnonstationary_statistics","category":"page"},{"location":"advanced/#TransitionPathTheory.stationary_statistics-advanced","page":"Advanced Usage","title":"TransitionPathTheory.stationary_statistics","text":"stationary_statistics(tpt)\n\nCompute and return the following statistics in a NamedTuple:\n\nstationary_distribution\nforward_committor\nbackward_committor\nreactive_density\nreactive_current\nforward_current\nreactive_rate\nreactive_time\nremaining_time\nhitting_location_distribution\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.nonstationary_statistics-advanced","page":"Advanced Usage","title":"TransitionPathTheory.nonstationary_statistics","text":"nonstationary_statistics(tpt, horizon)\n\nArguments\n\ntpt: The TPTProblem.\nhorizon: The time step at which to cut off the calculation. Note that the horizon value does NOT enforce that trajectories leaving A hit B by that time.\n\nCompute and return the following statistics in a NamedTuple:\n\ndensity\nreactive_density\nt_cdf\nt_cdf_AB\n\n\n\n\n\n","category":"function"},{"location":"advanced/#Low-level","page":"Advanced Usage","title":"Low level","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"stationary_distribution\n𝒫_backwards\nforward_committor\nbackward_committor\n𝒮_plus\n𝒫_plus\nremaining_time\nhitting_location_distribution","category":"page"},{"location":"advanced/#TransitionPathTheory.stationary_distribution-advanced","page":"Advanced Usage","title":"TransitionPathTheory.stationary_distribution","text":"stationary_distribution(tpt)\n\nCompute the stationary distribution of tpt.transition_matrix.P.\n\nErrors when P is not strongly connected.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.𝒫_backwards-advanced","page":"Advanced Usage","title":"TransitionPathTheory.𝒫_backwards","text":"𝒫_backwards(tpt)\n\nCompute the \"backwards\" transition matrix; return a Matrix.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.forward_committor-advanced","page":"Advanced Usage","title":"TransitionPathTheory.forward_committor","text":"forward_committor(tpt)\n\nCompute the forward committor.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.backward_committor-advanced","page":"Advanced Usage","title":"TransitionPathTheory.backward_committor","text":"backward_committor(tpt)\n\nCompute the backward committor.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.𝒮_plus-advanced","page":"Advanced Usage","title":"TransitionPathTheory.𝒮_plus","text":"𝒮_plus(tpt)\n\nThe set of indices i such that\n\nif i is in B_true, then i is in S_plus\nif i is outside B_true, then i is in S_plus if both \nq_minus[i] > 0.0\nsum(P[i, j]*qp[j] for j in sets.S) > 0.0\n\nIntuitively: i is either in B_true, or could have come from A_true and is connected to a state that is reactively connected to B.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.𝒫_plus-advanced","page":"Advanced Usage","title":"TransitionPathTheory.𝒫_plus","text":"𝒫_plus(tpt)\n\nThe forward-reactive analogue of 𝒫.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.remaining_time-advanced","page":"Advanced Usage","title":"TransitionPathTheory.remaining_time","text":"remaining_time(tpt)\n\nCompute the remaining time, aka the lead time.\n\n\n\n\n\n","category":"function"},{"location":"advanced/#TransitionPathTheory.hitting_location_distribution-advanced","page":"Advanced Usage","title":"TransitionPathTheory.hitting_location_distribution","text":"hitting_location_distribution(tpt)\n\nCompute the distribution of target-hitting locations. \n\nReturn a matrix R such that R[i, j] is the probability that - starting in state i - the first visit to B occurs in state j. Note that  R[i, j] == 0 for all j ∉ B.\n\n\n\n\n\n","category":"function"},{"location":"#TransitionPathTheory.jl","page":"Home","title":"TransitionPathTheory.jl","text":"","category":"section"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is a package for computing transition path theory[1] (TPT) statistics a Markov chain on a discrete set mathbbS.  All that is required to apply TPT is ","category":"page"},{"location":"","page":"Home","title":"Home","text":"A strongly connected, transition probability matrix P.\nA \"source\" set mathbbA subset mathbbS and a \"target\" set mathbbB subset mathbbS. The set mathbbA cap mathbbB is automatically avoided by TPT.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The main function of TPT is to compute statistics of \"reactive\" trajectories, namely trajectories which travel directly from the source to the target with no intermediate visit to either.","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Automatic validation and fast computation with simple interface.\nSupports discrete, homogeneous (time-independent) Markov chains.\nStationary and non-stationary TPT statistics.","category":"page"},{"location":"#Documentation","page":"Home","title":"Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is in the Julia General Registry. In the Julia REPL, run the following code and follow the prompts:","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg\nPkg.add(\"TransitionPathTheory\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Access the functionality of the package in your code by including the following line:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TransitionPathTheory","category":"page"},{"location":"#Quickstart","page":"Home","title":"Quickstart","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We will use a randomly generated stochastic transition matrix with 10 states and source mathbbA = 1 2 3  and target mathbbB = 3 4 5.","category":"page"},{"location":"","page":"Home","title":"Home","text":"P = TransitionMatrix(10)  \nA = [1, 2, 3]                           \nB = [3, 4, 5]                           ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Then, we set up HomogeneousTPTProblem and use stationary_statistics and nonstationary_statistics to compute all the relevant statistics in a NamedTuple.","category":"page"},{"location":"","page":"Home","title":"Home","text":"tpt = HomogeneousTPTProblem(P, A, B)\nstats = stationary_statistics(tpt)\nstats_ns = nonstationary_statistics(tpt, 100)  # 100 time step horizon","category":"page"},{"location":"#Citation","page":"Home","title":"Citation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"note: Note\nPlease use the following citation if you use this package in your research.@article{bonner2023improving,\ntitle={Improving the stability of temporal statistics in transition path theory with sparse data},\nauthor={Bonner, Gage and Beron-Vera, FJ and Olascoaga, MJ},\njournal={Chaos: An Interdisciplinary Journal of Nonlinear Science},\nvolume={33},\nnumber={6},\nyear={2023},\npublisher={AIP Publishing}\n}","category":"page"},{"location":"","page":"Home","title":"Home","text":"Initial development of this package was supported by the National Science Foundation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"[1]: Vanden-Eijnden, Eric. \"Transition path theory.\" Computer Simulations in Condensed Matter Systems: From Materials to Chemical Biology Volume 1. Springer, Berlin, Heidelberg, 2006. 453-493.","category":"page"}]
}
