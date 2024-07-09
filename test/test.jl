using Test, TransitionPathTheory

A = [1, 2, 3]
B = [3, 4, 5]
P = TransitionMatrix(10)

tpt = HomogeneousTPTProblem(P, A, B)
st_stats = stationary_statistics(tpt)
ns_stats = nonstationary_statistics(tpt, 10)