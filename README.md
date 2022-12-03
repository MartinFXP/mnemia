# mnemia

## Files:

src/mnemia.jl contains all functions for inference and a function to simulate data

## Install/usage:

mnemia is not available as a full julia package, but you can activate it within
its directory (cd("path")) in the pkg environment (] key).

```
activate .
```

## Usage:

Create data from four networks with five genes, ten effectors per gene and 100 samples.

```julia
import mnemia as mn
using Plots

n = 5
K = 4
m = 10
s = 100

truth = mn.sim(K, n, m, s, 0.5)
```

We extract the data object and the perturbation map (which gene has been perturbed
in which sample) from the simulation object. We also add a bit of noise to the data
and convert them to log ratios.

```julia
D = truth[4]
rho = truth[3]
R = mn.addNoise(D)
```

Finally, we use the mnem algorithm to compute the mixture and can compare to the
ground truth. We use ten different starts each initialised with a random set of
networks.

```julia
best = mn.mnem(R, rho, K)
bestll = findmax(best[2])[1]
lll = [best[2][1,:]]
lls = bestll
for i in 1:10
    result = mn.mnem(R, rho, K)
    ll = findmax(result[2])[1]
    lll = push!(lll, result[2][1,:])
    lls = hcat(lls,ll)
    if ll > bestll
        best = result
        bestll = ll
    end
    print(i)
end
best[1]
truth[1]
p1 = mn.plotMix(best)
p2 = plot(lll)
plot(p1, p2, layout = (1,2))
```

Depending on noise and complexity of the ground truth even more runs are advised. If
the number of components K is unknown, result for several different Ks need to be
compared. However, the log likelihood cannot be directly used for this. Some form of
Information Criterion like BIC or AIC (not implemented) is necessary to select the best K.

The method infers a mixture of dags by default, but can be restricted to trees, e.g., if
we are interested in mutational history (Pirkl et al., 2022).

```julia
trees = mn.sim(K, n, m, s, 0.5, true) # with edge probability 0.5 for denser networks
Rt = mn.addNoise(trees[4])
result = mn.mnem(Rt, trees[3], K, 100, true) # the local optimisation of the tree mixture is stopped after 100 (default) iterations
mn.plotMix(result)
```

Genotypic information is often provided in binary form with mutation observed (1) or not (0). In
this case you can use the following function to convert the data to log ratios. However, you must
know error probabilities for false positive and false negative observations. The
optimisation of both parameters is not integrated in mnemia.

```julia
Rdt = mn.transformData(trees[4], 0.01, 0.1) # the data from the simulation is actually provided in binary
result = mn.mnem(Rdt, trees[3], K, 100, true)
mn.plotMix(result)
```

DISCLAIMER:
-----------

This program is fully functional and computes a mixture of directed graphs based on
the M&NEM algorithm (Pirkl & Beerenwinkel. 2018).

However, this program is written in julia as a programming exercise and should
only be used for testing purposes (feedback is appreciated). For real data analysis
I recommend the more thoroughly tested R/Bioconductor package mnem (https://github.com/cbg-ethz/mnem).

## References:

Pirkl, M., Beerenwinkel, N.; Single cell network analysis with a mixture
of Nested Effects Models, Bioinformatics, Volume 34, Issue 17, 1 September
2018,
Pages i964-i971, https://doi.org/10.1093/bioinformatics/bty602.

Pirkl, M., Büch, J., Devaux, C., Böhm, M., Sönnerborg, A., Incardona, F., Abecasis, A., Vandamme, A.-M., Zazzi, M., Kaiser, R., Lengauer, T. and the EuResist Network study group. 2022. Analysis of mutational history of multi-drug resistant genotypes with a mutagenetic tree model. Journal of Medical Virology. _accepted_