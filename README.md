# Hit and Run

A Julia package based on the paper, "Hit and Run as a Unifying Device." by P. Diaconis & H. Andersen, (2007), *Journal de la Société Française de Statistique*, 148(4):5-28, [[PDF](http://statweb.stanford.edu/~cgates/PERSI/papers/hitandrun062207.pdf)] and its references.

This package is primarily a personal experiment, an attempt to gauge my own interest (and ability) in subjects which the first author, above, regards as revolutionary [[PDF](http://statweb.stanford.edu/~cgates/PERSI/papers/MCMCRev.pdf)]. Since Julia is potentially revolutionary as well<sup>1</sup>, the experiment seems apt.

I've chosen the Hit and Run paper mainly because it includes and references accessible examples such as generating random contingency tables, i.e., arrays with non-negative integer entries and fixed row and column sums<sup>2</sup>.

The paper's promise of unification, "an attempt to make sense of the zoo of Monte Carlo techniques," also appeals. The appeal is not only the prospect of avoiding a steep learning curve, but because generality facilitates variation and experiment, hence (possibly) progress.

1. [Julia: A Fresh Approach to Numerical Computing](http://arxiv.org/abs/1411.1607) Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah: "[Julia] allows programmers to write clear, high-level, generic and abstract code that closely resembles mathematical formulas ... yet produces fast, low-level machine code that has traditionally only been generated by static languages.”

2. I had originally expected to use models as specified in [Mamba](https://mambajl.readthedocs.org/en/latest/), an excellent, purely Julia platform for Monte Carlo Markov Chains. However, that specification, which is very natural for Bayesian networks, is somewhat unnatural for the examples to be considered here. Hence, I expect to use custom code to produce instances of the `Chains` type, to which Mamba's convergence diagnostics and other features apply.
