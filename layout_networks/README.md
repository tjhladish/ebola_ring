## Observation Process Model Outline

We represent the process of observing the network as a respondent-driven sample.

Each node in the network will name (or not) it's neighbors in the network with probability $p$ ($1-p$).  The identified neighbors will either be findable (on non-findable) name, with probability $q$ ($1-q$).

Explaining by way of example, if $L0$ has 3 contacts--A, B, and C--when surveyed, it might identify 2 (forgetting, say, C), and of those two findably identify one (say, A, with B' [confusingly named] being unfound).

Repeating this process for the found nodes, we can build up the observed network and compare it to the simulated reference network.  In the explanatory example, we would think the network is $(L0-A, L0-B')$ if we did no further sampling.  However, sampling beyond $L0$ requires the model represent how we might deal with supplemental information.  To extend our example, assume that in fact L0, A, B, and C form a completely connected network (a clique).  When we survey A (the only L0 contact we can survey initially), A identifies L0, B (not B'), but still not C.  At this point we think the network is $(L0-A, L0-B', A-B)$.  We then sample B, who identifies A' (not A), L0, and C.  These ambiguous identifications (A', and L0 [since L0 identified B', not B]) require a model choice to resolve.  We could:

 - assume that the prime people are real people that we cannot find - i.e., the example observed network is $(L0-A, L0-B', A-B, B-A', B-L0', B-C)$.
 - or that prime people are not real, when the real target correctly identifies the source(s) of a prime edge - i.e., the example observed network is $(L0-A, L0-B, A-B, B-C)$.

We will use the latter mechanism.  Lastly, we come to C, which identifies A, B, and L0'.  Because L0 did not identify C, we assume L0' is a new (unfindable) node, but we otherwise connect up A (despite A not naming C).