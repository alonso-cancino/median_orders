# Median orders:

A median order of a tournament $T=(V,E)$ is an ordering of the vertices $\sigma:V\to [n]$ that maximices the quantity:

$$e(\sigma,T) = |\{ e\in E(T) : \sigma(e_{1}) \leq \sigma(e_{2}) \}|$$

This is equivalent to the following optimization problem:

$$\max_{P\in \mathbb{N}^{(n,n)}} \sum\limits_{i =1}^{n} \sum_{j=1}^n (P^{t}AP)_{i,j}$$

subject to $(\vec{1})^{t}P=(\vec{1})^{t}$ and $P\hspace{0.7mm}\vec{1} = \vec{1}$, this yields an optimal permutation matrix $P$, which in turn induces an ordering given by $P (0,1,...,n-1)^t$.

The function ```find_median_order_rt(n,verbose)``` initialices a random tournament $T$ and its adjacency matrix $A_T$ in $n$ vertices, and then computes a linearization of the problem described above.
