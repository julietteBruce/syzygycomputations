doc ///
   Key 
      HirzebruchSyzygies
   Headline 
      syzygy data for Veronese embeddings of $\mathbb P^1$ and $\mathbb P^2$
   Description
    Text
      The authors of package used a combination of high throughput
      high perfomance computing and sparse numerical linear algebra to compute 
      the syzygies of $\mathbb{P}^{2}$ under the $d$-uple Veronese embedding
      for a number of values of $d$. See the paper ``Conjectures and 
      Computations about Veronese Syzygies'' by Bruce, Erman, Goldstein and Yang, which
      we refer to as [BEGY] throughout the documentation for this package.
      In addition, much of the data generated from
      these computations (graded Betti numbers, multigraded Betti numbers, 
      Schur functor decompositions, etc.) is currently available online via syzygydata.com.
      The goal of this package is to make this data more accessible and easy to 
      use by providing a way to access it via Macaulay2.

      Most functions have been implemented with three parameters $(d,n,b)$ where
      the goal is to compute the syzygies of the pushforward of
      the line bundle $\mathcal O_{\mathbb{P}^2}(b)$ under the $d$-uple
      embedding.  However, we have produced data for $n=1,2$ and for $b$ between $0$ and $d$ and
      for a limited range of values of $d$.  Other inputs will produce an error message.  Our hope is
      that as we (or others) are able to compute new data, we will be able to update the package.

      
      One of the main functions is totalBettiTally which produces the standard
      Betti tables.  Other main functions refine the data in the Betti table
      by providing the multigraded Betti number or the Schur functor decompositions,
      or by computing statistics related to the Betti table (e.g. the BoijSoederberg
      coefficients) or related to the SchurFunctor decomposition.      

///

