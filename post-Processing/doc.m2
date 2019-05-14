doc ///
   Key 
      SchurVeronese
   Headline 
      syzygy data for Veronese embeddings of $\mathbb P^1$ and $\mathbb P^2$
   Description
    Text
      This package includes all the syzygy data for a big computation...
///

doc ///
   Key 
    makeBettiTally
    (makeBettiTally,HashTable)
   Headline
    converts a hash table reprenting a Betti table to a Betti tally
   Usage
    makeBettiTally(H)
   Inputs
    H: HashTable
   Outputs
    : BettiTally
   Description
    Text
      Given a hash table $H$ whose keys are pairs of integers $(p,q)$ 
      such that $H#(p,q)$ is $\beta_{p,p+q}(M)$ for some module $M$
      this function outputs the BettiTally corresponding to the given
      hash table.     
    Example
      H = totalBetti(3,2,0);
      makeBettiTally H      
///

doc ///
   Key 
    totalBetti
    (totalBetti,ZZ,ZZ,ZZ)
   Headline
    a hash table containing the graded Betti numbers of a Veronese embedding
   Usage
    totalBetti(d,n,b)
   Inputs
    d: ZZ
    n: ZZ
    b: ZZ 
   Outputs
    : HashTable
   Description
    Text
      This is a hash table for the total numbers of $\mathcal{O}_{\mathbb{P}^n}(b)$
      under the embedding by $\mathcal{O}_{\mathbb{P}^n}(d)$. The keys of the hash
      table $H$ are pairs $(p,q)$ where $H#(p,q)$ gives the rank of 
      $K_{p,q}(\mathbb{P}^n, d;b)$. This equals the Betti number
      $\beta_{p,p+q}(\mathbb{P}^n, d;b)$.Some tables are incomplete and 
      we mark unknown entries with infinity.
      
      Note that totalBetti only differs from totalBetti2 in that the output is 
      a HashTable instead of a BettiTally. One can convert the output of
      totalBetti into a BettiTally via the makeBettiTally function.

      
      In example below we generate a hash table showing the total graded Betti numbers
      of $\mathbb{P}^{2}$ embedded by $\mathcal{O}_{\mathbb{P}^{2}}(3)$. 
    Example
      B = totalBetti(3,2,0)
    Text  
      If we wish to view these graded Betti numbers in the usual fashion, we can use
      makeBettiTally to convert the hash table above to a BettiTally.
    Example
      makeBettiTally B  
///

doc ///
   Key 
    totalBetti2
    (totalBetti2,ZZ,ZZ,ZZ)
   Headline
    a Betti tally containing the graded Betti numbers of a Veronese embedding
   Usage
    totalBetti2(d,n,b)
   Inputs
    d: ZZ
    n: ZZ
    b: ZZ 
   Outputs
    : BettiTally
   Description
    Text
      This function outputs a BettiTally for the total graded Betti numbers 
      of $\mathcal{O}_{\mathbb{P}^n}(b)$ under the embedding by 
      $\mathcal{O}_{\mathbb{P}^n}(d)$. Some tables are incomplete and we mark
      unknown entries with infinity.
      
      Note that totalBetti2 only differs from totalBetti in that the output is 
      a BettiTaly instead of a HashTable. 
      
      In example below we generate a hash table showing the total graded Betti numbers
      of $\mathbb{P}^{2}$ embedded by $\mathcal{O}_{\mathbb{P}^{2}}(3)$. 
    Example
      totalBetti2(3,2,0)
            
///


doc ///
   Key 
    schurBetti
    (schurBetti,ZZ,ZZ,ZZ)
   Headline
    a hash table for Schur module decomposition of Veronese Betti tables
   Usage
    schurBetti(d,n,b)
   Inputs
    d: ZZ
    n: ZZ
    b: ZZ 
   Outputs
    : HashTable
   Description
    Text
      This function ouptputs a hash table with the Schur functor decompositions 
      of the syzygies of $\mathcal{O}_{\mathbb{P}^n}(b)$ under the embedding by
      $\mathcal{O}_{\mathbb{P}^n}(d)$. The keys of the hash
      table $H$ are pairs $(p,q)$ where $H#(p,q)$ gives the the Schur functor 
      decomposition of $K_{p,q}(\mathbb{P}^n, d;b)$. We record the 
      Schur functor decomosition as a list of tuples $({a,b,c},b)$ where
      ${a,b,c}$ specifes the weight of the Schur functor and $m$ the multiplicty
      with which that particular Schur functor appears in the decomposition
      of $K_{p,q}(\mathbb{P}^n, d;b)$. 
      
      Some tables are incomplete and we mark unknown entries with 
      ({0,0,0},infinity).
      
    Example
      schurBetti(3,2,0)
///

doc ///
   Key 
    multiBetti
    (multiBetti,ZZ,ZZ,ZZ)
   Headline
    a hash table containing the multigraded Betti numbers of a Veronese embedding 
   Usage
    multiBetti(d,n,b)
   Inputs
    d: ZZ
    n: ZZ
    b: ZZ 
   Outputs
    : HashTable
   Description
    Text
      This function ouptputs a hash table $H$ containing the multigraded Betti numbers
      for $\mathcal{O}_{\mathbb{P}^n}(b)$ on $\mathbb{P}^{n}$ under the embedding by
      d'uple Veronese embedding given by $\mathcal{O}_{\mathbb{P}^n}(d)$. The keys 
      of the outputed hash table $H$ are pairs $(p,q)$ where $H#(p,q)$ gives the the
      multigraded Betti decomposition of $K_{p,q}(\mathbb{P}^n, d;b)$. We record the 
      multigraded Betti numbers via a multigraded Hilbert series. See Section 1.1 of 
      [BEGY] for a more precise descripion of the multigraded Betti numbers of a 
      Veronese embedding. 
      
    Example
      multiBetti(3,2,0)
///

doc ///
   Key 
    dominantWeightsBetti
    (dominantWeightsBetti,ZZ,ZZ,ZZ)
   Headline
    a hash table containing the dominant Schur functors of a Veronese embedding 
   Usage
    dominantWeightsBetti(d,n,b)
   Inputs
    d: ZZ
    n: ZZ
    b: ZZ 
   Outputs
    : HashTable
   Description
    Text
      This function ouptputs a hash table $H$ whose keys are pairs $(p,q)$ such that
      the corresponding value $H#(p,q)$ is a list of the dominant weights appearing
      in the Schur functor decomposition of $K_{p,q}(\mathbb{P}^n, d;b)$. The Schur
      functors are recorded via their weights. See Section 1.3 of [BEGY].
///

doc ///
   Key 
    lexWeightsBetti
    (lexWeightsBetti,ZZ,ZZ,ZZ)
   Headline
    a hash table containing the Lex leading weight Schur functors of a Veronese embedding 
   Usage
    lexWeightsBetti(d,n,b)
   Inputs
    d: ZZ
    n: ZZ
    b: ZZ 
   Outputs
    : HashTable
   Description
    Text
      This function ouptputs a hash table $H$ whose keys are pairs $(p,q)$ such that
      the corresponding value $H#(p,q)$ is a list of the lex-leading weight Schur
      functors appearing in the Schur functor decomposition of 
      $K_{p,q}(\mathbb{P}^n, d;b)$. The Schur functors are recorded via their weights.
      See Section 1.3 of [BEGY].
      
    Example
      lexWeightsBetti(3,2,0)
   
///

doc ///
   Key 
    numDistinctRepsBetti
    (numDistinctRepsBetti,ZZ,ZZ,ZZ)
   Headline
    a hash table containing the number of distinct Schur functors of a Veronese embedding 
   Usage
    numDistinctRepsBetti(d,n,b)
   Inputs
    d: ZZ
    n: ZZ
    b: ZZ 
   Outputs
    : HashTable
   Description
    Text
      This function ouptputs a hash table $H$ whose keys are pairs $(p,q)$ such that
      the corresponding value $H#(p,q)$ is the number of distinct Schur functors
      appearing in the Schur functor decomposition of $K_{p,q}(\mathbb{P}^n, d;b)$. 
      
    Example
      numDistinctRepsBetti(3,2,0)
      
///

doc ///
   Key 
    numRepsBetti
    (numRepsBetti,ZZ,ZZ,ZZ)
   Headline
    a hash table containing the number of Schur functors of a Veronese embedding 
   Usage
    numRepsBetti(d,n,b)
   Inputs
    d: ZZ
    n: ZZ
    b: ZZ 
   Outputs
    : HashTable
   Description
    Text
      This function ouptputs a hash table $H$ whose keys are pairs $(p,q)$ such that
      the corresponding value $H#(p,q)$ is the number of Schur functors appearing
      in the Schur functor decomposition of $K_{p,q}(\mathbb{P}^n, d;b)$ counted
      with multiplicity.  
      
    Example
      numRepsBetti(3,2,0)
      
///

doc ///
   Key 
    errorBetti
    (errorBetti,ZZ,ZZ,ZZ)
   Headline
    a hash table containing the errors encountered when computing sysygies of a Veronese embedding 
   Usage
    errorBetti(d,n,b)
   Inputs
    d: ZZ
    n: ZZ
    b: ZZ 
   Outputs
    : HashTable
   Description
    Text
      As the methods used to compute the multigraded Betti numbers that are at the 
      heart of the SchurVeronese package are numerical in natur there is room for
      error. Thus, we have implemented methods in post processing to catch and correct
      errors. This function ouptputs a hash table $H$ whose keys are pairs $(p,q)$ such that
      the corresponding value $H#(p,q)$ is a multigraded Hilbert series recording the
      $K_{p,q}(\mathbb{P}^n, d;b)$. See [Sec 5.2, BEGY] for a discussion on error 
      processing.
      
      Note that we did not enoucter errors in any cases when $d\leq 5$. However, there
      were errors for some cases when $d>5$. 
    Example
      errorBetti(3,2,0)
      errorBetti(6,2,3)
///

doc ///
   Key 
    bsCoeffs
    (bsCoeffs,ZZ,ZZ,ZZ)
   Headline
    a list of the Boij-Soederberg coefficents of a Veronese embedding 
   Usage
    bsCoeffs(d,n,b)
   Inputs
    d: ZZ
    n: ZZ
    b: ZZ 
   Outputs
    : List
   Description
    Text
      This function returns a list of the Boij-Soederberg coefficents for the 
      decomposition of the Betti table of $\mathcal{O}_{\mathbb{P}^n}(b)$ on 
      $\mathbb{P}^{n}$ under the embedding by d'uple Veronese embedding 
      given by $\mathcal{O}_{\mathbb{P}^n}(d)$. See Section 6.3 of [BEGY].
      
    Example
      bsCoeffs(3,2,0)
///