[![GitHub stars](https://img.shields.io/github/stars/sgtlaugh/algovault.svg?style=flat&label=star)](https://github.com/sgtlaugh/algovault/)
[![license](https://img.shields.io/github/license/sgtlaugh/algovault.svg?style=flat-square)](https://github.com/sgtlaugh/algovault/blob/master/LICENSE)

## Algovault
### A collection of algorithms, data structures and templates for competitive programming

<li>Useful in online competitions like <a href="https://codeforces.com">CodeForces</a>, <a href="https://codingcompetitions.withgoogle.com/codejam">Google Code Jam</a></li>
<li>Simple to use as a black-box without compromising performance</li>
<li>Example usage and sufficient documentation</li>

<br>
Codes are mostly written in C++. majority should work with C++11 and some might require C++14 or higher. Some algorithms are written in Python. For Python, use Python 3. Most of them should be compatible with PyPy as well.
</br>

<br>
Implementations are usually stress-tested and cross-checked against various problems. Nonetheless, they are not guranteed to be flawless and work in all cases.
</br>

<br>
For bugs, refactoring and improvements, feel free to file an issue or a pull request as contributions are always welcome.
</br>

## Index

C++ files compile on their own: copy one into a solution and use it as a black box. Every file ends with a `main()` of self-tests, and has a brute-force stress test under [`stress_tests/`](stress_tests) at the same relative path.

### Data Structures
- [Coordinate Compression](code_library/data_structures/coordinate_compression.cpp) - in-place compression, optionally order preserving
- [Disjoint Set Union](code_library/data_structures/disjoint_set.cpp) - path compression and union by size
- [Disjoint Sparse Table](code_library/data_structures/disjoint_sparse_table.cpp) - O(1) static range queries for any associative operation
- [Fenwick Tree](code_library/data_structures/fenwick_tree.cpp) - point/range update with point/range query
- [Fenwick Tree 2D](code_library/data_structures/fenwick_tree_2D.cpp) - point/range update with point/range query on a grid
- [Fenwick Tree 2D, Implicit](code_library/data_structures/fenwick_tree_2D_implicit.cpp) - huge grids, columns as on-demand segment trees
- [Fenwick Tree 2D, Sparse](code_library/data_structures/fenwick_tree_2D_sparse.cpp) - grids up to 1e9 x 1e9 with sparse updates, hashed
- [Fenwick Tree 3D](code_library/data_structures/fenwick_tree_3D.cpp) - point/range update with point/range query in 3D
- [Hash Map](code_library/data_structures/hashmap.cpp) - anti-hack hash map with expected O(1) operations
- [Li Chao Tree](code_library/data_structures/li_chao_tree.cpp) - min/max of lines and segments at a point
- [Link-Cut Tree](code_library/data_structures/link_cut_tree.cpp) - dynamic forest with link, cut, re-rooting, LCA and path aggregates
- [Merge Sort Tree](code_library/data_structures/merge_sort_tree.cpp) - count below a threshold and k-th smallest in a range
- [Mo's Algorithm](code_library/data_structures/mo.cpp) - offline range queries on arrays and tree paths
- [Monotonic Stack](code_library/data_structures/monotonic_stack.cpp) - nearest smaller element on each side, largest histogram rectangle
- [Ordered Set](code_library/data_structures/ordered_set.cpp) - GNU policy-based set and multiset with order statistics
- [Persistent Segment Tree](code_library/data_structures/persistent_segment_tree.cpp) - versioned point add/set, range query with any associative merge, k-th on version differences, range k-th smallest, persistent array
- [Rope](code_library/data_structures/rope.cpp) - persistent text, every edit makes a new readable version
- [Segment Tree](code_library/data_structures/segment_tree.cpp) - lazy propagation, customizable merge and update
- [Sparse Table](code_library/data_structures/sparse_table.cpp) - O(1) static range min, plus an O(n) build variant
- [Square Root Decomposition](code_library/data_structures/sqrt_decomposition.cpp) - range add and count of values below x in a range
- [Treap](code_library/data_structures/treap.cpp) - ordered multiset and implicit-key sequence with reversals
- [Trie](code_library/data_structures/trie.cpp) - prefix tree with pass-through and end counts

### Trees
- [Centroid Decomposition](code_library/trees/centroid_decomposition.cpp) - iterative centroid tree with a per-centroid visitor over branch distances, path counting example
- [DSU on Tree](code_library/trees/dsu_on_tree.cpp) - small-to-large answers for every subtree
- [Heavy Light Decomposition](code_library/trees/hld.cpp) - tree paths and subtrees as O(log n) ranges for any range structure
- [Lowest Common Ancestor](code_library/trees/lca.cpp) - binary lifting with k-th ancestor, and an O(n) / O(1) variant

### Graphs
- [2-SAT](code_library/graphs/2SAT_kosaraju.cpp) - satisfiability and an assignment via Kosaraju
- [2-SAT, Lexicographic](code_library/graphs/2SAT_lexicographic.cpp) - lexicographically smallest satisfying assignment
- [Articulation Points and Biconnected Components](code_library/graphs/articulation_points.cpp) - cut vertices, blocks (edge ids) and the block-cut tree, parallel edges allowed
- [Bridges](code_library/graphs/bridge.cpp) - bridges and the bridge tree, parallel edges allowed
- [Directed MST](code_library/graphs/directed_mst.cpp) - minimum arborescence, O(m log n)
- [Euler Path](code_library/graphs/euler_path.cpp) - Hierholzer, directed and undirected multigraphs, iterative, O(n + m)
- [General Graph Matching](code_library/graphs/graph_matching.cpp) - maximum matching size via the Tutte matrix rank, Edmonds blossom for the matched pairs in O(n^3)
- [Gomory-Hu Tree](code_library/graphs/gomory_hu.cpp) - all-pairs min cut and min cut partitions of an undirected graph from n - 1 max flows (Gusfield)
- [Hopcroft Karp](code_library/graphs/hopcroft_karp.cpp) - maximum bipartite matching, O(m sqrt(n))
- [Hungarian Algorithm](code_library/graphs/hungarian_algorithm.cpp) - minimum cost assignment, O(min(n, m)^2 max(n, m))
- [Johnson's Algorithm](code_library/graphs/johnsons_algorithm.cpp) - all-pairs shortest paths with negative edges
- [Maximum Flow](code_library/graphs/maxflow.cpp) - Dinic, with node capacities and a dense variant
- [Min Cost Max Flow, Dijkstra](code_library/graphs/mcmf_dijkstra.cpp) - successive shortest paths with potentials, the default
- [Min Cost Max Flow, SPFA](code_library/graphs/mcmf_spfa.cpp) - successive shortest paths with SPFA
- [Minimum Path Cover](code_library/graphs/minimum_path_cover.cpp) - disjoint and shared path covers of a DAG, maximum antichain
- [Minimum Spanning Tree](code_library/graphs/minimum_spanning_tree.cpp) - Kruskal and Boruvka (struct with add_edge, template for implicit graphs) plus dense O(n^2) Prim on a weight matrix, spanning forest of an undirected graph
- [Offline Dynamic Connectivity](code_library/graphs/offline_dynamic_connectivity.cpp) - connectivity and component count under edge insertions and deletions, O((m + q) log(m + q) log n)
- [Steiner Tree](code_library/graphs/steiner_tree.cpp) - minimum tree connecting k terminals, O(3^k n + 2^k m log n)
- [Strongly Connected Components](code_library/graphs/scc.cpp) - Kosaraju, components in topological order

### Strings
- [2D Pattern Matcher](code_library/strings/2D_pattern_matcher.cpp) - every occurrence of a 2D pattern in a 2D text by hashing
- [Aho-Corasick](code_library/strings/aho_corasick.cpp) - multi-pattern matching automaton
- [Aho-Corasick, Dynamic](code_library/strings/dynamic_aho_corasick.cpp) - multi-pattern matching with online pattern insertion
- [Bit-String LCS](code_library/strings/bit_string_lcs.cpp) - longest common subsequence with bitsets, O(nm / 64)
- [Dynamic String Hash](code_library/strings/dynamic_string_hash.cpp) - substring hashes under range assignment
- [Hashing](code_library/strings/hashing.cpp) - forward and reverse polynomial hash of any segment
- [Hunt-Szymanski](code_library/strings/hunt_szymanski.cpp) - LCS in O((r + n) log n) for r matching pairs
- [KMP](code_library/strings/kmp.cpp) - failure function and pattern search
- [Manacher](code_library/strings/manacher.cpp) - longest palindrome at every center
- [Minimum Rotation](code_library/strings/minimum_rotation.cpp) - start of the lexicographically smallest rotation
- [Palindromic Tree](code_library/strings/palindromic_tree.cpp) - eertree, one node per distinct palindrome
- [Suffix Array](code_library/strings/suffix_array.cpp) - DC3 in O(n) with the LCP array
- [Suffix Automaton](code_library/strings/suffix_automaton.cpp) - online DFA of all substrings with distinct counts, occurrence counts and longest common substring
- [Z Algorithm](code_library/strings/z_algorithm.cpp) - longest common prefix of every suffix with the string

### Number Theory
- [All Divisors](code_library/number_theory/all_divisors.cpp) - divisor lists of every number up to a limit
- [Chinese Remainder Theorem](code_library/number_theory/chinese_remainder_theorem.cpp) - solution of a system of congruences
- [Digits of Factorial](code_library/number_theory/digits_of_factorial.cpp) - number of digits of n! in any base
- [Discrete Logarithm](code_library/number_theory/discrete_logarithm.cpp) - smallest x with a^x = b mod m, any m
- [Divisors](code_library/number_theory/divisors.cpp) - sorted divisors of one number, or from its prime factors
- [Fast Fibonacci](code_library/number_theory/fast_fibonacci.cpp) - n-th Fibonacci number by fast doubling
- [Fast Prime Counting](code_library/number_theory/fast_prime_counting.cpp) - pi(n) in about O(n^(2/3)), Meissel-Lehmer
- [Fast Prime Sums](code_library/number_theory/fast_prime_sums.cpp) - sum of primes up to n, Meissel-Lehmer
- [Fast Sieve](code_library/number_theory/fast_sieve.cpp) - optimized sieve of Eratosthenes up to 2^31 - 1
- [Floor Sum](code_library/number_theory/floor_sum.cpp) - sum of floor((a * i + b) / m), min of (a * i + b) mod m, count and first x with a * x mod m in [l, r]
- [Integer Root](code_library/number_theory/integer_root.cpp) - exact floor of square, cube and k-th roots of 64-bit integers
- [Linear Sieve](code_library/number_theory/linear_sieve.cpp) - smallest prime factors, phi, mu and any multiplicative function
- [Maximum Divisors](code_library/number_theory/maximum_divisors.cpp) - number with the most divisors up to a limit
- [Miller Rabin](code_library/number_theory/miller_rabin.cpp) - deterministic primality test for 64-bit integers
- [Pisano Period](code_library/number_theory/pisano_period.cpp) - period of Fibonacci numbers modulo n
- [Pollard Rho](code_library/number_theory/pollard_rho.cpp) - factorization of 64-bit integers
- [Power Tower](code_library/number_theory/power_tower.cpp) - a1^(a2^(...^an)) modulo m
- [Segmented Sieve](code_library/number_theory/segmented_sieve.cpp) - primes in a window [L, R] far from zero

### Combinatorics
- [Bernoulli Numbers](code_library/combinatorics/bernoulli_numbers.cpp) - Bernoulli numbers modulo a prime
- [Binomial Coefficients](code_library/combinatorics/binomial_coefficients.cpp) - n choose k modulo any m, Lucas and prime powers
- [Combinatorics](code_library/combinatorics/combinatorics.cpp) - extended gcd, modular inverse, diophantine equations, nCr and nPr
- [Eulerian Numbers](code_library/combinatorics/eulerian_numbers.cpp) - permutations with k ascents, modulo any m
- [Faulhaber's Formula](code_library/combinatorics/faulhaber_formula.cpp) - sum of k-th powers of 1..n
- [Gambler's Ruin](code_library/combinatorics/gamblers_ruin.cpp) - probability that the first player goes broke in the gambler's ruin game
- [Josephus Problem](code_library/combinatorics/josephus_problem.cpp) - survivor of the Josephus elimination
- [Permutation Index](code_library/combinatorics/permutation_index.cpp) - lexicographic rank of a permutation, a perfect hash
- [Permutation Rank](code_library/combinatorics/permutation_rank.cpp) - rank and unrank of permutations
- [Stirling Numbers](code_library/combinatorics/stirling_numbers.cpp) - a whole row of either kind modulo m using FFT

### Algebra
- [Big Integer](code_library/algebra/bignum.cpp) - arbitrary precision signed integers
- [FFT](code_library/algebra/fft.cpp) - polynomial multiplication, exact modular and 64-bit products
- [Fraction](code_library/algebra/fraction.cpp) - exact rational arithmetic, always reduced
- [Lagrange Interpolation](code_library/algebra/lagrange_interpolation.cpp) - value of a polynomial from consecutive samples
- [Linear Recurrence](code_library/algebra/linear_recurrence.cpp) - find a recurrence from terms and compute its n-th term
- [NTT](code_library/algebra/ntt.cpp) - polynomial multiplication modulo an NTT prime, any modulus, or exact in 64 bits
- [Polynomial](code_library/algebra/polynomial.cpp) - NTT power series: inverse, division, log, exp, sqrt, pow, multipoint evaluation, interpolation and Taylor shift
- [Walsh Hadamard Transform](code_library/algebra/walsh_hadamard.cpp) - xor, or and and convolutions

### Linear Algebra
- [Determinant](code_library/linear_algebra/determinant.cpp) - integer determinant modulo any m, or exact
- [Freivalds' Algorithm](code_library/linear_algebra/freivalds_algorithm.cpp) - randomized check of a matrix product
- [Gauss, Band Matrix](code_library/linear_algebra/gauss_band_matrix.cpp) - linear systems whose equations touch nearby variables only
- [Gauss, Bitset](code_library/linear_algebra/gauss_bitset.cpp) - linear systems over GF(2)
- [Gauss Jordan](code_library/linear_algebra/gauss_jordan.cpp) - linear systems over the reals
- [Gauss, Prime Modulus](code_library/linear_algebra/gauss_prime_mod.cpp) - linear systems modulo a prime
- [Matrix](code_library/linear_algebra/matrix.cpp) - modular multiplication and exponentiation
- [Maximum XOR Subset](code_library/linear_algebra/max_xor_subset.cpp) - largest xor of any subset via a linear basis
- [Simplex](code_library/linear_algebra/simplex.cpp) - linear programming
- [Thomas Algorithm](code_library/linear_algebra/thomas_algorithm.cpp) - tridiagonal linear systems in O(n)

### Geometry
- [Circle Geometry](code_library/geometry/circle.cpp) - circle-line/circle-circle intersections, tangents, circumcircle, minimum enclosing circle, circle-polygon, lens and union areas, max points covered by a radius-r circle
- [Closest Pair of Points](code_library/geometry/closest_pair.cpp) - O(n log n) sweep with exact integer squared distances
- [Convex Hull](code_library/geometry/convex_hull.cpp) - monotone chain, O(n log n)
- [Geometry](code_library/geometry/geometry.cpp) - exact integer predicates, polygon queries and angular sort, floating point line and segment intersection, projection, reflection, half-plane polygon cut and centroid
- [Half-plane Intersection](code_library/geometry/halfplane_intersection.cpp) - convex region of half-planes clipped to a box, sort and deque, O(n log n)
- [Pick's Theorem](code_library/geometry/picks_theorem.cpp) - lattice points inside and on a polygon
- [Polygon Union](code_library/geometry/polygon_union.cpp) - area of the union of simple polygons, convex or not, O(N^2 log N)
- [Segment Intersection Sweep](code_library/geometry/segment_intersection_sweep.cpp) - Shamos-Hoey, finds a pair among n closed segments that share a point, O(n log n), exact __int128 predicates

### Dynamic Programming
- [Aliens Trick (WQS Binary Search, Lagrangian Relaxation)](code_library/dp/aliens_trick.cpp) - exactly-k optimum of a convex cost via binary search on the per-item penalty, tie-safe
- [Blocks](code_library/dp/blocks_dp.cpp) - interval DP for UVA 10559 Blocks
- [CKY](code_library/dp/cky.cpp) - context-free grammar membership in Chomsky normal form
- [Concave 1D1D DP](code_library/dp/concave_1d1d_dp.cpp) - dp[x] = min dp[i] + w(i, x) in O(n log n) under the quadrangle inequality
- [Divide and Conquer DP Optimization](code_library/dp/divide_conquer_dp.cpp)
- [Knuth Optimization](code_library/dp/knuth_optimization.cpp) - O(n^2) interval DP dp[i][j] = min dp[i][k] + dp[k][j] + cost(i, j) for monotone quadrangle-inequality costs (stone merging, optimal BST)
- [Longest Common Increasing Subsequence](code_library/dp/lcis.cpp) - LCIS in O(nm)
- [Longest Increasing Subsequence](code_library/dp/lis.cpp) - LIS and LDS lengths with any comparator
- [Maximum Square](code_library/dp/maximum_square.cpp) - largest filled square and diamond ending at every cell

### Miscellaneous
- [15 Puzzle Solver](code_library/misc/15_puzzle_solver.cpp) - solvability check and IDA* solver
- [Assembly](code_library/misc/assembly.cpp) - inline x86 popcount, leading zeros, bit scan, gcd and square root
- [Bit Twiddling](code_library/misc/bit_twiddling.cpp) - bit manipulation tricks
- [Comb Sort](code_library/misc/combsort.cpp) - comb sort on any random access range
- [Dancing Links](code_library/misc/dancing_links.cpp) - exact cover with Algorithm X
- [Fast I/O](code_library/misc/fast_io.cpp) - buffered input and output with fread and fwrite
- [Gray Codes](code_library/misc/gray_codes.cpp) - Gray code and its inverse
- [Knight Distance](code_library/misc/knight_distance.cpp) - fewest knight moves on an infinite board in O(1)
- [N Queens](code_library/misc/n_queen.cpp) - number of N queens solutions
- [Next Palindrome](code_library/misc/next_palindrome.cpp) - smallest palindromic number above a given one
- [Radix Sort](code_library/misc/radix_sort.cpp) - LSD radix sort of 32-bit integers
- [Radix Sort 64](code_library/misc/radix_sort_64.cpp) - LSD radix sort of 64-bit integers
- [Simulated Annealing](code_library/misc/simulated_annealing.cpp) - generic minimizer with a geometric time or iteration cooling schedule

### Hacking
- [Anti Double Hash](code_library/hacking/anti_double_hash.cpp) - two strings colliding under a double polynomial hash
- [Anti Polynomial Hash](code_library/hacking/anti_polyonmial_hash.cpp) - two strings colliding under a polynomial hash

### Python
- [Alpha Beta Pruning](code_library/python/alpha_beta_pruning.py) - minimax game search with alpha-beta pruning
- [Berlekamp Massey](code_library/python/berlekamp_massey.py) - shortest linear recurrence of a sequence modulo a prime
- [Derangements](code_library/python/derangements.py) - generalized derangement counts
- [Dijkstra](code_library/python/dijkstra.py) - single-source shortest paths
- [Discrete Log](code_library/python/discrete_log.py) - smallest x with a^x = b mod m, any m
- [Fast Fibonacci](code_library/python/fast_fibonacci.py) - n-th Fibonacci number by fast doubling
- [Fast I/O](code_library/python/fast_io.py) - token scanner over stdin
- [General Graph Matching](code_library/python/graph_matching.py) - maximum matching via the Tutte matrix rank
- [Lagrange Interpolation](code_library/python/lagrange_polynomial_interpolation.py) - value of a polynomial from any samples
- [Miller Rabin](code_library/python/miller_rabin.py) - deterministic primality test for 64-bit integers
- [N Queens](code_library/python/n_queen.py) - number of N queens solutions
- [Sieve](code_library/python/sieve.py) - sieve of Eratosthenes
- [String Hash](code_library/python/string_hash.py) - polynomial substring hashes
- [Z Algorithm](code_library/python/z_algorithm.py) - longest common prefix of every suffix with the string

## The Zen Of Contributing
Inspired from [The Zen Of Python](https://www.python.org/dev/peps/pep-0020/#id2)

```python
Beautiful is better than ugly

Simple is way better than complex

Consistency matters

Flat is preferred over nested

Typing is better than incomprehensible macros

Because readability counts

But not as much as speed

Some documentation is better than no documentation

No documentation is better than extensive documentation

But not as important as ease of reusing as a black box

Four spaces are better than tabs

Tabs are better than no spaces

If the implementation is hard to explain, it's a bad idea

If the implementation is easy to explain, it may be a good idea

Structs are one honking great idea, let's do more of those!
```

## Future Work
  <ol>
  <li>This is still a work in progress so I'll port more code from my template over the time</li>
  <li>Refactor and simplify old implementations</li>
  <li>Add practice problems</li>
  </ol>

## License
The project is licensed under the [MIT License](https://github.com/sgtlaugh/algovault/blob/master/LICENSE)
