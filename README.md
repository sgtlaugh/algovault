# Algovault

**Algorithms, data structures and templates for competitive programming**

[![CI](https://github.com/sgtlaugh/algovault/actions/workflows/ci.yml/badge.svg)](https://github.com/sgtlaugh/algovault/actions/workflows/ci.yml)
[![license](https://img.shields.io/github/license/sgtlaugh/algovault.svg?style=flat)](https://github.com/sgtlaugh/algovault/blob/master/LICENSE)

A black box you can copy easily - simple, fast, and well-tested. Templates are C++17 for GCC, with a few in Python 3 that mostly run on PyPy too, for contests like [Codeforces](https://codeforces.com), [AtCoder](https://atcoder.jp) and ICPC.

## How To Use

Every C++ file compiles on its own. Copy it into your solution, delete its `main()`, which only holds self-tests, and use the struct:

```cpp
DSU dsu(n);
dsu.connect(a, b);
bool same = dsu.find_root(a) == dsu.find_root(b);
```

The header at the top of each file lists its API, complexity and limits.

## Testing

Every template is checked four ways:

| Layer | Where | Checks |
|---|---|---|
| Self-tests | `main()` in each file | Usage examples with known answers |
| Stress tests | [`stress_tests/`](stress_tests), same path | Random inputs against an independent brute force |
| Repo checks | [`run_tests.sh`](.github/scripts/run_tests.sh) | Every file has a stress test and an index entry, consistent style |
| Judge tests | [`judge_tests/`](judge_tests) | Official [Library Checker](https://judge.yosupo.jp) and [Aizu](https://onlinejudge.u-aizu.ac.jp) test data at full size, within the real time limit |

Self-tests and stress tests build with GCC 14 under AddressSanitizer and UndefinedBehaviorSanitizer with an 8 MB stack, warnings as errors. Judge tests build with `-O2` like a judge.

- Every push runs the self, stress and judge tests of the files it changes, plus the repo checks
- Every night all self and stress tests run, and the stress tests again at 20x the iterations with a fresh seed
- Every week all judge tests run from scratch

A ✔ in the index marks a template that passes a judge test; it links to one of them. The commands to run everything locally are in [CONTRIBUTING.md](CONTRIBUTING.md#checks).

## Index

[Data Structures](#data-structures) | [Trees](#trees) | [Graphs](#graphs) | [Strings](#strings) | [Number Theory](#number-theory) | [Combinatorics](#combinatorics) | [Algebra](#algebra) | [Linear Algebra](#linear-algebra) | [Geometry](#geometry) | [Dynamic Programming](#dynamic-programming) | [Miscellaneous](#miscellaneous) | [Hacking](#hacking) | [Python](#python)

### Data Structures

[Practice problems](code_library/data_structures/README.md)

- [Coordinate Compression](code_library/data_structures/coordinate_compression.cpp) - in-place compression, optionally order preserving
- [Disjoint Set Union](code_library/data_structures/disjoint_set.cpp) - union by size with path compression, rollback, weighted and persistent variants [✔](judge_tests/data_structures/persistent_unionfind.cpp)
- [Distinct Subarray Aggregates](code_library/data_structures/distinct_subarray_aggregates.cpp) - distinct gcd/or/and of subarrays ending at each index, `O(log A)`
- [Fenwick Tree](code_library/data_structures/fenwick_tree.cpp) - point/range update with point/range query, `O(log n)` lower bound on prefix sums [✔](judge_tests/data_structures/aoj_dsl_2_e.cpp)
  - [2D](code_library/data_structures/fenwick_tree_2D.cpp) - the same on a grid [✔](judge_tests/data_structures/aoj_dsl_5_b.cpp)
  - [2D, Implicit](code_library/data_structures/fenwick_tree_2D_implicit.cpp) - huge grids, columns as on-demand segment trees
  - [2D, Sparse](code_library/data_structures/fenwick_tree_2D_sparse.cpp) - grids up to 1e9 x 1e9 with sparse updates, hashed [✔](judge_tests/data_structures/point_add_rectangle_sum.cpp)
  - [3D](code_library/data_structures/fenwick_tree_3D.cpp) - the same in 3D
- [Hash Map](code_library/data_structures/hashmap.cpp) - anti-hack hash map with expected `O(1)` operations [✔](judge_tests/data_structures/associative_array.cpp)
- [Interval Set](code_library/data_structures/interval_set.cpp) - disjoint half-open intervals with merging add, splitting remove and coverage queries [✔](judge_tests/data_structures/aoj_dsl_4_a.cpp)
- [Leftist Heap, Persistent](code_library/data_structures/leftist_heap.cpp) - meldable priority queue with `O(log n)` push/pop/meld, persistent [✔](judge_tests/data_structures/k_shortest_walk.cpp)
- [Li Chao Tree](code_library/data_structures/li_chao_tree.cpp) - min/max of lines and segments at a point, plus a persistent variant [✔](judge_tests/data_structures/line_add_get_min.cpp)
- [Link-Cut Tree](code_library/data_structures/link_cut_tree.cpp) - dynamic forest with link, cut, re-rooting, LCA and path aggregates [✔](judge_tests/data_structures/dynamic_tree_vertex_add_path_sum.cpp)
- [Merge Sort Tree and Wavelet](code_library/data_structures/merge_sort_tree.cpp) - count below a threshold, range frequency, range k-th smallest [✔](judge_tests/data_structures/range_kth_smallest.cpp)
- [Mo's Algorithm](code_library/data_structures/mo.cpp) - offline range queries on arrays, with updates, with rollback, and on tree paths [✔](judge_tests/data_structures/point_set_range_frequency.cpp)
- [Monotonic Stack](code_library/data_structures/monotonic_stack.cpp) - nearest smaller element on each side, largest histogram rectangle [✔](judge_tests/data_structures/aoj_dpl_3_c.cpp)
- [Ordered Set](code_library/data_structures/ordered_set.cpp) - GNU policy-based set and multiset with order statistics [✔](judge_tests/data_structures/ordered_set.cpp)
- [Permutation Tree](code_library/data_structures/permutation_tree.cpp) - `O(n log n)` decomposition into common intervals, counts all of them [✔](judge_tests/data_structures/common_interval_decomposition_tree.cpp)
- [Rope](code_library/data_structures/rope.cpp) - persistent text, every edit makes a new readable version
- [Segment Tree](code_library/data_structures/segment_tree.cpp) - lazy propagation, customizable merge and update, max_right / min_left descent [✔](judge_tests/data_structures/point_set_range_composite.cpp)
  - [Beats](code_library/data_structures/segment_tree_beats.cpp) - range chmin, chmax, add and assign with range sum, min and max queries [✔](judge_tests/data_structures/range_chmin_chmax_add_range_sum.cpp)
  - [Merging](code_library/data_structures/segment_tree_merge.cpp) - multisets as dynamic segment trees with merge, split by key or rank, k-th
  - [Persistent](code_library/data_structures/persistent_segment_tree.cpp) - versioned updates, range queries, k-th smallest and k-th on version diffs [✔](judge_tests/data_structures/persistent_queue.cpp)
  - [XOR](code_library/data_structures/xor_segment_tree.cpp) - range sum of `a[p ^ x]` over `[l, r]` for any x, with point add
- [Sliding Window Aggregation](code_library/data_structures/sliding_window_aggregation.cpp) - fold under any associative operation and a monotonic min/max queue [✔](judge_tests/data_structures/aoj_dsl_3_d.cpp)
- [Sparse Table](code_library/data_structures/sparse_table.cpp) - `O(1)` static range min, plus an `O(n)` build variant [✔](judge_tests/data_structures/staticrmq.cpp)
  - [2D](code_library/data_structures/sparse_table_2D.cpp) - `O(1)` static rectangle min or max on a grid
  - [Disjoint](code_library/data_structures/disjoint_sparse_table.cpp) - `O(1)` static range queries for any associative operation [✔](judge_tests/data_structures/static_range_sum.cpp)
- [Square Root Decomposition](code_library/data_structures/sqrt_decomposition.cpp) - range add and count of values below x in a range
- [Treap](code_library/data_structures/treap.cpp) - ordered multiset and implicit-key sequence with reversals [✔](judge_tests/data_structures/double_ended_priority_queue.cpp)
- [Trie](code_library/data_structures/trie.cpp) - prefix tree with pass-through and end counts, and a binary trie for xor queries [✔](judge_tests/data_structures/aoj_alds1_4_c.cpp)

### Trees

[Practice problems](code_library/trees/README.md)

- [Centroid Decomposition](code_library/trees/centroid_decomposition.cpp) - iterative centroid tree with a visitor over distances to each centroid [✔](judge_tests/trees/frequency_table_of_tree_distance.cpp)
- [DSU on Tree](code_library/trees/dsu_on_tree.cpp) - small-to-large answers for every subtree
- [Dynamic Tree Diameter](code_library/trees/dynamic_diameter.cpp) - weighted tree diameter under edge weight updates in `O(log n)` [✔](judge_tests/trees/aoj_grl_5_a.cpp)
- [Heavy Light Decomposition](code_library/trees/hld.cpp) - tree paths and subtrees as `O(log n)` ranges for any range structure [✔](judge_tests/trees/aoj_grl_5_d.cpp)
- [Long-Path Decomposition](code_library/trees/long_path_decomposition.cpp) - `O(n)` depth-indexed subtree DP and `O(1)` k-th ancestor via ladders [✔](judge_tests/trees/jump_on_tree_long_path_decomposition.cpp)
- [Lowest Common Ancestor](code_library/trees/lca.cpp) - k-th ancestor, path jump and intersection, plus an `O(n)` / `O(1)` LCA [✔](judge_tests/trees/jump_on_tree.cpp)
- [Prufer Code](code_library/trees/prufer_code.cpp) - `O(n)` encode/decode between labeled trees and Prufer sequences
- [Rerooting DP](code_library/trees/rerooting.cpp) - a tree DP answered for every node as the root in `O(n)` [✔](judge_tests/trees/aoj_grl_5_b.cpp)
- [Tree Isomorphism](code_library/trees/tree_isomorphism.cpp) - AHU canonical ids for rooted and unrooted trees, exact via a map dictionary [✔](judge_tests/trees/rooted_tree_isomorphism_classification.cpp)
- [Virtual Tree](code_library/trees/virtual_tree.cpp) - compress a tree to k vertices and their LCAs in `O(k log k)`

### Graphs

[Practice problems](code_library/graphs/README.md)

- [2-SAT](code_library/graphs/2SAT_kosaraju.cpp) - satisfiability and an assignment via Kosaraju, with at-most-one constraints in `O(k)` clauses [✔](judge_tests/graphs/two_sat.cpp)
- [2-SAT, Lexicographic](code_library/graphs/2SAT_lexicographic.cpp) - lexicographically smallest satisfying assignment [✔](judge_tests/graphs/two_sat_lexicographic.cpp)
- [3-Edge-Connected Components](code_library/graphs/three_edge_connected_components.cpp) - Tsin, iterative, parallel edges and self loops allowed, `O(n + m)` [✔](judge_tests/graphs/three_edge_connected_components.cpp)
- [Articulation Points and Blocks](code_library/graphs/articulation_points.cpp) - cut vertices, blocks and the block-cut tree, parallel edges allowed [✔](judge_tests/graphs/aoj_grl_3_a.cpp)
- [Bellman Ford](code_library/graphs/bellman_ford.cpp) - shortest paths with negative edges, -inf marking, negative cycle retrieval, `O(n m)` [✔](judge_tests/graphs/aoj_grl_1_b.cpp)
- [Bridges](code_library/graphs/bridge.cpp) - bridges and the bridge tree, parallel edges allowed [✔](judge_tests/graphs/aoj_grl_3_b.cpp)
- [Cactus Graph](code_library/graphs/cactus.cpp) - edge/vertex cactus recognition, cycle listing and cactus tree, `O(n + m)`
- [Chromatic Number](code_library/graphs/chromatic_number.cpp) - minimum vertex coloring by inclusion-exclusion, `O(2^n n)` for `n <= 24` [✔](judge_tests/graphs/chromatic_number.cpp)
- [Complement Graph BFS](code_library/graphs/complement_graph_bfs.cpp) - BFS and components in the complement of a sparse graph, `O(n + m)` [✔](judge_tests/graphs/connected_components_of_complement_graph.cpp)
- [Dijkstra](code_library/graphs/dijkstra.cpp) - shortest paths with parents, heap `O((n + m) log m)`, dense `O(n^2)`, 0-1 BFS [✔](judge_tests/graphs/aoj_alds1_11_c.cpp)
- [Directed MST](code_library/graphs/directed_mst.cpp) - minimum arborescence with parents, negative weights, `O(m log m)` [✔](judge_tests/graphs/directedmst.cpp)
- [Dominator Tree](code_library/graphs/dominator_tree.cpp) - Lengauer-Tarjan immediate dominators from a root, `O((n + m) log n)` [✔](judge_tests/graphs/dominatortree.cpp)
- [Edge Coloring](code_library/graphs/edge_coloring.cpp) - bipartite with D colors (Konig), simple graphs with D + 1 colors (Vizing) [✔](judge_tests/graphs/bipartite_edge_coloring.cpp)
- [Euler Path](code_library/graphs/euler_path.cpp) - Hierholzer, directed and undirected multigraphs, iterative, `O(n + m)` [✔](judge_tests/graphs/eulerian_trail_directed.cpp)
- [Flow with Lower Bounds](code_library/graphs/lower_bound_flow.cpp) - feasible circulation, feasible / max / min s-t flow with demands on Dinic [✔](judge_tests/graphs/aoj_grl_6_a_lower_bound.cpp)
- [Floyd Warshall](code_library/graphs/floyd_warshall.cpp) - all pairs shortest paths with paths, negative edges and negative cycle marking [✔](judge_tests/graphs/aoj_grl_1_c.cpp)
- [General Graph Matching](code_library/graphs/graph_matching.cpp) - maximum matching by Tutte rank, pairs by blossom, `O(n^3)` [✔](judge_tests/graphs/general_matching.cpp)
- [Global Minimum Cut](code_library/graphs/global_min_cut.cpp) - Stoer-Wagner, minimum cut weight and one side, `O(n^3)`
- [Gomory-Hu Tree](code_library/graphs/gomory_hu.cpp) - all-pairs min cuts of an undirected graph from n - 1 max flows
- [Hopcroft Karp](code_library/graphs/hopcroft_karp.cpp) - maximum bipartite matching, `O(m sqrt(n))` [✔](judge_tests/graphs/bipartitematching.cpp)
- [Hungarian Algorithm](code_library/graphs/hungarian_algorithm.cpp) - minimum cost assignment, `O(min(n, m)^2 max(n, m))` [✔](judge_tests/graphs/assignment.cpp)
- [Johnson's Algorithm](code_library/graphs/johnsons_algorithm.cpp) - all-pairs shortest paths with negative edges [✔](judge_tests/graphs/aoj_grl_1_c_johnson.cpp)
- [K Shortest Walks](code_library/graphs/k_shortest_walks.cpp) - Eppstein with a persistent leftist heap, `O((n + m) log m + k log k)` [✔](judge_tests/graphs/k_shortest_walk_k_shortest_walks.cpp)
- [Kruskal Reconstruction Tree](code_library/graphs/kruskal_reconstruction_tree.cpp) - vertices reachable from v using edges <= w, as a range, `O(log n)`
- [Manhattan MST](code_library/graphs/manhattan_mst.cpp) - minimum spanning tree of points under L1 distance, `O(n log n)` [✔](judge_tests/graphs/manhattanmst.cpp)
- [Matroid Intersection](code_library/graphs/matroid_intersection.cpp) - largest common independent set, graphic/partition/xor, `O(n r^2)` queries [✔](judge_tests/graphs/aoj_grl_7_a.cpp)
- [Maximum Clique](code_library/graphs/max_clique.cpp) - bitset branch and bound to ~150 vertices, independent set, clique enumeration [✔](judge_tests/graphs/maximum_independent_set.cpp)
- [Maximum Flow](code_library/graphs/maxflow.cpp) - Dinic with node capacities and a dense variant, min cut and closure recipes [✔](judge_tests/graphs/aoj_grl_6_a.cpp)
- [Min Cost Circulation](code_library/graphs/min_cost_circulation.cpp) - cost scaling with negative costs and cycles, `O(n^3 log(n C))` [✔](judge_tests/graphs/aoj_grl_6_b_circulation.cpp)
- [Min Cost Max Flow, Dijkstra](code_library/graphs/mcmf_dijkstra.cpp) - successive shortest paths with potentials and the cost slope, the default [✔](judge_tests/graphs/aoj_grl_6_b.cpp)
- [Min Cost Max Flow, SPFA](code_library/graphs/mcmf_spfa.cpp) - successive shortest paths with SPFA [✔](judge_tests/graphs/aoj_grl_6_b_spfa.cpp)
- [Minimum Mean Cycle (Karp)](code_library/graphs/minimum_mean_cycle.cpp) - exact minimum mean cycle as a fraction, with its edges, `O(n (n + m))`
- [Minimum Path Cover](code_library/graphs/minimum_path_cover.cpp) - disjoint and shared path covers of a DAG, maximum antichain [✔](judge_tests/graphs/aoj_2251.cpp)
- [Minimum Spanning Tree](code_library/graphs/minimum_spanning_tree.cpp) - Kruskal, Boruvka for implicit graphs, dense `O(n^2)` Prim [✔](judge_tests/graphs/aoj_alds1_12_a.cpp)
- [Offline Dynamic Connectivity](code_library/graphs/offline_dynamic_connectivity.cpp) - connectivity and components under edge insertions and deletions
- [Online Bridges](code_library/graphs/online_bridges.cpp) - bridges and 2-edge components under insertions, `O((n log n + m) α(n))` [✔](judge_tests/graphs/two_edge_connected_components.cpp)
- [Range Edge Dijkstra](code_library/graphs/range_edge_dijkstra.cpp) - shortest paths with vertex-to-range and range-to-vertex edges
- [Small Cycle Counting](code_library/graphs/small_cycle_counting.cpp) - triangles and 4-cycles of a simple graph in `O(n + m sqrt(m))`
- [Stable Marriage](code_library/graphs/stable_marriage.cpp) - Gale-Shapley, proposer optimal stable matching, `O(n^2)`
- [Steiner Tree](code_library/graphs/steiner_tree.cpp) - minimum tree connecting k terminals, `O(3^k n + 2^k m log n)`
- [Strongly Connected Components](code_library/graphs/scc.cpp) - Kosaraju, components in topological order [✔](judge_tests/graphs/scc.cpp)

### Strings

[Practice problems](code_library/strings/README.md)

- [2D Pattern Matcher](code_library/strings/2D_pattern_matcher.cpp) - every occurrence of a 2D pattern in a 2D text by hashing [✔](judge_tests/strings/aoj_alds1_14_c.cpp)
- [Aho-Corasick](code_library/strings/aho_corasick.cpp) - multi-pattern matching automaton [✔](judge_tests/strings/aho_corasick.cpp)
- [Aho-Corasick, Dynamic](code_library/strings/dynamic_aho_corasick.cpp) - multi-pattern matching with online pattern insertion
- [Bit-String LCS](code_library/strings/bit_string_lcs.cpp) - longest common subsequence with bitsets, `O(n m / 64)` [✔](judge_tests/strings/aoj_alds1_10_c.cpp)
- [De Bruijn Sequence](code_library/strings/de_bruijn_sequence.cpp) - lexicographically smallest `B(k, n)` by Lyndon words, `O(k^n)`
- [Dynamic String Hash](code_library/strings/dynamic_string_hash.cpp) - substring hashes and palindrome checks under range assignment [✔](judge_tests/strings/aoj_alds1_14_b_dynamic_string_hash.cpp)
- [Hashing](code_library/strings/hashing.cpp) - forward and reverse polynomial hash of any segment [✔](judge_tests/strings/enumerate_palindromes_hashing.cpp)
- [Hunt-Szymanski](code_library/strings/hunt_szymanski.cpp) - LCS in `O((r + n) log n)` for r matching pairs [✔](judge_tests/strings/aoj_alds1_10_c_hunt_szymanski.cpp)
- [KMP](code_library/strings/kmp.cpp) - failure function, pattern search and the prefix function automaton [✔](judge_tests/strings/aoj_alds1_14_b.cpp)
- [Lyndon Factorization](code_library/strings/lyndon_factorization.cpp) - Duval factorization into Lyndon words, plus a least rotation [✔](judge_tests/strings/lyndon_factorization.cpp)
- [Manacher](code_library/strings/manacher.cpp) - longest palindrome at every center [✔](judge_tests/strings/enumerate_palindromes.cpp)
- [Minimum Rotation](code_library/strings/minimum_rotation.cpp) - start of the lexicographically smallest rotation
- [Palindromic Tree](code_library/strings/palindromic_tree.cpp) - eertree, one node per distinct palindrome [✔](judge_tests/strings/eertree.cpp)
- [Suffix Array](code_library/strings/suffix_array.cpp) - DC3 in `O(n)` with LCP array, `O(1)` LCP queries and pattern ranges [✔](judge_tests/strings/aoj_alds1_14_d.cpp)
- [Suffix Automaton](code_library/strings/suffix_automaton.cpp) - substring counts (distinct, occurrences) and longest common substring [✔](judge_tests/strings/number_of_substrings.cpp)
- [Tandem Repeats](code_library/strings/tandem_repeats.cpp) - Main-Lorentz, every square as `O(n log n)` compact triples [✔](judge_tests/strings/runenumerate.cpp)
- [Wildcard Matching](code_library/strings/wildcard_matching.cpp) - matching with wildcards in text or pattern via NTT, `O(n log n)` [✔](judge_tests/strings/wildcard_pattern_matching.cpp)
- [Z Algorithm](code_library/strings/z_algorithm.cpp) - longest common prefix of every suffix with the string [✔](judge_tests/strings/zalgorithm.cpp)

### Number Theory

[Practice problems](code_library/number_theory/README.md)

- [All Divisors](code_library/number_theory/all_divisors.cpp) - divisor lists of every number up to a limit
- [Chinese Remainder Theorem](code_library/number_theory/chinese_remainder_theorem.cpp) - systems of congruences with any moduli, and Garner's mixed radix [✔](judge_tests/number_theory/convolution_mod_1000000007_chinese_remainder_theorem.cpp)
- [Continued Fractions](code_library/number_theory/continued_fractions.cpp) - convergents, best rational approximation, Stern-Brocot search [✔](judge_tests/number_theory/rational_approximation.cpp)
- [Digits of Factorial](code_library/number_theory/digits_of_factorial.cpp) - number of digits of n! in any base
- [Discrete Logarithm](code_library/number_theory/discrete_logarithm.cpp) - smallest x with `a^x = b mod m`, any m [✔](judge_tests/number_theory/discrete_logarithm_mod.cpp)
- [Divisors](code_library/number_theory/divisors.cpp) - sorted divisors of one number, or from its prime factors
- [Du's Sieve](code_library/number_theory/du_sieve.cpp) - prefix sums of Euler's phi and Mobius mu up to 1e11 in `O(n^(2/3))` [✔](judge_tests/number_theory/counting_squarefrees.cpp)
- [Fast Fibonacci](code_library/number_theory/fast_fibonacci.cpp) - n-th Fibonacci number by fast doubling [✔](judge_tests/number_theory/aoj_alds1_10_a.cpp)
- [Fast Prime Counting](code_library/number_theory/fast_prime_counting.cpp) - `pi(n)` in about `O(n^(2/3))`, Meissel-Lehmer [✔](judge_tests/number_theory/counting_primes.cpp)
- [Fast Prime Sums](code_library/number_theory/fast_prime_sums.cpp) - sum of primes up to n, Meissel-Lehmer
- [Fast Sieve](code_library/number_theory/fast_sieve.cpp) - optimized sieve of Eratosthenes up to `2^31 - 1` [✔](judge_tests/number_theory/enumerate_primes.cpp)
- [Floor Sum](code_library/number_theory/floor_sum.cpp) - `sum floor((a i + b) / m)`, `min (a i + b) mod m`, first x with `a x mod m` in `[l, r]` [✔](judge_tests/number_theory/min_of_mod_of_linear.cpp)
- [Integer Root](code_library/number_theory/integer_root.cpp) - exact floor of square, cube and k-th roots of 64-bit integers [✔](judge_tests/number_theory/kth_root_integer.cpp)
- [Linear Sieve](code_library/number_theory/linear_sieve.cpp) - smallest prime factors, phi, mu and any multiplicative function
- [Maximum Divisors](code_library/number_theory/maximum_divisors.cpp) - number with the most divisors up to a limit
- [Miller Rabin](code_library/number_theory/miller_rabin.cpp) - deterministic primality test for 64-bit integers [✔](judge_tests/number_theory/primality_test.cpp)
- [Min_25 Sieve](code_library/number_theory/min_25_sieve.cpp) - prefix sums of a multiplicative function in `O(deg n^(3/4) / log n)`, any modulus [✔](judge_tests/number_theory/sum_of_multiplicative_function.cpp)
- [Mobius Function](code_library/number_theory/mobius_function.cpp) - Mobius mu for `1..n` in `O(n)` at 1 byte per entry
- [Modular Roots](code_library/number_theory/modular_roots.cpp) - multiplicative order, primitive root, square root and k-th roots modulo a prime [✔](judge_tests/number_theory/kth_root_mod.cpp)
- [Pisano Period](code_library/number_theory/pisano_period.cpp) - period of Fibonacci numbers modulo n
- [Pollard Rho](code_library/number_theory/pollard_rho.cpp) - factorization of 64-bit integers [✔](judge_tests/number_theory/factorize.cpp)
- [Power Tower](code_library/number_theory/power_tower.cpp) - `a1^(a2^(...^an))` modulo m [✔](judge_tests/number_theory/tetration_mod.cpp)
- [Segmented Sieve](code_library/number_theory/segmented_sieve.cpp) - primes in a window `[L, R]` far from zero [✔](judge_tests/number_theory/enumerate_primes_segmented_sieve.cpp)

### Combinatorics

[Practice problems](code_library/combinatorics/README.md)

- [Bernoulli Numbers](code_library/combinatorics/bernoulli_numbers.cpp) - Bernoulli numbers modulo a prime
- [Binomial Coefficients](code_library/combinatorics/binomial_coefficients.cpp) - n choose k modulo any m, Lucas and prime powers [✔](judge_tests/combinatorics/binomial_coefficient.cpp)
- [Combinatorics](code_library/combinatorics/combinatorics.cpp) - extended gcd, modular inverse, diophantine equations, nCr and nPr [✔](judge_tests/combinatorics/aoj_ntl_1_e.cpp)
- [Eulerian Numbers](code_library/combinatorics/eulerian_numbers.cpp) - permutations with k ascents, modulo any m
- [Faulhaber's Formula](code_library/combinatorics/faulhaber_formula.cpp) - sum of k-th powers of 1..n
- [Gambler's Ruin](code_library/combinatorics/gamblers_ruin.cpp) - probability that the first player goes broke in the gambler's ruin game
- [Josephus Problem](code_library/combinatorics/josephus_problem.cpp) - survivor of the Josephus elimination [✔](judge_tests/combinatorics/aoj_0085.cpp)
- [Partition Numbers](code_library/combinatorics/partition_numbers.cpp) - `p(0..n)` modulo any m in `O(n sqrt(n))` via the pentagonal number theorem [✔](judge_tests/combinatorics/partition_function.cpp)
- [Permutation Cycles](code_library/combinatorics/permutation_cycles.cpp) - cycles, k-th power and k-th root in `O(n)`, parity and minimum swaps [✔](judge_tests/combinatorics/aoj_alds1_6_d.cpp)
- [Permutation Index](code_library/combinatorics/permutation_index.cpp) - lexicographic rank of a permutation, a perfect hash
- [Permutation Rank](code_library/combinatorics/permutation_rank.cpp) - rank and unrank of permutations [✔](judge_tests/combinatorics/aoj_itp2_5_c.cpp)
- [Stirling Numbers](code_library/combinatorics/stirling_numbers.cpp) - a whole row of either kind modulo m using exact NTT [✔](judge_tests/combinatorics/stirling_number_of_the_first_kind.cpp)

### Algebra

[Practice problems](code_library/algebra/README.md)

- [Big Integer](code_library/algebra/bignum.cpp) - arbitrary precision signed integers [✔](judge_tests/algebra/addition_of_big_integers.cpp)
- [FFT](code_library/algebra/fft.cpp) - polynomial multiplication, exact modular and 64-bit products [✔](judge_tests/algebra/convolution_mod_1000000007_fft.cpp)
- [Fraction](code_library/algebra/fraction.cpp) - exact rational arithmetic, always reduced
- [GCD and LCM Convolution](code_library/algebra/gcd_lcm_convolution.cpp) - `c[k] = sum of a[i] b[j]` over `gcd(i, j) = k` or `lcm(i, j) = k` [✔](judge_tests/algebra/gcd_convolution.cpp)
- [Lagrange Interpolation](code_library/algebra/lagrange_interpolation.cpp) - value of a polynomial from consecutive samples [✔](judge_tests/algebra/sum_of_exponential_times_polynomial.cpp)
- [Linear Recurrence](code_library/algebra/linear_recurrence.cpp) - find a recurrence from terms and compute its n-th term [✔](judge_tests/algebra/find_linear_recurrence.cpp)
- [Min-Plus Convolution](code_library/algebra/min_plus_convolution.cpp) - convex-convex in `O(n + m)`, convex-arbitrary via monotone minima [✔](judge_tests/algebra/min_plus_convolution_convex_arbitrary.cpp)
- [NTT](code_library/algebra/ntt.cpp) - polynomial multiplication modulo an NTT prime, any modulus, or exact in 64 bits [✔](judge_tests/algebra/convolution_mod.cpp)
- [Polynomial](code_library/algebra/polynomial.cpp) - NTT power series: inv, div, log, exp, sqrt, pow, multipoint eval, interpolation, Taylor shift [✔](judge_tests/algebra/division_of_polynomials.cpp)
- [Subset Convolution and SOS DP](code_library/algebra/subset_convolution.cpp) - zeta and Mobius transforms, subset convolution in `O(2^n n^2)` [✔](judge_tests/algebra/bitwise_and_convolution.cpp)
- [Walsh Hadamard Transform](code_library/algebra/walsh_hadamard.cpp) - xor, or and and convolutions [✔](judge_tests/algebra/bitwise_and_convolution_fwht.cpp)

### Linear Algebra

[Practice problems](code_library/linear_algebra/README.md)

- [Characteristic Polynomial](code_library/linear_algebra/characteristic_polynomial.cpp) - `det(xI - A)` modulo a prime in `O(n^3)` via Hessenberg reduction [✔](judge_tests/linear_algebra/characteristic_polynomial.cpp)
- [Determinant](code_library/linear_algebra/determinant.cpp) - integer determinant modulo any m or exact, Kirchhoff and Tutte tree counts [✔](judge_tests/linear_algebra/counting_spanning_tree_directed.cpp)
- [Freivalds' Algorithm](code_library/linear_algebra/freivalds_algorithm.cpp) - randomized check of a matrix product
- [Gauss Jordan](code_library/linear_algebra/gauss_jordan.cpp) - linear systems, rank and matrix inverse over the reals
  - [Band Matrix](code_library/linear_algebra/gauss_band_matrix.cpp) - linear systems whose equations touch nearby variables only
  - [Bitset](code_library/linear_algebra/gauss_bitset.cpp) - linear systems, rank and matrix inverse over GF(2) [✔](judge_tests/linear_algebra/inverse_matrix_mod_2.cpp)
  - [Prime Modulus](code_library/linear_algebra/gauss_prime_mod.cpp) - linear systems with kernel basis, rank and matrix inverse modulo a prime [✔](judge_tests/linear_algebra/inverse_matrix.cpp)
- [Matrix](code_library/linear_algebra/matrix.cpp) - modular multiplication and exponentiation [✔](judge_tests/linear_algebra/matrix_product.cpp)
- [Matrix Permanent](code_library/linear_algebra/permanent.cpp) - Ryser's formula with Gray code, modulo any m or exact, `O(2^n n)`
- [Simplex](code_library/linear_algebra/simplex.cpp) - linear programming
- [Thomas Algorithm](code_library/linear_algebra/thomas_algorithm.cpp) - tridiagonal linear systems in `O(n)`
- [XOR Basis](code_library/linear_algebra/xor_basis.cpp) - 64-bit linear basis with max/min xor, k-th smallest, count below x and intersection [✔](judge_tests/linear_algebra/intersection_of_f2_vector_spaces.cpp)

### Geometry

[Practice problems](code_library/geometry/README.md)

- [Circle Geometry](code_library/geometry/circle.cpp) - intersections, tangents, enclosing circle, polygon/lens/union areas, max cover [✔](judge_tests/geometry/aoj_0090.cpp)
- [Closest Pair of Points](code_library/geometry/closest_pair.cpp) - `O(n log n)` sweep with exact integer squared distances [✔](judge_tests/geometry/aoj_cgl_5_a.cpp)
- [Convex Hull](code_library/geometry/convex_hull.cpp) - monotone chain, `O(n log n)` [✔](judge_tests/geometry/aoj_cgl_3_b.cpp)
- [Geometry](code_library/geometry/geometry.cpp) - exact predicates, polygon queries/cuts, angular sort, intersections, centroid [✔](judge_tests/geometry/aoj_0177.cpp)
- [Half-plane Intersection](code_library/geometry/halfplane_intersection.cpp) - convex region of half-planes clipped to a box, sort and deque, `O(n log n)` [✔](judge_tests/geometry/aoj_1283.cpp)
- [Pick's Theorem](code_library/geometry/picks_theorem.cpp) - lattice points inside and on a polygon
- [Polygon Union](code_library/geometry/polygon_union.cpp) - area of a union of simple polygons, `O(N^2 log N)` for N vertices
- [Rectangle Union](code_library/geometry/rectangle_union.cpp) - area and perimeter of a union of axis-aligned rectangles, `O(n log n)` [✔](judge_tests/geometry/area_of_union_of_rectangles.cpp)
- [Segment Intersection Sweep](code_library/geometry/segment_intersection_sweep.cpp) - Shamos-Hoey, finds two intersecting segments exactly, `O(n log n)`

### Dynamic Programming

[Practice problems](code_library/dp/README.md)

- [Aliens Trick (WQS)](code_library/dp/aliens_trick.cpp) - exactly-k optimum of a convex cost by binary search on a penalty, tie-safe
- [Blocks](code_library/dp/blocks_dp.cpp) - interval DP for UVA 10559 Blocks
- [Bounded Knapsack](code_library/dp/bounded_knapsack.cpp) - max value with up to `counts[i]` copies of item i, `O(n capacity)` [✔](judge_tests/dp/aoj_dpl_1_g.cpp)
- [Bounded Subset Sum](code_library/dp/bounded_subset_sum.cpp) - reachable sums of a bounded multiset in `O(W sqrt(W))`, with a witness
- [CKY](code_library/dp/cky.cpp) - context-free grammar membership in Chomsky normal form
- [Concave 1D1D DP](code_library/dp/concave_1d1d_dp.cpp) - `dp[x] = min dp[i] + w(i, x)` in `O(n log n)` under the quadrangle inequality
- [Divide and Conquer DP](code_library/dp/divide_conquer_dp.cpp) - `dp[k][i] = min dp[k-1][j] + cost(j, i)` with monotone optimum
- [Hamiltonian Path and Cycle](code_library/dp/hamiltonian_dp.cpp) - shortest Hamiltonian path and cycle with the order, `O(2^n n^2)` [✔](judge_tests/dp/aoj_dpl_2_a.cpp)
- [Knuth Optimization](code_library/dp/knuth_optimization.cpp) - `dp[i][j] = min dp[i][k] + dp[k][j] + cost(i, j)` in `O(n^2)` [✔](judge_tests/dp/aoj_alds1_10_d.cpp)
- [Longest Common Increasing Subsequence](code_library/dp/lcis.cpp) - LCIS in `O(n m)`
- [Longest Increasing Subsequence](code_library/dp/lis.cpp) - LIS and LDS with one optimal subsequence, any comparator [✔](judge_tests/dp/longest_increasing_subsequence.cpp)
- [Maximum Square](code_library/dp/maximum_square.cpp) - largest filled square and diamond ending at every cell [✔](judge_tests/dp/aoj_dpl_3_a.cpp)
- [Slope Trick](code_library/dp/slope_trick.cpp) - convex piecewise linear functions with hinge adds, prefix/suffix min and shifts

### Miscellaneous

[Practice problems](code_library/misc/README.md)

- [15 Puzzle Solver](code_library/misc/15_puzzle_solver.cpp) - solvability check and IDA* solver
- [Assembly](code_library/misc/assembly.cpp) - inline x86 popcount, leading zeros, bit scan, gcd and square root [✔](judge_tests/misc/aoj_alds1_1_b.cpp)
- [Bit Twiddling](code_library/misc/bit_twiddling.cpp) - bit manipulation tricks [✔](judge_tests/misc/aoj_itp2_11_d.cpp)
- [Comb Sort](code_library/misc/combsort.cpp) - comb sort on any random access range [✔](judge_tests/misc/aoj_alds1_6_a.cpp)
- [Dancing Links](code_library/misc/dancing_links.cpp) - exact cover with Algorithm X [✔](judge_tests/misc/aoj_alds1_13_a.cpp)
- [Fast I/O](code_library/misc/fast_io.cpp) - buffered input and output with fread and fwrite [✔](judge_tests/misc/many_aplusb.cpp)
- [Gray Codes](code_library/misc/gray_codes.cpp) - Gray code and its inverse
- [Hackenbush](code_library/misc/hackenbush.cpp) - green Hackenbush Grundy values on graphs, exact red-blue values on trees
- [Knight Distance](code_library/misc/knight_distance.cpp) - fewest knight moves on an infinite board in `O(1)`
- [N Queens](code_library/misc/n_queen.cpp) - number of N queens solutions
- [Next Palindrome](code_library/misc/next_palindrome.cpp) - smallest palindromic number above a given one
- [Numerical Integration](code_library/misc/numerical_integration.cpp) - composite Simpson, adaptive Simpson and Romberg
- [Radix Sort](code_library/misc/radix_sort.cpp) - LSD radix sort of 32-bit integers [✔](judge_tests/misc/aoj_alds1_6_a_radix_sort.cpp)
- [Radix Sort 64](code_library/misc/radix_sort_64.cpp) - LSD radix sort of 64-bit integers [✔](judge_tests/misc/aoj_alds1_6_a_radix_sort_64.cpp)
- [Simulated Annealing](code_library/misc/simulated_annealing.cpp) - generic minimizer with a geometric time or iteration cooling schedule
- [Ternary and Golden-Section Search](code_library/misc/ternary_search.cpp) - integer argmax of a unimodal function, real argmin

### Hacking

[Practice problems](code_library/hacking/README.md)

- [Anti Double Hash](code_library/hacking/anti_double_hash.cpp) - two strings colliding under a double polynomial hash
- [Anti Polynomial Hash](code_library/hacking/anti_polyonmial_hash.cpp) - two strings colliding under a polynomial hash

### Python

[Practice problems](code_library/python/README.md)

- [Alpha Beta Pruning](code_library/python/alpha_beta_pruning.py) - minimax game search with alpha-beta pruning
- [Berlekamp Massey](code_library/python/berlekamp_massey.py) - shortest linear recurrence of a sequence modulo a prime
- [Derangements](code_library/python/derangements.py) - generalized derangement counts
- [Dijkstra](code_library/python/dijkstra.py) - single-source shortest paths
- [Discrete Log](code_library/python/discrete_log.py) - smallest x with `a^x = b mod m`, any m
- [Fast Fibonacci](code_library/python/fast_fibonacci.py) - n-th Fibonacci number by fast doubling
- [Fast I/O](code_library/python/fast_io.py) - token scanner over stdin
- [General Graph Matching](code_library/python/graph_matching.py) - maximum matching via the Tutte matrix rank
- [Lagrange Interpolation](code_library/python/lagrange_polynomial_interpolation.py) - value of a polynomial from any samples
- [Miller Rabin](code_library/python/miller_rabin.py) - deterministic primality test for 64-bit integers
- [N Queens](code_library/python/n_queen.py) - number of N queens solutions
- [Sieve](code_library/python/sieve.py) - sieve of Eratosthenes
- [String Hash](code_library/python/string_hash.py) - polynomial substring hashes
- [Z Algorithm](code_library/python/z_algorithm.py) - longest common prefix of every suffix with the string

## Contributing

Bug reports, fixes and new templates are welcome as issues or pull requests. See [CONTRIBUTING.md](CONTRIBUTING.md) for the conventions, the checks CI runs and The Zen Of Contributing.

## License

MIT, see [LICENSE](LICENSE)
