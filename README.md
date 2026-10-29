# Algovault

**Algorithms, data structures and templates for competitive programming**

[![CI](https://github.com/sgtlaugh/algovault/actions/workflows/ci.yml/badge.svg)](https://github.com/sgtlaugh/algovault/actions/workflows/ci.yml)
[![license](https://img.shields.io/github/license/sgtlaugh/algovault.svg?style=flat)](https://github.com/sgtlaugh/algovault/blob/master/LICENSE)

- A black box you can copy easily: simple, fast, and well-tested
- C++17 for GCC, with a few in Python 3 that mostly run on PyPy too
- Built for contests like [Codeforces](https://codeforces.com), [AtCoder](https://atcoder.jp) and ICPC

## How To Use

1. Copy a file into your solution, every C++ file compiles on its own
2. Delete its `main()`, which only holds self-tests and a worked example
3. Read the header at the top for the API, complexity and limits

For example, a [Suffix Automaton](code_library/strings/suffix_automaton.cpp) answers substring questions in three calls:

```cpp
SuffixAutomaton<> sa("banana");

sa.distinct_substrings();   // 15
sa.occurrences("ana");      // 2, overlapping
sa.contains("nab");         // false
```

## Index

- ✔ marks a template that passes a judge test, and links to one
- ✔ next to a practice problem links an accepted solution built on the template, kept in [`practice/`](practice)

| **[Topics](#topics)** | | | | |
|---|---|---|---|---|
| [Algebra](#algebra) | [Combinatorics](#combinatorics) | [Data Structures](#data-structures) | [Dynamic Programming](#dynamic-programming) | [Geometry](#geometry) |
| [Graphs](#graphs) | [Linear Algebra](#linear-algebra) | [Number Theory](#number-theory) | [Strings](#strings) | [Trees](#trees) |

| **[More](#more)** | | |
|---|---|---|
| [Hacking](#hacking) | [Miscellaneous](#miscellaneous) | [Python](#python) |

## Topics

### Algebra

[Practice problems](code_library/algebra/README.md)

- [Big Integer](code_library/algebra/bignum.cpp) - arbitrary precision signed integers [✔](judge_tests/algebra/addition_of_big_integers.cpp)
- [FFT](code_library/algebra/fft.cpp) - polynomial products, fastest for small values [✔](judge_tests/algebra/convolution_mod_1000000007_fft.cpp)
- [Fraction](code_library/algebra/fraction.cpp) - exact rational arithmetic, always reduced
- [GCD and LCM Convolution](code_library/algebra/gcd_lcm_convolution.cpp) - convolution over `gcd(i, j)` or `lcm(i, j)` [✔](judge_tests/algebra/gcd_convolution.cpp)
- [Lagrange Interpolation](code_library/algebra/lagrange_interpolation.cpp) - value of a polynomial from consecutive samples [✔](judge_tests/algebra/sum_of_exponential_times_polynomial.cpp)
- [Linear Recurrence](code_library/algebra/linear_recurrence.cpp) - linear recurrence of a sequence and its n-th term [✔](judge_tests/algebra/find_linear_recurrence.cpp)
- [Min-Plus Convolution](code_library/algebra/min_plus_convolution.cpp) - fast when one or both sides are convex [✔](judge_tests/algebra/min_plus_convolution_convex_arbitrary.cpp)
- [NTT](code_library/algebra/ntt.cpp) - polynomial products, fastest modulo NTT primes [✔](judge_tests/algebra/convolution_mod.cpp)
- [Online Convolution](code_library/algebra/online_convolution.cpp) - terms built from a convolution of earlier ones [✔](judge_tests/algebra/exp_of_formal_power_series_online.cpp)
- [Polynomial](code_library/algebra/polynomial.cpp) - power series ops and multipoint evaluation [✔](judge_tests/algebra/division_of_polynomials.cpp)
- [Subset Convolution and SOS DP](code_library/algebra/subset_convolution.cpp) - zeta and Mobius transforms, subset convolution [✔](judge_tests/algebra/bitwise_and_convolution.cpp)
- [Walsh-Hadamard Transform](code_library/algebra/walsh_hadamard.cpp) - `xor`, `or` and `and` convolutions [✔](judge_tests/algebra/bitwise_and_convolution_fwht.cpp)

### Combinatorics

[Practice problems](code_library/combinatorics/README.md)

- [Bell Numbers](code_library/combinatorics/bell_numbers.cpp) - set partition counts, also one `B_n` for huge n [✔](judge_tests/combinatorics/bell_number.cpp)
- [Bernoulli Numbers](code_library/combinatorics/bernoulli_numbers.cpp) - `B_0..B_n` modulo a prime
- [Binomial Coefficients](code_library/combinatorics/binomial_coefficients.cpp) - n choose k modulo any m, prime or not [✔](judge_tests/combinatorics/binomial_coefficient.cpp)
- [Combinatorics](code_library/combinatorics/combinatorics.cpp) - factorial tables, nCr and nPr modulo a prime [✔](judge_tests/combinatorics/binomial_coefficient_prime_mod.cpp)
- [Eulerian Numbers](code_library/combinatorics/eulerian_numbers.cpp) - permutations with k ascents, modulo any m
- [Faulhaber's Formula](code_library/combinatorics/faulhaber_formula.cpp) - sum of k-th powers of `1..n`
- [Fast Factorial](code_library/combinatorics/fast_factorial.cpp) - `n!` modulo a prime p in `O(sqrt(p) log p)` [✔](judge_tests/combinatorics/factorial.cpp)
- [Gambler's Ruin](code_library/combinatorics/gamblers_ruin.cpp) - probability that the first player goes broke
- [Josephus Problem](code_library/combinatorics/josephus_problem.cpp) - the m-th one out when every k-th person leaves [✔](judge_tests/combinatorics/aoj_0085.cpp)
- [Necklaces](code_library/combinatorics/necklaces.cpp) - ring colorings up to rotation or reflection
- [Partition Numbers](code_library/combinatorics/partition_numbers.cpp) - `p(0..n)` modulo any m [✔](judge_tests/combinatorics/partition_function.cpp)
- [Permutation Cycles](code_library/combinatorics/permutation_cycles.cpp) - cycles, k-th power and root, minimum swaps [✔](judge_tests/combinatorics/aoj_alds1_6_d.cpp)
- [Permutation Index](code_library/combinatorics/permutation_index.cpp) - `O(n)` rank of a permutation, a perfect hash
- [Permutation Rank](code_library/combinatorics/permutation_rank.cpp) - rank and unrank of permutations [✔](judge_tests/combinatorics/aoj_itp2_5_c.cpp)
- [Stirling Numbers](code_library/combinatorics/stirling_numbers.cpp) - a whole row of either kind modulo m [✔](judge_tests/combinatorics/stirling_number_of_the_first_kind.cpp)

### Data Structures

[Practice problems](code_library/data_structures/README.md)

- [Convex Hull Trick](code_library/data_structures/convex_hull_trick.cpp) - min or max of lines added in slope order
- [Coordinate Compression](code_library/data_structures/coordinate_compression.cpp) - replaces values with small ids
- [Disjoint Set Union](code_library/data_structures/disjoint_set.cpp) - merges sets and tells if two nodes share one [✔](judge_tests/data_structures/unionfind.cpp)
  - [Partially Persistent](code_library/data_structures/partially_persistent_dsu.cpp) - connectivity and component sizes at any past time
  - [Potential](code_library/data_structures/potential_dsu.cpp) - differences `x[b] - x[a]` within each set [✔](judge_tests/data_structures/unionfind_with_potential.cpp)
  - [Rollback](code_library/data_structures/rollback_dsu.cpp) - unions can be undone, for offline connectivity [✔](judge_tests/data_structures/persistent_unionfind.cpp)
- [Distinct Subarray Aggregates](code_library/data_structures/distinct_subarray_aggregates.cpp) - distinct `gcd`/`or`/`and` of subarrays ending at i
- [Euler Tour Tree](code_library/data_structures/euler_tour_tree.cpp) - dynamic forest with subtree sums and subtree add [✔](judge_tests/data_structures/dynamic_tree_subtree_add_subtree_sum.cpp)
- [Fenwick Tree](code_library/data_structures/fenwick_tree.cpp) - point or range add with point or range sums [✔](judge_tests/data_structures/aoj_dsl_2_e.cpp)
  - [2D](code_library/data_structures/fenwick_tree_2D.cpp) - the same on a grid [✔](judge_tests/data_structures/aoj_dsl_5_b.cpp)
  - [2D, Implicit](code_library/data_structures/fenwick_tree_2D_implicit.cpp) - rectangle add and cell value on huge grids
  - [2D, Sparse](code_library/data_structures/fenwick_tree_2D_sparse.cpp) - point add and rectangle sum on huge grids [✔](judge_tests/data_structures/point_add_rectangle_sum.cpp)
  - [3D](code_library/data_structures/fenwick_tree_3D.cpp) - the same in 3D
- [Hash Map](code_library/data_structures/hashmap.cpp) - anti-hack hash table for integer keys [✔](judge_tests/data_structures/associative_array.cpp)
- [Interval Set](code_library/data_structures/interval_set.cpp) - interval union with add, remove, coverage [✔](judge_tests/data_structures/aoj_dsl_4_a.cpp)
- [Leftist Heap, Persistent](code_library/data_structures/leftist_heap.cpp) - persistent meldable priority queue [✔](judge_tests/data_structures/k_shortest_walk.cpp)
- [Li Chao Tree](code_library/data_structures/li_chao_tree.cpp) - min or max of lines and segments at a point [✔](judge_tests/data_structures/line_add_get_min.cpp)
- [Link-Cut Tree](code_library/data_structures/link_cut_tree.cpp) - dynamic forest with link, cut and path queries [✔](judge_tests/data_structures/dynamic_tree_vertex_add_path_sum.cpp)
- [Merge Sort Tree and Wavelet](code_library/data_structures/merge_sort_tree.cpp) - range k-th smallest and count below x [✔](judge_tests/data_structures/range_kth_smallest.cpp)
- [Mo's Algorithm](code_library/data_structures/mo.cpp) - offline range queries on an array [✔](judge_tests/data_structures/static_range_count_distinct.cpp)
  - [On Tree Paths](code_library/data_structures/tree_path_mo.cpp) - offline path queries on a tree
  - [With Rollback](code_library/data_structures/mo_with_rollback.cpp) - state that can add but not remove [✔](judge_tests/data_structures/static_range_mode_query.cpp)
  - [With Updates](code_library/data_structures/mo_with_updates.cpp) - range queries with point assignments [✔](judge_tests/data_structures/point_set_range_frequency.cpp)
- [Monotonic Stack](code_library/data_structures/monotonic_stack.cpp) - nearest smaller or greater and histogram rectangle [✔](judge_tests/data_structures/aoj_dpl_3_c.cpp)
- [Ordered Set](code_library/data_structures/ordered_set.cpp) - set and multiset with k-th smallest and rank [✔](judge_tests/data_structures/ordered_set.cpp)
- [Permutation Tree](code_library/data_structures/permutation_tree.cpp) - decomposes a permutation into common intervals [✔](judge_tests/data_structures/common_interval_decomposition_tree.cpp)
- [Rope](code_library/data_structures/rope.cpp) - persistent text, every version stays readable
- [Segment Tree](code_library/data_structures/segment_tree.cpp) - range update and range query, with binary search [✔](judge_tests/data_structures/point_set_range_composite.cpp)
  - [Beats](code_library/data_structures/segment_tree_beats.cpp) - range chmin and chmax with range sum [✔](judge_tests/data_structures/range_chmin_chmax_add_range_sum.cpp)
  - [Merging](code_library/data_structures/segment_tree_merge.cpp) - multisets with merge, split and k-th
  - [Persistent](code_library/data_structures/persistent_segment_tree.cpp) - every version stays queryable, with k-th search [✔](judge_tests/data_structures/persistent_queue.cpp)
  - [XOR](code_library/data_structures/xor_segment_tree.cpp) - range fold of `a[p ^ x]` for any x
- [Sliding Window Aggregation](code_library/data_structures/sliding_window_aggregation.cpp) - queue fold under any associative operation [✔](judge_tests/data_structures/aoj_dsl_3_d.cpp)
- [Sparse Table](code_library/data_structures/sparse_table.cpp) - `O(1)` range min, max or gcd on a static array [✔](judge_tests/data_structures/staticrmq.cpp)
  - [2D](code_library/data_structures/sparse_table_2D.cpp) - `O(1)` static rectangle min or max on a grid
  - [Disjoint](code_library/data_structures/disjoint_sparse_table.cpp) - `O(1)` query for any associative operation [✔](judge_tests/data_structures/static_range_sum.cpp)
- [Square Root Decomposition](code_library/data_structures/sqrt_decomposition.cpp) - range add and count of values below x in a range
- [Treap](code_library/data_structures/treap.cpp) - ordered multiset and sequence with reversals [✔](judge_tests/data_structures/double_ended_priority_queue.cpp)
  - [Run-Length Implicit](code_library/data_structures/run_length_implicit_treap.cpp) - huge sequences stored as runs of equal values
- [Trie](code_library/data_structures/trie.cpp) - prefix counts and binary trie xor queries [✔](judge_tests/data_structures/aoj_alds1_4_c.cpp)

### Dynamic Programming

[Practice problems](code_library/dp/README.md)

- [Aliens Trick (WQS)](code_library/dp/aliens_trick.cpp) - optimum with exactly k items for a convex cost
- [Block Removal DP](code_library/dp/blocks_dp.cpp) - max score clearing runs of equal blocks
- [Bounded Knapsack](code_library/dp/bounded_knapsack.cpp) - limited copies per item in `O(n capacity)` [✔](judge_tests/dp/aoj_dpl_1_g.cpp)
- [Bounded Subset Sum](code_library/dp/bounded_subset_sum.cpp) - reachable sums of a multiset, with a witness
- [CKY](code_library/dp/cky.cpp) - context-free grammar membership
- [Concave 1D1D DP](code_library/dp/concave_1d1d_dp.cpp) - `dp[x] = min dp[i] + w(i, x)` for Monge w
- [Divide and Conquer DP](code_library/dp/divide_conquer_dp.cpp) - split into k segments of minimum total cost
  - [Knuth Partition DP](code_library/dp/knuth_partition_dp.cpp) - same split in `O(n^2)`, faster for large k
- [Garsia-Wachs](code_library/dp/garsia_wachs.cpp) - min cost to merge adjacent piles in `O(n log n)`
- [Hamiltonian Path and Cycle](code_library/dp/hamiltonian_dp.cpp) - shortest Hamiltonian path and cycle [✔](judge_tests/dp/aoj_dpl_2_a.cpp)
- [Knuth Optimization](code_library/dp/knuth_optimization.cpp) - interval DP over split points in `O(n^2)` [✔](judge_tests/dp/aoj_alds1_10_d.cpp)
- [Longest Common Increasing Subsequence](code_library/dp/lcis.cpp) - LCIS and one witness in `O(n m)`
- [Longest Increasing Subsequence](code_library/dp/lis.cpp) - LIS and LDS with one optimal subsequence [✔](judge_tests/dp/longest_increasing_subsequence.cpp)
- [Maximum Square](code_library/dp/maximum_square.cpp) - largest filled square and diamond at each cell [✔](judge_tests/dp/aoj_dpl_3_a.cpp)
- [Plug DP](code_library/dp/plug_dp.cpp) - count or optimize cycles and paths on a grid
- [Slope Trick](code_library/dp/slope_trick.cpp) - maintains a convex piecewise linear function

### Geometry

[Practice problems](code_library/geometry/README.md)

- [Circle Geometry](code_library/geometry/circle.cpp) - intersections, tangents, enclosing circle, areas [✔](judge_tests/geometry/aoj_0090.cpp)
- [Closest Pair of Points](code_library/geometry/closest_pair.cpp) - nearest two points by exact squared distance [✔](judge_tests/geometry/aoj_cgl_5_a.cpp)
- [Convex Hull](code_library/geometry/convex_hull.cpp) - static and online hull, Minkowski sum [✔](judge_tests/geometry/aoj_cgl_3_b.cpp)
  - [3D](code_library/geometry/convex_hull_3D.cpp) - triangulated hull in space, volume and area
- [Delaunay Triangulation](code_library/geometry/delaunay.cpp) - triangulation with empty circumcircles [✔](judge_tests/geometry/euclidean_mst.cpp)
- [Geometry](code_library/geometry/geometry.cpp) - exact predicates and polygon queries [✔](judge_tests/geometry/aoj_0177.cpp)
- [Half-Plane Intersection](code_library/geometry/halfplane_intersection.cpp) - convex region common to all half-planes [✔](judge_tests/geometry/aoj_1283.cpp)
- [KD-Tree](code_library/geometry/kd_tree.cpp) - nearest point to a query, points in a rectangle [✔](judge_tests/geometry/aoj_dsl_2_c.cpp)
- [Pick's Theorem](code_library/geometry/picks_theorem.cpp) - lattice points inside and on a polygon
- [Polygon Union](code_library/geometry/polygon_union.cpp) - area of a union of simple polygons
- [Rectangle Union](code_library/geometry/rectangle_union.cpp) - area and perimeter of a union of rectangles [✔](judge_tests/geometry/area_of_union_of_rectangles.cpp)
- [Segment Intersection Sweep](code_library/geometry/segment_intersection_sweep.cpp) - whether any two segments intersect

### Graphs

[Practice problems](code_library/graphs/README.md)

- [2-SAT](code_library/graphs/2SAT_kosaraju.cpp) - satisfiability with an assignment and at-most-one [✔](judge_tests/graphs/two_sat.cpp)
  - [Lexicographic](code_library/graphs/2SAT_lexicographic.cpp) - satisfying assignment preferring true from x1 [✔](judge_tests/graphs/two_sat_lexicographic.cpp)
- [3-Edge-Connected Components](code_library/graphs/three_edge_connected_components.cpp) - vertices joined by 3 edge-disjoint paths [✔](judge_tests/graphs/three_edge_connected_components.cpp)
- [Absolute Center](code_library/graphs/absolute_center.cpp) - graph center and minimum diameter spanning tree [✔](judge_tests/graphs/minimum_diameter_spanning_tree.cpp)
- [Articulation Points and Blocks](code_library/graphs/articulation_points.cpp) - cut vertices, blocks and the block-cut tree [✔](judge_tests/graphs/aoj_grl_3_a.cpp)
- [Bellman-Ford](code_library/graphs/bellman_ford.cpp) - shortest paths with negative edges and cycles [✔](judge_tests/graphs/aoj_grl_1_b.cpp)
- [Bridges](code_library/graphs/bridge.cpp) - bridges and 2-edge-connected components [✔](judge_tests/graphs/aoj_grl_3_b.cpp)
- [Cactus Graph](code_library/graphs/cactus.cpp) - cactus recognition, cycles and the cactus tree
- [Chordal Graph](code_library/graphs/chordal_graph.cpp) - recognition, elimination order or induced cycle [✔](judge_tests/graphs/chordal_graph_recognition.cpp)
- [Chromatic Number](code_library/graphs/chromatic_number.cpp) - minimum colors in a vertex coloring [✔](judge_tests/graphs/chromatic_number.cpp)
- [Complement Graph BFS](code_library/graphs/complement_graph_bfs.cpp) - BFS and components of the complement graph [✔](judge_tests/graphs/connected_components_of_complement_graph.cpp)
- [Dijkstra](code_library/graphs/dijkstra.cpp) - shortest paths with non-negative edges, 0-1 BFS [✔](judge_tests/graphs/aoj_alds1_11_c.cpp)
- [Directed MST](code_library/graphs/directed_mst.cpp) - minimum arborescence from a root [✔](judge_tests/graphs/directedmst.cpp)
- [Dominator Tree](code_library/graphs/dominator_tree.cpp) - immediate dominators from a root [✔](judge_tests/graphs/dominatortree.cpp)
- [Edge Coloring](code_library/graphs/edge_coloring.cpp) - D colors if bipartite, D + 1 for simple graphs [✔](judge_tests/graphs/bipartite_edge_coloring.cpp)
- [Euler Path](code_library/graphs/euler_path.cpp) - walk using every edge exactly once [✔](judge_tests/graphs/eulerian_trail_directed.cpp)
- [Flow with Lower Bounds](code_library/graphs/lower_bound_flow.cpp) - circulation and s-t flow with edge lower bounds [✔](judge_tests/graphs/aoj_grl_6_a_lower_bound.cpp)
- [Floyd-Warshall](code_library/graphs/floyd_warshall.cpp) - all pairs shortest paths with negative cycles [✔](judge_tests/graphs/aoj_grl_1_c.cpp)
- [General Graph Matching](code_library/graphs/graph_matching.cpp) - maximum matching in general graphs [✔](judge_tests/graphs/general_matching.cpp)
  - [Weighted](code_library/graphs/weighted_graph_matching.cpp) - edge weights, maximizes the total weight [✔](judge_tests/graphs/general_weighted_matching.cpp)
- [Global Minimum Cut](code_library/graphs/global_min_cut.cpp) - minimum cut weight and one side
- [Gomory-Hu Tree](code_library/graphs/gomory_hu.cpp) - all-pairs min cuts of an undirected graph
- [Hopcroft-Karp](code_library/graphs/hopcroft_karp.cpp) - maximum bipartite matching, `O(m sqrt(n))` [✔](judge_tests/graphs/bipartitematching.cpp)
- [Hungarian Algorithm](code_library/graphs/hungarian_algorithm.cpp) - minimum or maximum cost assignment [✔](judge_tests/graphs/assignment.cpp)
- [Johnson's Algorithm](code_library/graphs/johnsons_algorithm.cpp) - all-pairs shortest paths with negative edges [✔](judge_tests/graphs/aoj_grl_1_c_johnson.cpp)
- [K Shortest Walks](code_library/graphs/k_shortest_walks.cpp) - lengths of the k shortest walks from s to t [✔](judge_tests/graphs/k_shortest_walk_k_shortest_walks.cpp)
- [Kruskal Reconstruction Tree](code_library/graphs/kruskal_reconstruction_tree.cpp) - reachable vertices using edges of weight at most w
- [Manhattan MST](code_library/graphs/manhattan_mst.cpp) - minimum spanning tree under L1 distance [✔](judge_tests/graphs/manhattanmst.cpp)
- [Matroid Intersection](code_library/graphs/matroid_intersection.cpp) - largest common independent set of two matroids [✔](judge_tests/graphs/aoj_grl_7_a.cpp)
  - [Weighted](code_library/graphs/weighted_matroid_intersection.cpp) - the max weight common independent set
- [Maximum Clique](code_library/graphs/max_clique.cpp) - max clique, max independent set, maximal cliques [✔](judge_tests/graphs/maximum_independent_set.cpp)
- [Maximum Flow](code_library/graphs/maxflow.cpp) - max flow and min cut, with node capacities [✔](judge_tests/graphs/aoj_grl_6_a.cpp)
  - [HLPP](code_library/graphs/hlpp.cpp) - push-relabel, faster where Dinic needs many phases [✔](judge_tests/graphs/aoj_grl_6_a_hlpp.cpp)
- [Min Cost Circulation](code_library/graphs/min_cost_circulation.cpp) - cheapest circulation with negative cycles [✔](judge_tests/graphs/aoj_grl_6_b_circulation.cpp)
- [Min Cost Max Flow, Dijkstra](code_library/graphs/mcmf_dijkstra.cpp) - default min cost flow, with the cost slope [✔](judge_tests/graphs/aoj_grl_6_b.cpp)
- [Min Cost Max Flow, SPFA](code_library/graphs/mcmf_spfa.cpp) - same as Dijkstra, without the cost slope [✔](judge_tests/graphs/aoj_grl_6_b_spfa.cpp)
- [Minimum Mean Cycle (Karp)](code_library/graphs/minimum_mean_cycle.cpp) - exact minimum mean cycle and its edges
- [Minimum Path Cover](code_library/graphs/minimum_path_cover.cpp) - path covers of a DAG and maximum antichain [✔](judge_tests/graphs/aoj_2251.cpp)
- [Minimum Spanning Tree](code_library/graphs/minimum_spanning_tree.cpp) - minimum spanning forest, also implicit graphs [✔](judge_tests/graphs/aoj_alds1_12_a.cpp)
- [Network Simplex](code_library/graphs/network_simplex.cpp) - min cost flow with supplies and lower bounds [✔](judge_tests/graphs/min_cost_b_flow.cpp)
- [Offline Dynamic Connectivity](code_library/graphs/offline_dynamic_connectivity.cpp) - connectivity under edge insertions and deletions
  - [Online](code_library/graphs/online_dynamic_connectivity.cpp) - answers as it goes, with component sums [✔](judge_tests/graphs/dynamic_graph_vertex_add_component_sum.cpp)
- [Online Bridges](code_library/graphs/online_bridges.cpp) - bridges and 2-edge components under insertions [✔](judge_tests/graphs/two_edge_connected_components.cpp)
- [Range Edge Dijkstra](code_library/graphs/range_edge_dijkstra.cpp) - shortest paths with edges to and from ranges
- [Small Cycle Counting](code_library/graphs/small_cycle_counting.cpp) - counts triangles and 4-cycles
- [Stable Marriage](code_library/graphs/stable_marriage.cpp) - proposer-optimal stable matching
- [Steiner Tree](code_library/graphs/steiner_tree.cpp) - minimum tree connecting k terminals [✔](judge_tests/graphs/minimum_steiner_tree.cpp)
- [Strongly Connected Components](code_library/graphs/scc.cpp) - components in topological order [✔](judge_tests/graphs/scc.cpp)

### Linear Algebra

[Practice problems](code_library/linear_algebra/README.md)

- [Characteristic Polynomial](code_library/linear_algebra/characteristic_polynomial.cpp) - `det(xI - A)` modulo a prime [✔](judge_tests/linear_algebra/characteristic_polynomial.cpp)
- [Determinant](code_library/linear_algebra/determinant.cpp) - determinant modulo any m, spanning tree counts [✔](judge_tests/linear_algebra/counting_spanning_tree_directed.cpp)
- [Freivalds' Algorithm](code_library/linear_algebra/freivalds_algorithm.cpp) - randomized check of a matrix product
- [Gauss-Jordan](code_library/linear_algebra/gauss_jordan.cpp) - linear systems, rank and inverse over the reals
  - [Band Matrix](code_library/linear_algebra/gauss_band_matrix.cpp) - grid equations that only use nearby cells
  - [Bitset](code_library/linear_algebra/gauss_bitset.cpp) - the same over GF(2) [✔](judge_tests/linear_algebra/inverse_matrix_mod_2.cpp)
  - [Prime Modulus](code_library/linear_algebra/gauss_prime_mod.cpp) - the same modulo a prime, with kernel basis [✔](judge_tests/linear_algebra/inverse_matrix.cpp)
- [Matrix](code_library/linear_algebra/matrix.cpp) - modular multiplication and exponentiation [✔](judge_tests/linear_algebra/matrix_product.cpp)
- [Matrix Permanent](code_library/linear_algebra/permanent.cpp) - permanent modulo any m, or exact
- [Simplex](code_library/linear_algebra/simplex.cpp) - linear programming
- [Thomas Algorithm](code_library/linear_algebra/thomas_algorithm.cpp) - tridiagonal linear systems in `O(n)`
- [XOR Basis](code_library/linear_algebra/xor_basis.cpp) - max xor, k-th xor and span intersection [✔](judge_tests/linear_algebra/intersection_of_f2_vector_spaces.cpp)

### Number Theory

[Practice problems](code_library/number_theory/README.md)

- [All Divisors](code_library/number_theory/all_divisors.cpp) - divisor lists of every number up to a limit
- [Chinese Remainder Theorem](code_library/number_theory/chinese_remainder_theorem.cpp) - systems of congruences with any moduli [✔](judge_tests/number_theory/convolution_mod_1000000007_chinese_remainder_theorem.cpp)
- [Continued Fractions](code_library/number_theory/continued_fractions.cpp) - best rational approximation, Stern-Brocot search [✔](judge_tests/number_theory/rational_approximation.cpp)
- [Digits of Factorial](code_library/number_theory/digits_of_factorial.cpp) - number of digits of n! in any base
- [Discrete Logarithm](code_library/number_theory/discrete_logarithm.cpp) - smallest x with `a^x = b mod m`, any m [✔](judge_tests/number_theory/discrete_logarithm_mod.cpp)
- [Divisors](code_library/number_theory/divisors.cpp) - sorted divisors of one number
- [Du's Sieve](code_library/number_theory/du_sieve.cpp) - prefix sums of phi and mu in `O(n^(2/3))` [✔](judge_tests/number_theory/counting_squarefrees.cpp)
- [Extended GCD](code_library/number_theory/extended_gcd.cpp) - modular inverse, linear diophantine equations [✔](judge_tests/number_theory/aoj_ntl_1_e.cpp)
- [Fast Fibonacci](code_library/number_theory/fast_fibonacci.cpp) - n-th Fibonacci number by fast doubling [✔](judge_tests/number_theory/aoj_alds1_10_a.cpp)
- [Fast Prime Counting](code_library/number_theory/fast_prime_counting.cpp) - count of primes up to n in sublinear time [✔](judge_tests/number_theory/counting_primes.cpp)
- [Fast Prime Sums](code_library/number_theory/fast_prime_sums.cpp) - sum of primes up to n in sublinear time
- [Fast Sieve](code_library/number_theory/fast_sieve.cpp) - lists all primes up to n [✔](judge_tests/number_theory/enumerate_primes.cpp)
- [Floor Sum](code_library/number_theory/floor_sum.cpp) - floor sums and minima of `(a i + b) mod m` [✔](judge_tests/number_theory/min_of_mod_of_linear.cpp)
- [Integer Root](code_library/number_theory/integer_root.cpp) - exact floor of square, cube and k-th roots [✔](judge_tests/number_theory/kth_root_integer.cpp)
- [Linear Sieve](code_library/number_theory/linear_sieve.cpp) - smallest prime factors, phi, mu and more
- [Maximum Divisors](code_library/number_theory/maximum_divisors.cpp) - number with the most divisors up to a limit
- [Miller-Rabin](code_library/number_theory/miller_rabin.cpp) - deterministic primality test for 64-bit integers [✔](judge_tests/number_theory/primality_test.cpp)
- [Min_25 Sieve](code_library/number_theory/min_25_sieve.cpp) - prefix sums of a multiplicative function [✔](judge_tests/number_theory/sum_of_multiplicative_function.cpp)
- [Mobius Function](code_library/number_theory/mobius_function.cpp) - Mobius mu for `1..n` at 1 byte per entry
- [Modular Roots](code_library/number_theory/modular_roots.cpp) - k-th roots modulo a prime, primitive root [✔](judge_tests/number_theory/kth_root_mod.cpp)
- [Modular Square Root](code_library/number_theory/modular_sqrt.cpp) - Tonelli-Shanks square root modulo a prime [✔](judge_tests/number_theory/sqrt_mod.cpp)
- [Pell Equation](code_library/number_theory/pell_equation.cpp) - solutions of `x^2 - d y^2 = 1` or `-1`, also mod m
- [Pisano Period](code_library/number_theory/pisano_period.cpp) - period of Fibonacci numbers modulo n
- [Pollard Rho](code_library/number_theory/pollard_rho.cpp) - factorization and primality below `2^63` [✔](judge_tests/number_theory/factorize.cpp)
- [Power Tower](code_library/number_theory/power_tower.cpp) - `a1^(a2^(...^an))` modulo m [✔](judge_tests/number_theory/tetration_mod.cpp)
- [Segmented Sieve](code_library/number_theory/segmented_sieve.cpp) - primes in a window `[L, R]` far from zero [✔](judge_tests/number_theory/enumerate_primes_segmented_sieve.cpp)
- [Sum of Two Squares](code_library/number_theory/sum_of_two_squares.cpp) - every `x^2 + y^2 = n` up to `2^63` [✔](judge_tests/number_theory/two_square_sum.cpp)

### Strings

[Practice problems](code_library/strings/README.md)

- [2D Pattern Matcher](code_library/strings/2D_pattern_matcher.cpp) - every occurrence of a 2D pattern in a 2D text [✔](judge_tests/strings/aoj_alds1_14_c.cpp)
- [Aho-Corasick](code_library/strings/aho_corasick.cpp) - multi-pattern matching automaton [✔](judge_tests/strings/aho_corasick.cpp)
- [Hirschberg LCS](code_library/strings/hirschberg_lcs.cpp) - one longest common subsequence in linear memory [✔](judge_tests/strings/aoj_alds1_10_c_hirschberg_lcs.cpp)
  - [Dynamic](code_library/strings/dynamic_aho_corasick.cpp) - patterns can be inserted online
- [Bit-String LCS](code_library/strings/bit_string_lcs.cpp) - longest common subsequence in `O(n m / 64)` [✔](judge_tests/strings/aoj_alds1_10_c.cpp)
- [De Bruijn Sequence](code_library/strings/de_bruijn_sequence.cpp) - cyclic sequence containing every length n word
- [Dynamic String Hash](code_library/strings/dynamic_string_hash.cpp) - substring hashes under range assignment [✔](judge_tests/strings/aoj_alds1_14_b_dynamic_string_hash.cpp)
- [Hashing](code_library/strings/hashing.cpp) - forward and reverse polynomial hash of any segment [✔](judge_tests/strings/enumerate_palindromes_hashing.cpp)
- [Hunt-Szymanski](code_library/strings/hunt_szymanski.cpp) - LCS in `O((r + n) log n)` for r matching pairs [✔](judge_tests/strings/aoj_alds1_10_c_hunt_szymanski.cpp)
- [KMP](code_library/strings/kmp.cpp) - prefix function and pattern search [✔](judge_tests/strings/aoj_alds1_14_b.cpp)
- [Lyndon Factorization](code_library/strings/lyndon_factorization.cpp) - splits a sequence into Lyndon words [✔](judge_tests/strings/lyndon_factorization.cpp)
- [Manacher](code_library/strings/manacher.cpp) - longest palindrome at every center [✔](judge_tests/strings/enumerate_palindromes.cpp)
- [Minimum Rotation](code_library/strings/minimum_rotation.cpp) - start of the lexicographically smallest rotation
- [Palindromic Tree](code_library/strings/palindromic_tree.cpp) - eertree, one node per distinct palindrome [✔](judge_tests/strings/eertree.cpp)
- [Subsequence Automaton](code_library/strings/subsequence_automaton.cpp) - next occurrence table, distinct subsequences [✔](judge_tests/strings/number_of_subsequences.cpp)
- [Suffix Array](code_library/strings/suffix_array.cpp) - `O(n)` suffix array with LCP and pattern search [✔](judge_tests/strings/aoj_alds1_14_d.cpp)
- [Suffix Automaton](code_library/strings/suffix_automaton.cpp) - substring counts and longest common substring [✔](judge_tests/strings/number_of_substrings.cpp)
- [Tandem Repeats](code_library/strings/tandem_repeats.cpp) - every tandem repeat (square) in a string [✔](judge_tests/strings/runenumerate.cpp)
- [Wildcard Matching](code_library/strings/wildcard_matching.cpp) - matching with wildcards in text or pattern [✔](judge_tests/strings/wildcard_pattern_matching.cpp)
- [Z Algorithm](code_library/strings/z_algorithm.cpp) - LCP of every suffix with the string [✔](judge_tests/strings/zalgorithm.cpp)

### Trees

[Practice problems](code_library/trees/README.md)

- [Centroid Decomposition](code_library/trees/centroid_decomposition.cpp) - centroid tree with distances to each centroid [✔](judge_tests/trees/frequency_table_of_tree_distance.cpp)
- [DSU on Tree](code_library/trees/dsu_on_tree.cpp) - small-to-large answers for every subtree
- [Dynamic Tree Diameter](code_library/trees/dynamic_diameter.cpp) - tree diameter under edge weight updates [✔](judge_tests/trees/aoj_grl_5_a.cpp)
- [Heavy Light Decomposition](code_library/trees/hld.cpp) - tree paths and subtrees as array ranges [✔](judge_tests/trees/aoj_grl_5_d.cpp)
- [Long-Path Decomposition](code_library/trees/long_path_decomposition.cpp) - depth-indexed subtree DP and `O(1)` k-th ancestor [✔](judge_tests/trees/jump_on_tree_long_path_decomposition.cpp)
- [Lowest Common Ancestor](code_library/trees/lca.cpp) - LCA, distance, jumps and path folds on a tree [✔](judge_tests/trees/jump_on_tree.cpp)
- [Prufer Code](code_library/trees/prufer_code.cpp) - converts labeled trees to and from sequences
- [Rerooting DP](code_library/trees/rerooting.cpp) - a tree DP with every node as the root [✔](judge_tests/trees/aoj_grl_5_b.cpp)
- [Tree Isomorphism](code_library/trees/tree_isomorphism.cpp) - canonical ids for rooted and unrooted trees [✔](judge_tests/trees/rooted_tree_isomorphism_classification.cpp)
- [Virtual Tree](code_library/trees/virtual_tree.cpp) - compresses a tree to k vertices and their LCAs

## More

### Hacking

[Practice problems](code_library/hacking/README.md)

- [Anti Double Hash](code_library/hacking/anti_double_hash.cpp) - two strings colliding under a double hash
- [Anti Polynomial Hash](code_library/hacking/anti_polynomial_hash.cpp) - two strings colliding under a polynomial hash

### Miscellaneous

[Practice problems](code_library/misc/README.md)

- [15 Puzzle Solver](code_library/misc/15_puzzle_solver.cpp) - solvability check and IDA* solver
- [Assembly](code_library/misc/assembly.cpp) - inline x86 bit ops, gcd and square root [✔](judge_tests/misc/aoj_alds1_1_b.cpp)
- [Bit Twiddling](code_library/misc/bit_twiddling.cpp) - bit manipulation tricks [✔](judge_tests/misc/aoj_itp2_11_d.cpp)
- [Comb Sort](code_library/misc/combsort.cpp) - comb sort on any random access range [✔](judge_tests/misc/aoj_alds1_6_a.cpp)
- [Dancing Links](code_library/misc/dancing_links.cpp) - exact cover with Algorithm X [✔](judge_tests/misc/aoj_alds1_13_a.cpp)
- [Fast I/O](code_library/misc/fast_io.cpp) - fast buffered input and output [✔](judge_tests/misc/many_aplusb.cpp)
- [Gray Codes](code_library/misc/gray_codes.cpp) - Gray code and its inverse
- [Hackenbush](code_library/misc/hackenbush.cpp) - green values on graphs, red-blue on trees
- [Knight Distance](code_library/misc/knight_distance.cpp) - fewest knight moves on an infinite board in `O(1)`
- [Local Minimization](code_library/misc/local_minimization.cpp) - minimum of a function of a few real variables
- [N Queens](code_library/misc/n_queen.cpp) - number of N queens solutions
- [Next Palindrome](code_library/misc/next_palindrome.cpp) - smallest palindromic number above a given one
- [Nimber Arithmetic](code_library/misc/nimber.cpp) - nim product and inverse over GF(2^64) [✔](judge_tests/misc/nim_product_64.cpp)
- [Numerical Integration](code_library/misc/numerical_integration.cpp) - definite integrals by Simpson or Romberg
- [Radix Sort](code_library/misc/radix_sort.cpp) - LSD radix sort of 32-bit integers [✔](judge_tests/misc/aoj_alds1_6_a_radix_sort.cpp)
- [Radix Sort 64](code_library/misc/radix_sort_64.cpp) - LSD radix sort of 64-bit integers [✔](judge_tests/misc/aoj_alds1_6_a_radix_sort_64.cpp)
- [Simulated Annealing](code_library/misc/simulated_annealing.cpp) - generic randomized minimizer
- [Ternary and Golden-Section Search](code_library/misc/ternary_search.cpp) - integer argmax of a unimodal function, real argmin

### Python

[Practice problems](code_library/python/README.md)

- [Alpha Beta Pruning](code_library/python/alpha_beta_pruning.py) - minimax game search with alpha-beta pruning
- [Berlekamp-Massey](code_library/python/berlekamp_massey.py) - shortest linear recurrence and its n-th term
- [Derangements](code_library/python/derangements.py) - generalized derangement counts
- [Dijkstra](code_library/python/dijkstra.py) - single-source shortest paths
- [Discrete Log](code_library/python/discrete_log.py) - smallest x with `a^x = b mod m`, any m
- [Fast Fibonacci](code_library/python/fast_fibonacci.py) - n-th Fibonacci number by fast doubling
- [Fast I/O](code_library/python/fast_io.py) - token scanner over stdin
- [General Graph Matching](code_library/python/graph_matching.py) - maximum matching size in a general graph
- [Lagrange Interpolation](code_library/python/lagrange_polynomial_interpolation.py) - value of a polynomial from any samples
- [Miller-Rabin](code_library/python/miller_rabin.py) - deterministic primality test for 64-bit integers
- [N Queens](code_library/python/n_queen.py) - number of N queens solutions
- [Sieve](code_library/python/sieve.py) - primes up to n by Eratosthenes
- [String Hash](code_library/python/string_hash.py) - polynomial substring hashes
- [Z Algorithm](code_library/python/z_algorithm.py) - LCP of every suffix with the string

## Testing

Every template is checked four ways:

| Layer | Where | Checks |
|---|---|---|
| Self-tests | `main()` in each file | Usage examples with known answers |
| Stress tests | [`stress_tests/`](stress_tests) | Random inputs against an independent brute force |
| Repo checks | [`run_tests.sh`](.github/scripts/run_tests.sh) | Every file has a stress test and an index entry, consistent style |
| Judge tests | [`judge_tests/`](judge_tests) | Official [Library Checker](https://judge.yosupo.jp) and [Aizu](https://onlinejudge.u-aizu.ac.jp) test data at full size |

### How tests build

- **Self and stress tests:** GCC 14
  - `-fsanitize=address,undefined`: AddressSanitizer and UndefinedBehaviorSanitizer
  - an 8 MB stack, as on most judges
  - `-Werror`: warnings as errors
- **Judge tests:** `-O2`, like a judge

### When they run

- **Every push:** the self, stress and judge tests of the files it changes, plus the repo checks
- **Every night:** all self and stress tests, then the stress tests again at 20x the iterations with a fresh seed
- **Every week:** all judge tests from scratch

The commands to run everything locally are in [CONTRIBUTING.md](CONTRIBUTING.md#checks).

## Contributing

- Bug reports, fixes and new templates are welcome as issues or pull requests, see [CONTRIBUTING.md](CONTRIBUTING.md) for the conventions and the checks CI runs
- The skills in [`.claude/skills/`](.claude/skills) apply them:
  - `/crucible` adds or refactors a template
  - `/lucid` keeps the docs consistent

## License

MIT, see [LICENSE](LICENSE)
