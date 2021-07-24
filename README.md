# fenwick

A simple generic implementation for a Fenwick Tree.

Fenwick Trees are useful where a prefix sum is required but updating it is too 
costly. The Fenwick Tree balances out the costs of update and query operations,
each with `O(log n)` time complexity.

## Operations

 * `new(<size>)` - Create a new tree with an internal array of the given size.
 * `prefix_sum(<idx>)` - Get the prefix sum of all elements up to idx inclusive.
 * `end()` - The index of the array's last element.
 * `add(<idx>, <delta>)` - Add delta to element at idx.
 * `sub(<idx>, <delta>)` - Subtract delta from element at idx.
 * `set(<idx>, <value>)` - Set element at idx to value.
 * `get(<idx>)` - Get value of element at idx.
 * `range_sum(<idx_i>, <idx_j>)` - Get sum of elements from idx_i + 1 to idx_j.
 * `rank_query(<value>)` - Find the largest index with `.prefix_sum(index) <= value`.
 * `min_rank_query(<value>)` - Find the smallest index with `.prefix_sum(index) >= value`.
 
 
 
