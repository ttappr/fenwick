# fenwick

A simple generic implementation for a Fenwick Tree.

Fenwick Trees are useful where a prefix sum is required but updating it is too 
costly. The Fenwick Tree balances out the costs of update and query operations,
each with `O(log n)` time complexity.

Although the code is derived from the 1-based version of the tree, it's been
modified to be 0-based.

The debug build has `debug_assert`'s to ensure parameters are within range for
some functions and all values are positive for `.rank_query()` and 
`.min_rank_query()`. This may make the query functions seem slow for debug, but
the asserts aren't included in the release build, which will execute much faster.

## Operations

 * `new(<size>)` - Create a new tree with an internal array of the given size.
 * `from_slice(<slice>)` - Create a new tree from a slice. An `O(n)` operation.
 * `prefix_sum(<idx>)` - Get the prefix sum of all elements up to idx inclusive.
 * `end()` - The index of the array's last element.
 * `total()` - The prefix sum of all elements - an `O(1)` operation.
 * `add(<idx>, <delta>)` - Add delta to element at idx.
 * `sub(<idx>, <delta>)` - Subtract delta from element at idx.
 * `set(<idx>, <value>)` - Set element at idx to value.
 * `get(<idx>)` - Get value of element at idx.
 * `range_sum(<idx_i>, <idx_j>)` - Get sum of elements from idx_i to idx_j inclusive.
 * `rank_query(<value>)` - Find the largest index with `.prefix_sum(index) <= value`.
 * `min_rank_query(<value>)` - Find the smallest index with `.prefix_sum(index) >= value`.
 
## Example

Using a Fenwick Tree prefix sum to produce a weighted random sample without
replacement.

```rust
use rand::prelude::*;
use fenwick::*;

/// Produce a weighted sample from the population without replacement. 
/// Shows an application of a Fenwick Tree prefix sum.
///
/// # Arguments
/// * `population`  - The population to select from.
/// * `weights`     - The weights corresponding to each member of the 
///                   population.
/// * `k`           - The sample size.
///
/// # Returns
/// * A vector holding the members of the population randomly selected.
///
fn wsample(population: &[i32], weights: &[i32], k: usize) -> Vec<i32>
{
    let mut prefix_sum = Fenwick::from_slice(weights);  // Create from slice.
    let mut rng        = rand::thread_rng();
    let mut sample     = vec![];
    
    for _ in 0..k {
        let r   = (rng.gen::<f32>() * prefix_sum.total() as f32).ceil();
        
        let sel = prefix_sum.min_rank_query(r as i32);  // Querying by value.
        prefix_sum.set(sel, 0);                         // Eliminate selection.
        
        sample.push(population[sel]);
    }
    sample
}


fn main() {
    let sample = wsample(&[1, 2, 3, 4, 5], &[5, 80, 50, 1, 20], 2);
    
    println!("{:?}", sample);
}
```
 
In the code above, the tree acts as a prefix sum of the accumulated weights
of the population members. The greater the weight, the better the odds it will
be included in the sample taken.

A random number is generated from 0 to the prefix sum total (all weights added
together). Then this generated value is used to look up the index of the element 
in the tree that represents that portion of the range. Some indices represent
wider ranges than others which increases their odds of the random value falling
within their range.

After the selected item is found, its value is then set to 0 to remove it from
subsequent samples (sampling without replacement).
 
