# fenwick

A simple generic implementation of a 0-based Fenwick Tree (or BITS).

Fenwick Trees are useful where a prefix sum is required but updating it is too 
costly. The Fenwick Tree balances out the costs of update and query operations,
each with `O(log n)` time complexity.

Some functions include `debug_assert`'s to ensure parameters are within range.
`.rank_query()` and `.min_rank_query()` have debug asserts to ensure all 
elements are non-negative. These asserts are excluded when the code is compiled 
with optimizations, so release builds will have significantly faster 
performance.

## Operations

 * `new(<size>)` - Create a new tree with an internal array of at least the 
 given size (1 + a power of 2).
 * `from_slice(<slice>)` - Create a new tree from a slice. An `O(n)` operation.
 * `prefix_sum(<idx>)` - Get the prefix sum of all elements up to idx inclusive.
 * `end()` - The index of the array's last element.
 * `total()` - The prefix sum of all elements.
 * `add(<idx>, <delta>)` - Add delta to element at idx.
 * `sub(<idx>, <delta>)` - Subtract delta from element at idx.
 * `set(<idx>, <value>)` - Set element at idx to value.
 * `get(<idx>)` - Get value of element at idx.
 * `range_sum(<start>, <end>)` - Get sum of elements from start to end 
 inclusive.
 * `rank_query(<value>)` - Find the largest index with 
 `.prefix_sum(index) <= value`.
 * `min_rank_query(<value>)` - Find the smallest index with 
 `.prefix_sum(index) >= value`.
 
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
        
        // Query for first index >= r and effectively remove it.
        let sel = prefix_sum.min_rank_query(r as i32).unwrap();
        prefix_sum.set(sel, 0);
        
        sample.push(population[sel]);
    }
    sample
}


fn main() {
    let sample = wsample(&[1, 2, 3, 4, 5], &[5, 80, 50, 1, 20], 2);
    
    println!("{:?}", sample);
}
```
 
In the code above, the tree acts as a prefix sum for the cummulative weights
of the population members. For each member, the greater its weight, the better
the odds it will be included in the sample taken.

A random number is generated from 0 to the prefix sum total (all weights added
together). This generated value is then used to look up the index of the element 
corresponding to that portion of the overall range.

After a member's index is selected, its value is set to 0 to remove it from 
subsequent samples (sampling without replacement).
 
