
//! A generic implementation for a Fenwick Tree, useful for managing prefix
//! sums. Each operation has `O(log n)` time-complexity. Some operations like
//! `.end()` have `O(1)` time-complexity.
//!
//! I created this simple lib for use in solving coding challenge type problems
//! that benefit from prefix sums requiring fast operations. The code was 
//! originally taken from the Wikipedia article on Fenwick Trees and modified
//! to create a class. The code also has some fixes for indexing problems in
//! the original code.
//!
//! Wikipedia article: https://en.wikipedia.org/wiki/Fenwick_tree
//!

use std::ops::Add;
use std::ops::Sub;
use std::ops::AddAssign;
use std::ops::SubAssign;
use std::cmp::PartialOrd;

macro_rules! lsb {
    ($i:expr) => (($i) & -($i))
}

macro_rules! lsb_usize {
    ($i:expr) => { lsb!($i as isize) as usize }
}

/// Represents a prefix sum array with `O(log n)` update operations.
///
#[derive(Debug)]
pub struct Fenwick<T>
{
    data: Vec<T>,
    size: usize,
}

impl<T> Fenwick<T>
where 
    T: Add<Output = T> + Sub<Output = T> + AddAssign + SubAssign + PartialOrd + 
       Default + Copy
{
    /// Creates a new Fenwick Tree for use in calculating and updating
    /// prefix sums.
    ///
    pub fn new(size: usize) -> Self
    {
        // Ensure size is 1 plus a power of 2.
        let n_bits = (size as f64).log(2_f64).ceil();
        let size   = 2_usize.pow(n_bits as u32) + 1 as usize;
        
        Fenwick { data: vec![T::default(); size], size }
    }
    
    /// Returns the sum of the first `idx` elements (indices 0 to `idx`)
    /// Equivalent to `.range_sum(0, idx)`. Range inclusive, [0..idx].
    ///
    pub fn prefix_sum(&self, idx: usize) -> T
    {
        let mut sum = self.data[0];
        let mut i   = idx;
        while i != 0 { 
            sum += self.data[i];
            i   -= lsb_usize!(i);
        }
        sum
    }
    
    /// Returns the index of the last element. This can be used as a parameter
    /// to `.range_sum()` or other methods.
    ///
    pub fn end(&self) -> usize
    {
        self.size - 1
    }
    
    /// Add `delta` to element with index `idx` (zero-based).
    ///
    pub fn add(&mut self, idx: usize, delta: T)
    {
        if idx == 0 {
            self.data[0] += delta;
        } else {
            let mut i = idx;
            while i < self.size {
                self.data[i] += delta;
                i += lsb_usize!(i);
            }
        }
    }
    
    /// Subtract `delta` from element with index `idx`.
    /// 
    pub fn sub(&mut self, idx: usize, delta: T)
    {
        if idx == 0 {
            self.data[0] -= delta;
        } else {
            let mut i = idx;
            while i < self.size {
                self.data[i] -= delta;
                i += lsb_usize!(i);
            }
        }
    }
    
    /// Set (as opposed to adjust) a single element's value.
    ///
    pub fn set(&mut self, idx: usize, value: T)
    {
        self.add(idx, value - self.get(idx))
    }
    
    /// Return a single element's value.
    ///
    pub fn get(&self, idx: usize) -> T
    {
        if idx > 0 {
            self.range_sum(idx - 1, idx)
        } else {
            self.data[0]
        }
    }
    
    /// Returns the sum of elements from `idx_i + 1` to `idx_j`, Equivalent to 
    /// `.prefix_sum(idx_j) - .prefix_sum(idx_i)`, but faster.
    ///
    pub fn range_sum(&self, idx_i: usize, idx_j: usize) -> T
    {
        let mut sum = T::default();
        let mut i   = idx_i;
        let mut j   = idx_j;
        while j > i {
            sum += self.data[j];
            j -= lsb_usize!(j);
        }
        while i > j {
            sum -= self.data[i];
            i -= lsb_usize!(i);
        }
        sum
    }
    
    /// Find the largest index with `.prefix_sum(index) <= value`.
    /// NOTE: Requires all values are non-negative.
    ///
    pub fn rank_query(&self, value: T) -> usize
    {
        let mut i = 0;
        let mut j = self.size - 1;
        let mut v = value - self.data[0];
        
        while j > 0 {
            if i + j < self.size && self.data[i + j] <= v {
                v -= self.data[i + j];
                i += j;
            }
            j >>= 1;
        }
        i
    }
    
    /// Find the minimum index with `.prefix_sum(index) >= value`.
    /// Also requires all values non-negative.
    ///
    pub fn min_rank_query(&self, value: T) -> usize
    {
        let mut idx = self.rank_query(value);
        if self.prefix_sum(idx) < value && idx < self.end() {
            idx += 1;
        } else {
            while self.get(idx) == T::default() && idx > 0 {
                idx -= 1;
            }
        }
        idx
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    
    #[test]
    fn new() {
        let fw = Fenwick::<i32>::new(8);
        assert_eq!(fw.end(), 8);
    }
    
    #[test]
    fn add() {
        let mut fw = Fenwick::new(8);
        fw.add(1, 1);
        fw.add(2, 1);
        fw.add(1, 1);
        assert_eq!(fw.prefix_sum(2), 3); // Range is inclusive.
    }
    
    #[test]
    fn sub() {
        let mut fw = Fenwick::new(8);
        fw.add(1, 1);
        fw.add(2, 1);
        fw.add(2, 1);
        fw.sub(2, 1);
        assert_eq!(fw.prefix_sum(2), 2);
    }
    
    #[test] 
    fn prefix_sum() {
        let mut fw = Fenwick::new(8);
        fw.add(1, 1);  // sum = 2
        fw.add(0, 1);  // sum = 1
        fw.add(3, 1);  // sum = 3
        assert_eq!(fw.prefix_sum(3), 3);
        
        fw.sub(0, 1);
        assert_eq!(fw.prefix_sum(3), 2);
    }
    
    #[test]
    fn get() {
        let mut fw = Fenwick::new(8);
        fw.add(0, 5);
        fw.add(1, 4);
        assert_eq!(fw.get(1), 4);
        assert_eq!(fw.get(0), 5);
    }
    
    #[test]
    fn set() {
        let mut fw = Fenwick::new(8);
        fw.add(0, 5);
        fw.add(1, 4);
        assert_eq!(fw.get(1), 4);
        assert_eq!(fw.get(0), 5);
        
        fw.set(1, 0);
        assert_eq!(fw.get(1), 0);
        assert_eq!(fw.prefix_sum(8), 5);
    }
    
    #[test]
    fn rank_query() {
        let mut fw = Fenwick::new(8);
        fw.add(0, 1);  // sum = 1
        fw.add(1, 1);  // sum = 2 
        fw.add(2, 3);  // sum = 5
        fw.add(3, 1);  // sum = 6
        fw.add(4, 1);  // sum = 7
        
        assert_eq!(fw.rank_query(5), 2);
        assert_eq!(fw.rank_query(6), 3);
        
        fw.set(0, 0);
        assert_eq!(fw.rank_query(5), 3);
    }
    
    #[test]
    fn min_rank_query_pre_test() {
        let mut fw = Fenwick::new(8);
        fw.add(0, 1);
        
        assert_eq!(fw.rank_query(2), 8); // Note the greedy nature of query.
        assert_eq!(fw.end(), 8);         // Goes all the way to the end.
        
        // fw.add(0, 1);  // sum = 1
        fw.add(1, 1);  // sum = 2 
        fw.add(2, 3);  // sum = 5
        fw.add(3, 1);  // sum = 6
        fw.add(4, 1);  // sum = 7
        
        // Stops shy of finding min index >= 3.
        assert_eq!(fw.rank_query(3), 1); 
        
        // If we want to know what index we'll land on to satisfy a value of
        // 3, we'll have to test the prefix sum and increment.
        let mut idx = fw.rank_query(3);
        
        assert_eq!(fw.prefix_sum(idx), 2);
        
        if fw.prefix_sum(idx) < 3 {
            idx += 1;
        }
        assert_eq!(idx, 2);
        
        // If we choose a value larger than the prefix sum, we should end up
        // at the tail end of the array and have to backtrack to find the
        // index of an element with a non-0 value.
        idx = fw.rank_query(8);
        
        assert_eq!(idx, 8);
        
        while fw.get(idx) == 0 {
            idx -= 1;
        }
        assert_eq!(idx, 4);
    }
    
    #[test]
    fn min_rank_query() {
        let mut fw = Fenwick::new(8);
        fw.add(0, 1);  // sum = 1
        fw.add(1, 1);  // sum = 2 
        fw.add(2, 3);  // sum = 5
        fw.add(3, 1);  // sum = 6
        fw.add(4, 1);  // sum = 7
        
        assert_eq!(fw.min_rank_query(3), 2);  // Check basic.
        assert_eq!(fw.min_rank_query(8), 4);  // Check backtracking.
        
        fw.add(7, 3);  // sum = 10
        
        assert_eq!(fw.min_rank_query(7), 4);  // Should backtrack to 4.
        assert_eq!(fw.min_rank_query(8), 7);  // Should advance to 7.
        
    }
}







