
//! A generic implementation for a Fenwick Tree, useful for managing prefix
//! sums. Each operation has `O(log n)` time-complexity. Some operations like
//! `.end()` have `O(1)` time-complexity.
//!
//! The code was originally taken from the Wikipedia article on Fenwick Trees 
//! and modified to create a class. The code has also been modified so the 
//! tree is 0-based despite being taken from a 1-based implementation 
//! originally.
//!
//! Wikipedia article: https://en.wikipedia.org/wiki/Fenwick_tree
//!

use std::ops::Add;
use std::ops::Sub;
use std::ops::AddAssign;
use std::ops::SubAssign;
use std::cmp::Ord;

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
    T: Add<Output = T> + Sub<Output = T> + AddAssign + SubAssign + Ord + 
       Default + Copy
{
    /// Creates a new Fenwick Tree for use in calculating and updating
    /// prefix sums.
    ///
    pub fn new(size: usize) -> Self
    {
        // Ensure size is 1 plus a power of 2.
        let n_bits = (size as f64).log(2_f64).ceil();
        let size   = 2_usize.pow(n_bits as u32) + 1_usize;
        
        Fenwick { data: vec![T::default(); size], size }
    }
    
    /// Creates a new Fenwick instance from the provided slice. The data in 
    /// the slice itself doesn't need to be in accumulated prefix sum form.
    /// It should just be a slice of unsummed values.
    ///
    pub fn from_slice(slice: &[T]) -> Self
    {
        // Ensure size is 1 plus a power of 2.
        let n_bits = (slice.len() as f64).log(2_f64).ceil();
        let size   = 2_usize.pow(n_bits as u32) + 1_usize;
        
        let mut data = Vec::with_capacity(size);
        
        data.extend_from_slice(slice);        
        data.resize(size, T::default());
        
        for i in 1..size {
            let j = i + lsb_usize!(i);
            if j < size {
                let d = data[i];
                data[j] += d;
            }
        }
        Fenwick { data, size }
    }
    
    /// Returns the sum of the first `idx` elements (indices 0 to `idx`)
    /// Equivalent to `.range_sum(0, idx)`. Range inclusive, [0..idx].
    ///
    pub fn prefix_sum(&self, idx: usize) -> T
    {
        debug_assert!(idx <= self.end());
        let mut sum = self.data[0];
        let mut i   = idx;
        while i != 0 { 
            sum += self.data[i];
            i   -= lsb_usize!(i);
        }
        sum
    }
    
    /// Returns the total prefix sum of all the elements.
    ///
    pub fn total(&self) -> T
    {
        self.data[self.end()] + self.data[0]
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
        debug_assert!(idx <= self.end());
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
        debug_assert!(idx <= self.end());
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
        debug_assert!(idx <= self.end());
        let cur_val = self.get(idx);
        if cur_val <= value {
            self.add(idx, value - cur_val);
        } else {
            self.sub(idx, cur_val - value);
        }
    }
    
    /// Return a single element's value.
    ///
    pub fn get(&self, idx: usize) -> T
    {
        debug_assert!(idx <= self.end());
        if idx == 0 {
            self.data[0]
        } else {
            self.range_sum(idx, idx)
        }
    }
    
    /// Returns the sum of elements from `idx_i` to `idx_j` inclusive, Similar 
    /// to `.prefix_sum(idx_j) - .prefix_sum(idx_i - 1)`, but faster.
    ///
    pub fn range_sum(&self, idx_i: usize, idx_j: usize) -> T
    {
        debug_assert!(idx_i <= idx_j && idx_j <= self.end());
        let mut sum = if idx_i > 0 { T::default() } else { self.data[0] };
        let mut i   = if idx_i > 0 { idx_i - 1    } else { 0            };
        let mut j   = idx_j;
        while j > i {
            sum += self.data[j];
            j   -= lsb_usize!(j);
        }
        while i > j {
            sum -= self.data[i];
            i   -= lsb_usize!(i);
        }
        sum
    }
    
    /// Find the largest index with `.prefix_sum(index) <= value`.
    /// NOTE: Requires all values are non-negative.
    ///
    pub fn rank_query(&self, value: T) -> usize
    {
        debug_assert!(self.data.iter().all(|&n| n >= /* 0 */ T::default()));
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
    
    /// Find the smallest index with `.prefix_sum(index) >= value` - if there is
    /// an index where the prefix sum is >= value; however, if not the case,
    /// this method will return the index of the last element with a non-0 
    /// value.
    /// NOTE: This also requires all values non-negative.
    ///
    pub fn min_rank_query(&self, value: T) -> usize
    {
        debug_assert!(self.data.iter().all(|&n| n >= /* 0 */ T::default()));
        if value <= self.data[0] {
            0
        } else {
            let mut i = self.end();
            let mut d = self.data[i];
            let mut v = (value - self.data[0]).min(d);
            
            while i & 0x01 == 0 {
                if d < v {
                    v -= d;
                    i += lsb_usize!(i >> 1);
                } else {
                    i -= lsb_usize!(i >> 1);
                }
                d = self.data[i];
            }
            if v > d {
                i + 1
            } else {
                i
            }
        }
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
        
        assert_eq!(fw.prefix_sum(fw.end()), 3);
        fw.add(fw.end(), 3);
        assert_eq!(fw.prefix_sum(fw.end()), 6);
        assert_eq!(fw.prefix_sum(fw.end() - 1), 3);
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
        assert_eq!(fw.get(fw.end()), 0);
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
        
        assert_eq!(fw.get(0), 5);
        fw.set(0, 4);
        assert_eq!(fw.get(0), 4);
        
        fw.set(fw.end(), 3);
        assert_eq!(fw.get(fw.end()), 3);
        assert_eq!(fw.prefix_sum(fw.end()), 7);
        assert_eq!(fw.total(), 7);
        
        fw.set(0, 8);
        assert_eq!(fw.get(0), 8);
        fw.set(0, 0);
        assert_eq!(fw.get(0), 0);
    }
    
    #[test]
    fn range_sum() {
        let mut fw = Fenwick::new(8);
        fw.add(0, 1);  // sum = 1
        fw.add(1, 1);  // sum = 2 
        fw.add(2, 3);  // sum = 5
        fw.add(3, 1);  // sum = 6
        fw.add(4, 1);  // sum = 7
        
        assert_eq!(fw.range_sum(1, 3), 5);
        assert_eq!(fw.range_sum(0, 3), 6);
        
        assert_eq!(fw.range_sum(2, 4), 5);
        assert_eq!(fw.range_sum(0, 0), 1);
    }
    
    #[test]
    fn rank_query() {
        let mut fw = Fenwick::new(8);
        fw.add(0, 1);  // sum = 1
        fw.add(1, 1);  // sum = 2 
        fw.add(2, 3);  // sum = 5
        fw.add(3, 1);  // sum = 6
        fw.add(4, 1);  // sum = 7

        assert_eq!(fw.rank_query(0), 0);
        assert_eq!(fw.rank_query(1), 0);
        assert_eq!(fw.rank_query(2), 1);
        assert_eq!(fw.rank_query(7), 8);
        
        assert_eq!(fw.rank_query(5), 2);
        assert_eq!(fw.rank_query(6), 3);
        
        fw.set(0, 0);
        assert_eq!(fw.rank_query(5), 3);
    }
    
    #[test]
    fn min_rank_query() {
        let mut fw = Fenwick::new(8);
        fw.add(0, 1);  // sum = 1
        fw.add(1, 1);  // sum = 2 
        fw.add(2, 3);  // sum = 5
        fw.add(3, 1);  // sum = 6
        fw.add(4, 1);  // sum = 7

        assert_eq!(fw.min_rank_query(1), 0);
        
        assert_eq!(fw.min_rank_query(3), 2);  // Check basic.
        assert_eq!(fw.min_rank_query(8), 4);  // Check that it falls on min.
        
        fw.add(7, 3);  // sum = 10
        
        assert_eq!(fw.min_rank_query(7), 4);  // Should fall to 4.
        assert_eq!(fw.min_rank_query(8), 7);  // Should advance to 7.

        assert_eq!(fw.min_rank_query(10), 7);
        assert_eq!(fw.min_rank_query(11), 7);
    }
}







