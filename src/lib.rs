
//! A generic implementation for a Fenwick Tree, useful for managing prefix
//! sums. Each operation has `O(log n)` time-complexity. Some operations like
//! `.end()` have `O(1)` time-complexity.
//!
//! The code was originally taken from the Wikipedia article on Fenwick Trees 
//! and modified to create a class. The code has also been modified so the 
//! tree is 0-based despite being taken from a 1-based implementation 
//! originally.
//!
//! Wikipedia article: <https://en.wikipedia.org/wiki/Fenwick_tree>
//!

use std::iter::{FromIterator, IntoIterator};
use std::ops::Add;
use std::ops::Sub;
use std::ops::AddAssign;
use std::ops::SubAssign;
use std::cmp::Ord;

macro_rules! lsb {
    ($i:expr) => { $i & $i.wrapping_neg() }
}

fn msb(n: usize) -> usize {
    if n == 0 { 0 } else { 1 << n.ilog(2) }
}

/// Represents a prefix sum array with `O(log n)` update operations.
///
#[derive(Debug)]
pub struct Fenwick<T> {
    data         : Vec<T>,
    size         : usize,
    max_idx_msb  : usize,
}

impl<T> Fenwick<T>
where 
    T: Add<Output = T> + Sub<Output = T> + AddAssign + SubAssign + Ord + 
       Default + Copy, 
{
    /// Creates a new Fenwick Tree for use in calculating and updating
    /// prefix sums. The size is adjusted to be 1 + a power of 2 if it already
    /// isn't.
    ///
    pub fn new(size: usize) -> Self {
        Fenwick { 
            data        : vec![T::default(); size], 
            size,
            max_idx_msb : msb(size), 
        }
    }
    
    /// Creates a new Fenwick instance from the provided slice. The data in 
    /// the slice itself doesn't need to be in accumulated prefix sum form.
    /// It should just be a slice of unsummed values. This function has `O(n)`
    /// time-complexity.
    ///
    fn from_slice(slice: &[T]) -> Self {
        let mut data = slice.to_vec();
        
        for i in 1..=data.len() {
            let j = i + lsb!(i);
            if j <= data.len() {
                let d = data[i - 1];
                data[j - 1] += d;
            }
        }
        Fenwick { data, max_idx_msb: msb(slice.len()), size: slice.len() }
    }
    
    /// Creates a new Fenwick instance from the provided vector. The data in 
    /// the vector itself doesn't need to be in accumulated prefix sum form.
    /// It should just be a vector of unsummed values. This function has `O(n)`
    /// time-complexity. The vector passed in is incorporated directly into the
    /// tree without copying.
    ///
    fn from_vec(data: Vec<T>) -> Self {
        let max_idx_msb = msb(data.len());
        let len         = data.len();

        let mut data = data;
        
        for i in 1..=len {
            let j = i + lsb!(i);
            if j <= data.len() {
                let d = data[i - 1];
                data[j - 1] += d;
            }
        }
        Fenwick { data, max_idx_msb, size: len }
    }

    /// Returns a non-consuming iterator over the Fenwick Tree. The iterator 
    /// will return the prefix sum of each element in the tree. The iterator 
    /// iterates over elements with `O(log(n))` time-complexity each.
    /// 
    pub fn iter(&self) -> FenwickIter<T> {
        self.into_iter()
    }

    /// Returns the sum of the first `idx` elements (indices 0 to `idx`)
    /// Equivalent to `.range_sum(0, idx)`. Range inclusive, [0..idx].
    ///
    pub fn prefix_sum(&self, idx: usize) -> T {
        debug_assert!(idx < self.len());
        let mut idx = idx + 1;
        let mut sum = T::default();
        while idx > 0 {
            sum += self.data[idx - 1];
            idx -= lsb!(idx);
        }
        sum
    }
    
    /// Returns the total prefix sum of all the elements.
    ///
    pub fn total(&self) -> T {
        self.data[self.max_idx_msb - 1]
    }
    
    /// Returns the index of the last element. This can be used as a parameter
    /// to `.range_sum()` or other methods.
    ///
    pub fn len(&self) -> usize {
        self.size
    }
    
    /// Add `delta` to element with index `idx` (zero-based).
    ///
    pub fn add(&mut self, idx: usize, delta: T) {
        debug_assert!(idx < self.len());
        let mut idx = idx + 1;
        while idx <= self.data.len() {
            self.data[idx - 1] += delta;
            idx += lsb!(idx);
        }
    }
    
    /// Subtract `delta` from element with index `idx`.
    /// 
    pub fn sub(&mut self, idx: usize, delta: T) {
        debug_assert!(idx < self.len());
        let mut idx = idx + 1;
        while idx <= self.data.len() {
            self.data[idx - 1] -= delta;
            idx += lsb!(idx);
        }
    }
    
    /// Set (as opposed to adjust) a single element's value.
    ///
    pub fn set(&mut self, idx: usize, value: T) {
        debug_assert!(idx < self.len());
        let cur_val = self.get(idx);
        if cur_val <= value {
            self.add(idx, value - cur_val);
        } else {
            self.sub(idx, cur_val - value);
        }
    }
    
    /// Return a single element's value.
    ///
    pub fn get(&self, idx: usize) -> T {
        debug_assert!(idx < self.len());
        self.range_sum(idx, idx)
    }
    
    /// Returns the sum of elements from `idx_i` to `idx_j` inclusive, Similar 
    /// to `.prefix_sum(idx_j) - .prefix_sum(idx_i - 1)`, but faster.
    ///
    pub fn range_sum(&self, start: usize, end: usize) -> T {
        debug_assert!(start <= end && end < self.len());
        let mut sum = T::default();
        let mut i   = start;
        let mut j   = end + 1;
        while j > i {
            sum += self.data[j - 1];
            j   -= lsb!(j);
        }
        while i > j {
            sum -= self.data[i - 1];
            i   -= lsb!(i);
        }
        sum
    }

    /// Find the largest index with `.prefix_sum(index) <= value`.
    /// NOTE: Requires all values are non-negative.
    ///
    pub fn rank_query(&self, value: T) -> Option<usize> {
        debug_assert!(self.data.iter().all(|&n| n >= T::default()),
                      "All elements must be non-negative to use this feature.");
        let mut i = 0;
        let mut j = self.max_idx_msb;
        let mut v = value;
        
        while j > 0 {
            if i + j <= self.data.len() && self.data[i + j - 1] < v {
                v -= self.data[i + j - 1];
                i += j;
            }
            j >>= 1;
        }
        if i == 0 && self.data[0] > value { None } else { Some(i) }
    }
    
    /// Find the smallest index with `.prefix_sum(index) >= value` - if there is
    /// an index where the prefix sum is >= value; however, if not the case,
    /// this method will return the index of the last element with a non-0 
    /// value.
    /// NOTE: This also requires all values non-negative.
    ///
    pub fn min_rank_query(&self, value: T) -> Option<usize> {
        debug_assert!(self.data.iter().all(|&n| n >= T::default()), 
                      "All elements must be non-negative to use this feature.");
        let mut sum   = T::default();
        let mut index = 0;
        let mut step  = self.max_idx_msb;

        while step > 0 {
            let next_idx = index + step;
            
            if next_idx <= self.data.len() 
                && sum + self.data[next_idx - 1] < value
            {
                sum += self.data[next_idx - 1];
                index = next_idx;
            }
            step >>= 1;
        }
        if index < self.data.len() { Some(index) } else { None}
    }
}

impl<T> From<Vec<T>> for Fenwick<T>
where
    T: Add<Output = T> + Sub<Output = T> + AddAssign + SubAssign + Ord + 
       Default + Copy, 
{
    fn from(vec: Vec<T>) -> Self {
        Self::from_vec(vec)
    }
}

impl<T> From<&[T]> for Fenwick<T>
where
    T: Add<Output = T> + Sub<Output = T> + AddAssign + SubAssign + Ord + 
       Default + Copy, 
{
    fn from(slice: &[T]) -> Self {
        Self::from_slice(slice)
    }
}

impl<T> FromIterator<T> for Fenwick<T>
where
    T: Add<Output = T> + Sub<Output = T> + AddAssign + SubAssign + Ord + 
       Default + Copy, 
{
    fn from_iter<I>(iter: I) -> Self
    where
        I: IntoIterator<Item = T>,
    {
        Self::from_vec(iter.into_iter().collect::<Vec<T>>())
    }
}

pub struct FenwickIntoIter<T> {
    fw  : Fenwick<T>,
    idx : usize,
}

impl<T> Iterator for FenwickIntoIter<T>
where
    T: Add<Output = T> + Sub<Output = T> + AddAssign + SubAssign + Ord + 
       Default + Copy, 
{
    type Item = T;
    
    fn next(&mut self) -> Option<Self::Item> {
        if self.idx < self.fw.len() {
            self.idx += 1;
            Some(self.fw.get(self.idx - 1))
        } else {
            None
        }
    }
}

impl<T> IntoIterator for Fenwick<T>
where
    T: Add<Output = T> + Sub<Output = T> + AddAssign + SubAssign + Ord + 
       Default + Copy, 
{
    type Item = T;
    type IntoIter = FenwickIntoIter<T>;
    
    fn into_iter(self) -> Self::IntoIter {
        FenwickIntoIter { fw: self, idx: 0 }
    }
}

pub struct FenwickIter<'a, T> {
    fw  : &'a Fenwick<T>,
    idx : usize,
}

impl<'a, T> Iterator for FenwickIter<'a, T>
where
    T: Add<Output = T> + Sub<Output = T> + AddAssign + SubAssign + Ord + 
       Default + Copy, 
{
    type Item = T;
    
    fn next(&mut self) -> Option<Self::Item> {
        if self.idx < self.fw.len() {
            self.idx += 1;
            Some(self.fw.get(self.idx - 1))
        } else {
            None
        }
    }
}

impl<'a, T> IntoIterator for &'a Fenwick<T>
where
    T: Add<Output = T> + Sub<Output = T> + AddAssign + SubAssign + Ord + 
       Default + Copy, 
{
    type Item = T;
    type IntoIter = FenwickIter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        FenwickIter { fw: self, idx: 0 }
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    
    #[test]
    fn add() {
        let mut fw = Fenwick::new(8);
        fw.add(1, 1);
        fw.add(2, 1);
        fw.add(1, 1);
        assert_eq!(fw.prefix_sum(2), 3); // Range is inclusive.
        
        assert_eq!(fw.prefix_sum(fw.len() - 1), 3);
        fw.add(fw.len() - 1, 3);
        assert_eq!(fw.prefix_sum(fw.len() - 1), 6);
        assert_eq!(fw.prefix_sum(fw.len() - 2), 3);
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
        assert_eq!(fw.get(fw.len() - 1), 0);
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
        //assert_eq!(fw.prefix_sum(8), 5);
        
        assert_eq!(fw.get(0), 5);
        fw.set(0, 4);
        assert_eq!(fw.get(0), 4);

        fw.set(fw.len() - 1, 3);
        assert_eq!(fw.get(fw.len() - 1), 3);
        assert_eq!(fw.prefix_sum(fw.len() - 1), 7);
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

        assert_eq!(fw.total(), 7);
        
        assert_eq!(fw.range_sum(1, 3), 5);
        assert_eq!(fw.range_sum(0, 3), 6);
        
        assert_eq!(fw.range_sum(2, 4), 5);
        assert_eq!(fw.range_sum(0, 0), 1);
        assert_eq!(fw.range_sum(2, 2), 3);
    }

    #[test]
    fn rank_query() {
        let mut fw = Fenwick::new(8);
        fw.add(0, 1);  // sum = 1
        fw.add(1, 1);  // sum = 2 
        fw.add(2, 3);  // sum = 5
        fw.add(3, 1);  // sum = 6
        fw.add(4, 1);  // sum = 7

        assert_eq!(fw.rank_query(0), None);
        assert_eq!(fw.rank_query(1), Some(0));
        assert_eq!(fw.rank_query(2), Some(1));
        assert_eq!(fw.rank_query(7), Some(4));
        
        assert_eq!(fw.rank_query(5), Some(2));
        assert_eq!(fw.rank_query(6), Some(3));
        
        fw.set(0, 0);
        assert_eq!(fw.rank_query(5), Some(3));
    }
    
    #[test]
    fn min_rank_query() {
        let mut fw = Fenwick::new(8);
        fw.add(0, 1);  // sum = 1
        fw.add(1, 1);  // sum = 2 
        fw.add(2, 3);  // sum = 5
        fw.add(3, 1);  // sum = 6
        fw.add(4, 1);  // sum = 7

        assert_eq!(fw.min_rank_query(1), Some(0));
        
        assert_eq!(fw.min_rank_query(3), Some(2));  // Check basic.
        assert_eq!(fw.min_rank_query(8), None);  // Check that it falls on min.
        
        fw.add(7, 3);  // sum = 10
        
        assert_eq!(fw.min_rank_query(7), Some(4));  // Should fall to 4.
        assert_eq!(fw.min_rank_query(8), Some(7));  // Should advance to 7.

        assert_eq!(fw.min_rank_query(10), Some(7));
        assert_eq!(fw.min_rank_query(11), None);
    }

    #[test]
    fn iterators() {
        let mut fw = Fenwick::new(8);
        fw.add(0, 1);  // sum = 1
        fw.add(1, 1);  // sum = 2 
        fw.add(2, 3);  // sum = 5
        fw.add(3, 1);  // sum = 6
        fw.add(4, 1);  // sum = 7
        
        let mut iter = fw.iter();
        assert_eq!(iter.next(), Some(1));
        assert_eq!(iter.next(), Some(1));
        assert_eq!(iter.next(), Some(3));
        assert_eq!(iter.next(), Some(1));
        assert_eq!(iter.next(), Some(1));
        assert_eq!(iter.next(), Some(0));
        assert_eq!(iter.next(), Some(0));
        assert_eq!(iter.next(), Some(0));
        assert_eq!(iter.next(), None);
        
        let mut iter = fw.into_iter();
        assert_eq!(iter.next(), Some(1));
        assert_eq!(iter.next(), Some(1));
        assert_eq!(iter.next(), Some(3));
        assert_eq!(iter.next(), Some(1));
        assert_eq!(iter.next(), Some(1));
        assert_eq!(iter.next(), Some(0));
        assert_eq!(iter.next(), Some(0));
        assert_eq!(iter.next(), Some(0));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn from_iterator() {
        let fw = Fenwick::from_iter(vec![1, 1, 3, 1, 1, 0, 0, 0]);
        assert_eq!(fw.prefix_sum(3), 6);
        assert_eq!(fw.prefix_sum(4), 7);
        assert_eq!(fw.prefix_sum(5), 7);
        assert_eq!(fw.prefix_sum(6), 7);
        assert_eq!(fw.prefix_sum(7), 7);

        let fw = vec![1, 1, 3, 1, 1, 0, 0, 0, 0]
                    .into_iter().collect::<Fenwick<i32>>();
        assert_eq!(fw.prefix_sum(3), 6);
        assert_eq!(fw.prefix_sum(4), 7);
        assert_eq!(fw.prefix_sum(5), 7);
        assert_eq!(fw.prefix_sum(6), 7);
        assert_eq!(fw.prefix_sum(7), 7);

        let fw: Fenwick<_> = vec![1, 1, 3, 1, 1, 0, 0, 0, 0].into();
        assert_eq!(fw.prefix_sum(3), 6);
        assert_eq!(fw.prefix_sum(4), 7);
        assert_eq!(fw.prefix_sum(5), 7);
        assert_eq!(fw.prefix_sum(6), 7);
        assert_eq!(fw.prefix_sum(7), 7);

    }
}
