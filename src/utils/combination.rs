pub struct CombinationIterator {
    n: usize,
    i: usize,
    j: usize,
}

impl CombinationIterator {
    pub fn new(n: usize) -> Self {
        Self {
            n,
            i: 0,
            j: 0,
        }
    }
    pub fn len(&self) -> usize {
        self.n * self.n
    }
}

impl Iterator for CombinationIterator {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if self.i < self.n {
            let result = (self.i, self.j);
            self.j += 1;
            if self.j == self.n {
                self.i += 1;
                if self.i < self.n {
                    self.j = 0;
                    if self.i == self.j {
                        self.j += 1;
                    }
                }
            }
            Some(result)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_combination_iterator() {
        let mut comb_iter = CombinationIterator::new(3);
        assert_eq!(comb_iter.next(), Some((0, 0)));
        assert_eq!(comb_iter.next(), Some((0, 1)));
        assert_eq!(comb_iter.next(), Some((0, 2)));
        assert_eq!(comb_iter.next(), Some((1, 0)));
        assert_eq!(comb_iter.next(), Some((1, 1)));
        assert_eq!(comb_iter.next(), Some((1, 2)));
        assert_eq!(comb_iter.next(), Some((2, 0)));
        assert_eq!(comb_iter.next(), Some((2, 1)));
        assert_eq!(comb_iter.next(), Some((2, 2)));
        assert_eq!(comb_iter.next(), None);
    }
    #[test]
    fn test_long_combination_iterator() {
        let comb_iter = CombinationIterator::new(10);
        comb_iter.for_each(|(i, j)| {
            println!("{} {}", i, j);
        });
    }
}

#[derive(Hash, PartialEq, Eq)]
pub struct CombinationVecIterator {
    vec1: Vec<usize>,
    vec2: Vec<usize>,
    index1: usize,
    index2: usize,
}

impl CombinationVecIterator {
    pub fn new(vec1: Vec<usize>, vec2: Vec<usize>) -> Self {
        Self {
            vec1,
            vec2,
            index1: 0,
            index2: 0,
        }
    }
    pub fn is_empty(&self) -> bool {
        if self.vec1.is_empty() || self.vec2.is_empty() {
            true
        } else {
            false
        }
    }
}

impl Iterator for CombinationVecIterator {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_empty() {
            return None;
        }
        if self.index1 >= self.vec1.len() {
            return None;
        }
        if self.index2 >= self.vec2.len() {
            self.index1 += 1;
            self.index2 = 0;
        }
        if self.index1 >= self.vec1.len() {
            return None;
        }

        let result = Some((self.vec1[self.index1], self.vec2[self.index2]));
        self.index2 += 1;
        result
    }
}

// TODO: Make tests