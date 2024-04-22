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