## Minimizer optimizations

- [ ] Uses optimal DNA data structure, 2 bits per base
- [ ] Mby use multithreading
- [ ] Check for canonical kmer during creation of complement, not after
- [ ] In situation that we have a tie between two canonical kmers in one windows, we pick the one with smaller position. This will prodece different minimizers for reversed window. Chance for this is negligible tho, in context of whole process.
  - Test case should be added for this behavior.
- [ ] Reduce memory usage by not creating new kmers/rc_kmers for each windows, just remove first one and add last one
- [ ] Improve hash function
- [ ] Add Deque to not rerun all kmers