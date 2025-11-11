## Complex
This portion thus far explores the generation of Butson-type complex Hadamard matrices, along with their more specific relatives, permutation-core Butson-type complex Hadamard matrices.

The `.sage` program is the current iteration of my work, and the `.ipynb` is severely out of date.
The sage program is a bit harder to read, since a good portion of the logic is parallelized and therefore a bit less straightforward to read.

### Output

This contains the outputs of the program. It works essentially as intended, but generates more matrices than are given by the aalto.fi database. For example, butson-6-6.txt contains 6 CHMs, but the database suggests there should be only four, up to monomial equivalence. This is, however, a drastic improvement over the previous 20 that were generated before significant pruning measures were implemented.

## Real
At present, this is more of a shot-in-the-dark hobby project than anything deserving substantial recognition. `minimal-real-generator.py` is, in theory, an attempt to generate the lexicographically-minimal real Hadamard matrix in a given dimension.

The logic seems to be sound, but the combinatorial explosion seems to occur for a 4n x 4n real Hadamard matrix for n values greater than or equal to 5. Until and unless I come up with some flash of inspiration to radically optimize the program, it's essentially a moot project.

## TO DO
* Circulant Permutation Core: Program the routine to iterate through the set of permutations, apply the individual chosen permutation repeatedly, and see if the result is a Hadamard matrix.
* Efficiency Sorting: Parallelize the chunking based on more efficient computations. Since the earlier chosen rows have to be compared to everything in the list, whereas the last few rows only need to be compared to those after themselves, we could make the first chunk much smaller and parse them based on, say, the product between the number of rows in the chunk AND the number of rows exceeding that particular chunk. This may help the program run in parallel for longer.
* Equivalence Pruning: Write a separate function to traverse a list of complex Hadamard matrices and see if any of them are equivalent to others on the list. For example, since they're in ascending order based on sorted first row, such an algorithm could renormalize each matrix at all possible points, and then see if there exists a row in the new matrix that would have canonically appeared earlier in the list, thereby showing its equivalence to a prior matrix. This might not necessarily get rid of ALL duplicates, though I technically don't have a proof that it won't. It should also be applied separately from the main loop, both so that we can avoid excess computations slowing down at runtime, and because some of the potentially-redundant matrices are interesting, i.e. they may turn out to have symmetries not displayed in their minimal form.
* Further Optimization: Investigate how to further reduce computational overhead and RAM consumption so that we can run this algorithm for larger values of n. Since the RAM consumption seems to be primarily tied to the `generate_L1_vectors` function, maybe there's a way to more efficiently handle these vectors.