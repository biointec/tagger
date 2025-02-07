# Big-BWT

Tool to build the BWT and optionally the Suffix Array and LCP array for highly repetitive files using the approach described in *Prefix-Free Parsing for Building Big BWTs* by Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini [1].

Copyrights 2018- by the authors. 
 

## Prerequisites 


* GCC 4.8.1 or later

* Python version 3.8 or later

* [psutil](https://pypi.org/project/psutil/) (`pip install psutil`)



## Installation

* Download/Clone the repository
* `make` (create the C/C++ executables) 
* `bigbwt -h` (get usage instruction)
 

## Sample usage

The only requirement for the input file is that it does not contain the characters 0x00, 0x01, and 0x02 which are used internally by `bigbwt`. To build the BWT for file *yeast.fasta* just type

```bash
        bigbwt yeast.fasta
```
If no errors occurrs the BWT file yeast.fasta.bwt is created: it should be one character longer than the input; the extra character is the BWT eos symbol represented by the ASCII character 0x00. Two important command line options are the window size `-w` and the modulus `-m`. Their use is explained in the paper. 

Using option `-S` it is possible to compute the Suffix Array at the same time of the BWT. 
The Suffix Array is written in a file with extension `.sa` using 5 bytes per integer (5 is a compilation constant defined in utils.h). This is exactly the same format used by the [pSAscan/SAscan](https://www.cs.helsinki.fi/group/pads/) tools that indeed should produce exactly the same Suffix Array file. 

Using options `-s` and/or `-e` (in alternative to `-S`) it is possible to compute a sampled Suffix Array. Option `-s` generates a `.ssa` file containing the SA entries at the beginning of BWT runs. More precisely, the `.ssa` file contains pairs of 5 bytes integers *<j,SA[j]>* for all indices *j* such that *BWT[j]!=BWT[j-1]*. Option `-e` generates a `.esa` file containing the SA entries at the end of BWT runs: i.e. the pairs *<j,SA[j]>* for all indices *j* such that *BWT[j]!=BWT[j+1]*.

Using opton `-L` together with `-S` the LCP array is written to a file with extension `.lcp` using again 5 bytes per integer. The computation of the LCP array is done using the  [EM-SparsePhi 0.2.0 algorithm](https://www.cs.helsinki.fi/group/pads/better_em_laca.html) [2] which is an external memory algorithm and therefore works even when the input is larger than the available RAM. However, this algorithm does not take advantage of the repetitiveness in the input text: an implementation of the computation of the LCP values within the Prefix-Free Parsing algorithm is currently being developed.

The `bigbwt` tool has some support for multiple threads. Use option `-t` to specify the number of helper threads: in our tests `-t 4` reduced the running time by roughly a factor two.


You can check the correctness of the computation using the `-C` option. With this option, after the usual computation, `bigbwt` computes the SA using [divsufsort](https://github.com/y-256/libdivsufsort) and, if necessary, the LCP array using LCP9 [4]. These arrays are then used to check the correctness of the SA and LCP arrays previously computed using Prefix-Free Parsing (at the moment we do not support to check the correctness of partial SAs). Please note that divsufsort and LCP9 do not take advantage of the repetitiveness in the input text, so they are relatively slow and use 9n bytes of RAM for checking the BWT and SA, and 13n bytes of RAM when checking also the LCP array. Note that the checking of `bigbwt` output can be done also at a later time by invoking the tool `bwtcheck` as follows:

```bash
        bwtcheck textfile bwtfile [-S safile] [-L lcpfile] [-v]
```

As a (deprecated) alternative to the above checking procedure, you can run `bigbwt` with option `-c`. 
This will compute the BWT using the [SACAK algorithm](https://github.com/felipelouza/gsa-is) [3] and compare it with the one computed using Prefix-Free parsing. Be warned that SACAK, although being the most space economical among linear time algorithms, still needs up to *9n* bytes of RAM, since it first compute the Suffix Array and then outputs the BWT (with extension .Bwt).



## References

\[1\]  Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini 
*Prefix-Free Parsing for Building Big BWTs* [Proc. WABI '18](http://drops.dagstuhl.de/opus/volltexte/2018/9304/).

\[2\] Juha Kärkkäinen, Dominik Kempa. *LCP Array Construction in External Memory*,
ACM Journal of Experimental Algorithmics (JEA), [Volume 21(1), 2016, Article 1.7](http://doi.acm.org/10.1145/2851491)

\[3\] Ge Nong, 
*Practical linear-time O(1)-workspace suffix sorting for constant alphabets*, ACM Trans. Inform. Syst., [Vol. 31, no. 3, pp. 1-15, 2013](https://dl.acm.org/doi/10.1145/2493175.2493180)

\[4\] Giovanni Manzini,
*Two Space Saving Tricks for Linear Time LCP Array Computation*, [Proc. SWAT '04](https://link.springer.com/chapter/10.1007/978-3-540-27810-8_32)


