fast_align_ng
==========

`fast_align_ng` is a fast, unsupervised word aligner which can use POS tags to improve performance. It is based on [`fast_align`](https://github.com/clab/fast_align),
and incorporates some additional features described in:

* [Simple extensions for a reparameterised IBM Model 2](http://people.eng.unimelb.edu.au/tcohn/papers/gelling14acl.pdf). Douwe Gelling and Trevor Cohn. In *Proceedings of ACL Short papers*, 2014.

Please cite this paper if you use this software.

# Compiling fast_align_ng

fast_align_ng depends on the TCLAP library and a c++11 compiler. You can run make in the base folder to compile, which should give you
an executable.

# Using fast_align_ng

The following command will give you forward alignments, add the `-r` flag for backwards alignments:

    ./fast_align_ng -i data.zh-en -P data.pos.zh-en > forward.align

The input file (`-i`) should have 2 sentences per line, separated by `|||`, which are translations of each other. Sentences should be tokenized and tokens separated by whitespace.
The POS file (`-P`) is optional, but should have the same format as the input file, and each line should have as many tokens (POS tags) as the corresponding input line.

run `fast_align_ng` with no arguments for additional options

# License

The source code in this repository is provided under the terms of the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0.html).

