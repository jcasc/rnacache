RNACache Vocabulary
=====================
```
 term        explanation

 query       read or pair of reads to be classified
 target      (transcriptomic) reference sequence

 window      sequence of w consecutive characters used for constructing a sketch
 kmer        sequence of k < w consecutive characters within a window
 feature     hash of a kmer
 sketch      set of s <= w-k+1 features
 sketcher    takes a window and returns a sketch

 location    {target id , window index within target}

 database    feature->location hash table + some auxiliary data
```
