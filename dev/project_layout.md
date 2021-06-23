Configuration
---------------------------------------------------------------------
```
config.h                     all basic datatype definitions

version.h                    database and project version constants;
                             needed for database compatibility check
```

Mode Starting
---------------------------------------------------------------------
```
main.cpp                     main entry point; selects modes

modes.h                      declares mode starting functions
mode_query.cpp               start classification 
mode_build.cpp               database building 
mode_info.cpp                database property queries
mode_help.cpp                shows help files from /docs
```

Classification ("query" mode)
---------------------------------------------------------------------
```
query_options.h              all query/classification settings 
query_options.cpp            command line args -> query settings

classification.h             classification starting function
classify_common.h            query database and classify reads based on matches 
classify_rna.cpp             query database and classify reads based on matches 
candidates.h                 determine top hits in contiguous windows
matches_per_target.h         coverage map construction

querying.h                   multi-threaded, batched database query functions;
                             also paired-end read handling

printing.h/cpp               classification and analysis output functions
```

Database Operations / Sketching
---------------------------------------------------------------------
```                          
database.h            feature->location map + auxiliary data
hash_multimap.h              hashmap used by database
chunk_allocator.h            default allocator used by hash_multimap

bitmanip.h                   bit-level manipulation functions
dna_encoding.h               ASCII DNA strings -> 2-bit encoded kmers
hash_dna.h                   sketcher classes
hash_int.h                   integer hash functions
hash_family.h                family of hash functions; not used by default
                             (could probably be removed in the near future)

sequence_io.h/cpp            reference sequence readers for FASTA/FASTQ files
sequence_view.h              non-owning string view class
```

Analysis
---------------------------------------------------------------------
```
classification_statistics.h  classification statistics

alignment.h                  construct (semi-global) alignments

stat_confusion.h             thread-safe (binary) confusion statistics
stat_moments.h               accumulators for statistical moments
                             (mean, variance, skewness) and extrema (min/max)
stat_combined.h              combined statistical accumulators
```

Utilities
---------------------------------------------------------------------
```
args_handling.h/cpp          command line args parsing needed for all modes
args_parser.h                simple, ad-hoc command line arguments parser

filesys_utility.h/cpp
cmdline_utility.h/cpp

io_serialize.h               fundamental datatype (de-)serialization functions
io_options.h                 I/O related option types
io_error.h                   I/O exception definitions

timer.h                      simple std::chrono based timer class
string_utils.h               string processing utilities (trim)
typename.h                   name demangling functions
```
