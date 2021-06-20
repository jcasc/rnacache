/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 André Müller (muellan@uni-mainz.de)
 *                       & Robin Kobus  (kobus@uni-mainz.de)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *****************************************************************************/

#ifndef MC_CMDLINE_INTERFACE_H_
#define MC_CMDLINE_INTERFACE_H_


#include <string>
#include <vector>
#include <cstdint>
#include <cstddef>
#include <limits>
#include <string>
#include <thread>
#include <vector>

#include "config.h"
#include "taxonomy.h"
#include "cmdline_utility.h"
#include "sequence_io.h"
#include "io_options.h"


namespace mc {


/*************************************************************************//**
 *
 *
 *  S H A R E D
 *
 *
 *****************************************************************************/



/*************************************************************************//**
 * @brief sequence sketching parameters
 *****************************************************************************/
struct sketching_options
{
    int kmerlen = 16;

    // number of features (kmer hashes) in a sketch of ONE window
    int sketchlen = 16;

    // number of characters in one window
    int winlen = 46;


    // difference between two successive window start positions
    int winstride = -1;  // < 0 : automatic: winstride = (winlen - (kmerlen-1))
};




/*************************************************************************//**
 * @brief
 *****************************************************************************/
struct database_storage_options
{
    float maxLoadFactor = -1;  // < 0 : use database default

    // restrict number of locations per feature
    int maxLocationsPerFeature = -1;  // < 0: use database default
    bool removeOverpopulatedFeatures = true;


    // restrict number of taxa (on a given rank) per feature
    bool removeAmbigFeatures = false;
    int maxTaxaPerFeature = 1;
};





/*************************************************************************//**
 *
 *
 *  B U I L D   M O D E
 *
 *
 *****************************************************************************/

/*************************************************************************//**
 *
 * @brief database creation parameters
 *
 *****************************************************************************/
struct build_options
{
    std::string dbfile;
    std::vector<std::string> infiles;

    sketching_options sketching;
    database_storage_options dbconfig;

    info_level infoLevel = info_level::moderate;
};



/*************************************************************************//**
 * @brief command line args -> database creation parameters
 *****************************************************************************/
build_options get_build_options(const cmdline_args&, build_options
                                defaults = build_options{});


/*************************************************************************//**
 * @brief build mode documentation
 *****************************************************************************/
std::string build_mode_usage();
std::string build_mode_examples();
std::string build_mode_docs();





/*************************************************************************//**
 *
 *
 *  M O D I F Y   M O D E
 *
 *
 *****************************************************************************/
using modify_options = build_options;

/*************************************************************************//**
 * @brief command line args -> database creation parameters
 *****************************************************************************/
modify_options get_modify_options(const cmdline_args&, modify_options
                                  defaults = modify_options{});


/*************************************************************************//**
 * @brief modify mode documentation
 *****************************************************************************/
std::string modify_mode_usage();
std::string modify_mode_examples();
std::string modify_mode_docs();





/*************************************************************************//**
 *
 *
 *  Q U E R Y   M O D E
 *
 *
 *****************************************************************************/

/*************************************************************************//**
 * @brief pairing of queries
 *****************************************************************************/
enum class pairing_mode : unsigned char {
    none, files, sequences
};


/*************************************************************************//**
 * @brief alignment mode
 *****************************************************************************/
enum class align_mode : unsigned char {
    none, semi_global
};


/*************************************************************************//**
 * @brief how to show mapping
 *****************************************************************************/
enum class map_view_mode : unsigned char {
    none, mapped_only, all
};


/*************************************************************************//**
 *
 * @brief how taxon formatting will be done
 *
 *****************************************************************************/
struct taxon_print_style {
    bool showName = true;
    bool showId = false;
};


/*************************************************************************//**
 *
 * @brief how to process input queries
 *
 *****************************************************************************/
struct performance_tuning_options {
    int numThreads = std::thread::hardware_concurrency();
    std::size_t batchSize = 4096;
    //limits number of reads per sequence source (file)
    std::int_least64_t queryLimit = std::numeric_limits<std::int_least64_t>::max();

    #ifdef RC_BAM
    size_t bamBufSize = 25;
    int bamThreads = std::thread::hardware_concurrency();

    #endif
};


/*************************************************************************//**
 *
 * @brief classification options
 *
 *****************************************************************************/
enum class coverage_norm {
    none, max
};

enum class coverage_fill {
    matches, fill
};

struct classification_options
{

    //maximum range in sequence that read (pair) is expected to be in
    std::size_t insertSizeMax = 0;

    std::size_t maxNumCandidatesPerQuery = 2;

    std::uint16_t hitsMin = 4;
    double hitsCutoff = 0.8;
    double covMin = 0.9;
    coverage_norm covNorm = coverage_norm::max;
    coverage_fill covFill = coverage_fill::matches;
    int maxEditDist = -1;

};


/*************************************************************************//**
 *
 * @brief ground truth based testing
 *
 *****************************************************************************/
struct classification_evaluation_options
{
    //show ground thruth if available
    bool showGroundTruth = false;

    //test accuracy (ground truth must be available)
    bool accuracy = false;

    //show known taxon (or complete lineage if 'showLineage' on)
    bool determineGroundTruth = false;
};


/*************************************************************************//**
 *
 * @brief   tokens that separate mapping output fields
 *
 * @details Note that the choice of separators, prefixes and suffixes
 *          reflects the output format of the first public MetaCache
 *          version.
 *
 *****************************************************************************/
struct formatting_tokens {
    //prefix for each non-mapping line
    std::string comment = "# ";
    std::string none = "--";
    //column separator
    std::string column = "\t|\t";
    //taxon separator
    std::string taxSeparator = ",";
    //if both taxid AND taxon name are to be printed,
    //taxids will be enclosed by these:
    std::string taxidPrefix = "(";
    std::string taxidSuffix = ")";
};


/*************************************************************************//**
 *
 * @brief classification output formatting
 *
 *****************************************************************************/
struct classification_output_formatting
{
    //how to show classification (read mappings), if 'none', only summary will be shown
    map_view_mode mapViewMode = map_view_mode::all;

    bool showQueryIds = false;

    taxon_print_style taxonStyle;

    formatting_tokens tokens;
};


/*************************************************************************//**
 *
 * @brief classification analysis options
 *
 *****************************************************************************/
struct classification_analysis_options
{
    //show top candidate sequences and their associated k-mer hash hit count
    bool showTopHits = true;

    //show all k-mer-hash hits in database for each given read
    bool showAllHits = false;
    //show candidate position(s) in reference sequence(s)
    bool showLocations = false;

};


/*************************************************************************//**
 *
 * @brief classification output options
 *
 *****************************************************************************/
struct classification_output_options
{
    classification_analysis_options analysis;
    classification_output_formatting format;
    classification_evaluation_options evaluate;

    //show classification summary
    bool showQueryParams = true;
    bool showSummary = true;
    bool showDBproperties = false;
    bool showErrors = true;
};


enum class rna_mode {
    map,
    sam,
    bam,
    sam_align,
    search
};

/*************************************************************************//**
 *
 * @brief query options
 *
 *****************************************************************************/
struct query_options
{
    std::string dbfile;
    std::vector<std::string> infiles;

    // how to pair up reads
    pairing_mode pairing = pairing_mode::none;

    // make a separate output file for each input file
    bool splitOutputPerInput = false;
    // output filename for mappings per read
    std::string queryMappingsFile;

    database_storage_options dbconfig;

    // query sketching options (all set to -1 : use value from database)
    sketching_options sketching {-1,-1,-1,-1};

    performance_tuning_options performance;

    classification_options classify;
    classification_output_options output;

    rna_mode rnaMode = rna_mode::map;

};


/*************************************************************************//**
 * @brief command line args -> query mode options
 *****************************************************************************/
query_options get_query_options(const cmdline_args&, query_options
                                defaults = query_options{});



/*************************************************************************//**
 * @brief query mode documentation
 *****************************************************************************/
std::string query_mode_usage();
std::string query_mode_examples();
std::string query_mode_docs();





/*************************************************************************//**
 *
 *
 *  I N F O   M O D E
 *
 *
 *****************************************************************************/
enum class info_mode {
    basic,
    targets,
    db_config, db_statistics, db_feature_map, db_feature_counts
};

struct info_options {
    std::string dbfile;
    info_mode mode = info_mode::basic;
    std::vector<std::string> targetIds;
};



/*************************************************************************//**
 * @brief command line args -> info mode options
 *****************************************************************************/
info_options get_info_options(const cmdline_args&);



/*************************************************************************//**
 * @brief info mode documentation
 *****************************************************************************/
std::string info_mode_usage();
std::string info_mode_examples();
std::string info_mode_docs();

} // namespace mc


#endif
