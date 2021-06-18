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

/*************************************************************************//**
 *
 * @file Facilities that generate option objects from command line arguments.
 * Argument parsing and doc generation uses the 'CLIPP' library.
 * If an option has multiple associated flag strings like, e.g.,
 * "-no-map" or "-nomap" the first one will appear in the documentation
 * and is considered to be the canonical flag while following flags are
 * there to preserve backwards compatibility with older versions of MetaCach.
 *
 *****************************************************************************/

#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <regex>

#include "options.h"
#include "filesys_utility.h"
#include "database.h"

#include "../dep/clipp.h"


namespace mc {

using std::size_t;
using std::vector;
using std::string;
using std::to_string;

using namespace std::string_literals;



/*************************************************************************//**
 *
 * @brief collects all command line interface error messages
 *
 *****************************************************************************/
class error_messages {
public:
    error_messages& operator += (const string& message) {
        messages_.push_back(message);
        return *this;
    }
    error_messages& operator += (string&& message) {
        messages_.push_back(std::move(message));
        return *this;
    }

    bool any() const noexcept   { return !messages_.empty(); }

    string str() const {
        string s;
        for(const auto msg : messages_) {
            if(!msg.empty()) s += msg + '\n';
        }
        return s;
    }

//    auto begin() const noexcept { return messages_.begin(); }
//    auto end()   const noexcept { return messages_.end();   }

private:
    vector<string> messages_;
};




/*****************************************************************************
 *
 *  H E L P E R S
 *
 *****************************************************************************/
string taxon_rank_names(const string& separator = ", ")
{
    string s;
    for(auto r = taxon_rank::Sequence; r < taxon_rank::Domain; ++r) {
        s += taxonomy::rank_name(r);
        s += separator;
    }
    s += taxonomy::rank_name(taxon_rank::Domain);
    return s;
}



//-------------------------------------------------------------------
/// @return database filename with extension
string sanitize_database_name(string name)
{
    if(name.find(".db") == string::npos) {
        name += ".db";
    }
    return name;
}


//-------------------------------------------------------------------
/// @return replaces '\t' with tab char and remove other special chars
string sanitize_special_chars(const string& text)
{
    return std::regex_replace( std::regex_replace(text,
        // no newlines, vertical tabs, etc. allowed
        std::regex(R"((\\n)|(\\v)|(\\r)|(\\a))"), ""),
        // substitute literat "\t" with tab character
        std::regex(R"(\\t)"), "\t");
}


//-------------------------------------------------------------------
void replace_directories_with_contained_files(vector<string>& names)
{
    vector<string> result;
    result.reserve(names.size());

    for(const auto& name : names) {

        auto fnames = files_in_directory(name);

        if(!fnames.empty()) {
            result.insert(result.end(), std::make_move_iterator(fnames.begin()),
                                        std::make_move_iterator(fnames.end()));
        } else {
            result.push_back(name);
        }
    }

    std::swap(result, names);
}



//-------------------------------------------------------------------
auto cli_doc_formatting()
{
    return clipp::doc_formatting{}
        .first_column(0)
        .doc_column(22)
        .last_column(80)
        .indent_size(4)
        .line_spacing(1)
        .alternatives_min_split_size(2)
        .paragraph_spacing(2)
        .max_flags_per_param_in_usage(1)
        .max_flags_per_param_in_doc(1)
        ;
}

auto cli_usage_formatting()
{
    return cli_doc_formatting().first_column(4).line_spacing(0);
}



//-------------------------------------------------------------------
/// @brief adds 'parameter' that catches unknown args with '-' prefix
clipp::parameter
catch_unknown(error_messages& err) {
    return clipp::any(clipp::match::prefix{"-"},
        [&](const string& arg) { err += "unknown argument: "s + arg; });
}



//-------------------------------------------------------------------
taxon_rank rank_from_name(const string& name, error_messages& err)
{
    auto r = taxonomy::rank_from_name(name);

    if(r == taxon_rank::none) {
        err += "Unknown taxonomic rank '"s + name + "'!\n";
        err += "Valid rank names are:\n    " + taxon_rank_names("\n    ") + "\n";
    }

    return r;
}



//-------------------------------------------------------------------
void raise_default_error(const error_messages& err,
                         const string& mode = "",
                         const string& usage = "",
                         const string& examples = "")
{
    auto msg = err.str();

    if(!msg.empty())      msg += "\n";

    if(!usage.empty())    msg += "USAGE:\n" + usage + "\n\n";
    if(!examples.empty()) msg += "EXAMPLES:\n" + examples + "\n\n";

    if(!mode.empty()) {
        msg += "\nYou can view the full interface documentation of mode '"s
            + mode + "' with:\n    rnacache help " + mode + " | less";
    }

    throw std::invalid_argument{std::move(msg)};
}




/*****************************************************************************
 *
 *
 *  S H A R E D
 *
 *
 *****************************************************************************/

/// @brief
clipp::parameter
database_parameter(string& filename, error_messages& err)
{
    using namespace clipp;

    return value(match::prefix_not{"-"}, "database")
        .call([&](const string& arg){ filename = sanitize_database_name(arg); })
        .if_missing([&]{ err += "Database filename is missing!"; })
        .doc("database file name;\n"
             "An RNACache database contains "
             "min-hash signatures of reference sequences.\n");
}



//-------------------------------------------------------------------
clipp::group
info_level_cli(info_level& lvl, error_messages& err)
{
    using namespace clipp;

    return one_of (
        option("-silent").set(lvl, info_level::silent),
        option("-verbose").set(lvl, info_level::verbose)
        .if_conflicted([&]{
            err += "Info level must be either '-silent' or '-verbose'!";
        })
    )
        % "information level during build:\n"
          "silent => none / verbose => most detailed\n"
          "default: neither => only errors/important info";
}



//-------------------------------------------------------------------
/// @brief shared command-line options for sequence sketching
clipp::group
sketching_options_cli(sketching_options& opt, error_messages& err)
{
    using namespace clipp;
    return (
    (   option("-kmerlen") &
        integer("k", opt.kmerlen)
            .if_missing([&]{ err += "Number missing after '-kmerlen'!"; })
    )
        %("number of nucleotides/characters in a k-mer\n"
          "default: "s + (opt.kmerlen > 0 ? to_string(opt.kmerlen)
                                          : "determined by database"s))
    ,
    (   option("-sketchlen") &
        integer("s", opt.sketchlen)
            .if_missing([&]{ err += "Number missing after '-sketchlen'!"; })
    )
        %("number of features (k-mer hashes) per sampling window\n"
          "default: "s + (opt.sketchlen > 0 ? to_string(opt.sketchlen)
                                            : "determined by database"s))
    ,
    (   option("-winlen") &
        integer("w", opt.winlen)
            .if_missing([&]{ err += "Number missing after '-winlen'!"; })
    )
        %("number of letters in each sampling window\n"
          "default: "s + (opt.winlen > 0 ? to_string(opt.winlen)
                                         : "determined by database"s))
    ,
    (   option("-winstride") &
        integer("l", opt.winstride)
            .if_missing([&]{ err += "Number missing after '-winstride'!"; })
    )
        %("distance between window starting positions\n"
          "default: "s +
          (opt.winlen > 0 && opt.kmerlen > 0
              ? (to_string(opt.winlen - opt.kmerlen + 1) + " (w-k+1)")
              : "determined by database"s)
    )
    );
}



//-------------------------------------------------------------------
/// @brief shared command-line options for sequence sketching
clipp::group
database_storage_options_cli(database_storage_options& opt, error_messages& err)
{
    using namespace clipp;

    const database defaultDb;

    return (
    (   option("-max-locations-per-feature") &
        integer("#", opt.maxLocationsPerFeature)
            .if_missing([&]{ err += "Number missing after '-max-locations-per-feature'!"; })
    )
        %("maximum number of reference sequence locations to be stored per feature;\n"
          "If the value is too high it will significantly impact querying speed. "
          "Note that an upper hard limit is always imposed by the data type "
          "used for the hash table bucket size (set with "
          "compilation macro '-DMC_LOCATION_LIST_SIZE_TYPE')."
          "Can also be set in query mode.\n"
          "default: "s + to_string(defaultDb.max_locations_per_feature()))
    ,
    (
        option("-remove-overpopulated-features")
            .set(opt.removeOverpopulatedFeatures)
        %("Removes all features that have reached the maximum allowed "
          "amount of locations per feature. This can improve querying "
          "speed and can be used to remove non-discriminative features."
          "Can also be set in query mode.\n"
          "default: "s + (opt.removeOverpopulatedFeatures ? "on" : "off"))
    )
    ,
    (   option("-max-ambig-per-feature").set(opt.removeAmbigFeaturesOnRank, taxon_rank::Sequence) &
        integer("#", opt.maxTaxaPerFeature)
            .if_missing([&]{ err += "Number missing after '-max-ambig-per-feature'!"; })
    )
        % ("Maximum number of allowed different reference sequences per feature. "
           "Removes all features exceeding this limit from database.\n"
           "default: "s + (opt.removeAmbigFeaturesOnRank != taxon_rank::none ? to_string(opt.maxTaxaPerFeature) : "off"))
    ,
    (   option("-max-load-fac", "-max-load-factor") &
        number("factor", opt.maxLoadFactor)
            .if_missing([&]{ err += "Number missing after '-max-load-fac'!"; })
    )
        %("maximum hash table load factor;\n"
          "This can be used to trade off larger memory consumption for "
          "speed and vice versa. A lower load factor will improve speed, "
          "a larger one will improve memory efficiency.\n"
          "default: "s + to_string(defaultDb.max_load_factor())
    )
    );
}



/*****************************************************************************
 *
 *
 *  B U I L D   M O D E
 *
 *
 *****************************************************************************/
/// @brief build mode command-line options
clipp::group
build_mode_cli(build_options& opt, error_messages& err)
{
    using namespace clipp;

    return (
    "REQUIRED PARAMETERS" %
    (
        database_parameter(opt.dbfile, err)
        ,
        values(match::prefix_not{"-"}, "sequence file/directory", opt.infiles)
            .if_missing([&]{
                err += "No reference sequence files provided or found!";
            })
            % "FASTA or FASTQ files containing transcriptomic sequences.\n"
              "If directory names are given, they will be searched for "
              "sequence files (at most 10 levels deep).\n"
    ),
    "BASIC OPTIONS" %
    (
        info_level_cli(opt.infoLevel, err)
    ),
    "SKETCHING (SUBSAMPLING)" %
        sketching_options_cli(opt.sketching, err)
    ,
    "ADVANCED OPTIONS" %
    (
        database_storage_options_cli(opt.dbconfig, err)
    ),
    catch_unknown(err)
    );

}



//-------------------------------------------------------------------
build_options
get_build_options(const cmdline_args& args, build_options opt)
{
    error_messages err;

    auto cli = build_mode_cli(opt, err);

    auto result = clipp::parse(args, cli);

    if(!result || err.any()) {
        raise_default_error(err, "build", build_mode_usage());
    }

    replace_directories_with_contained_files(opt.infiles);

    if(opt.dbconfig.maxLocationsPerFeature < 0)
        opt.dbconfig.maxLocationsPerFeature = database::max_supported_locations_per_feature();

    auto& sk = opt.sketching;
    if(sk.winstride < 0) sk.winstride = sk.winlen - sk.kmerlen + 1;

    return opt;
}



//-------------------------------------------------------------------
string build_mode_usage() {
    return
    "    rnacache build <database> <sequence file/directory>... [OPTION]...\n\n"
    "    rnacache build <database> [OPTION]... <sequence file/directory>...";
}



//-------------------------------------------------------------------
string build_mode_examples() {
    return
    "    Build database 'mydb' from sequence file 'transcripts.fa':\n"
    "        rnacache build mydb transcripts.fa\n"
    "\n"
    "    Build database 'mydb' from two sequence files:\n"
    "        rnacache build mydb one.fa two.fa\n"
    "\n"
    "    Build database 'myTranscriptome' from folder containing sequence files:\n"
    "        rnacache build myTranscriptome all_transcripts\n";
}



//-------------------------------------------------------------------
string build_mode_docs() {

    build_options opt;
    error_messages err;

    auto cli = build_mode_cli(opt, err);

    string docs = "SYNOPSIS\n\n";

    docs += build_mode_usage();

    docs += "\n\n\n"
        "DESCRIPTION\n"
        "\n"
        "    Create a new database of reference sequences.\n"
        "\n\n";

    docs += clipp::documentation(cli, cli_doc_formatting()).str();

    docs += "\n\nEXAMPLES\n\n";
    docs += build_mode_examples();

    return docs;
}





/*****************************************************************************
 *
 *
 *  M O D I F Y   M O D E
 *
 *
 *****************************************************************************/

/// @brief build mode command-line options
clipp::group
modify_mode_cli(build_options& opt, error_messages& err)
{
    using namespace clipp;

    return (
    "REQUIRED PARAMETERS" %
    (
        database_parameter(opt.dbfile, err)
        ,
        values(match::prefix_not{"-"}, "sequence file/directory", opt.infiles)
            .if_missing([&]{
                err += "No reference sequence files provided or found!";
            })
            % "FASTA or FASTQ files containing transcriptomic sequences.\n"
              "If directory names are given, they will be searched for "
              "sequence files (at most 10 levels deep).\n"
    ),
    "BASIC OPTIONS" %
    (
        info_level_cli(opt.infoLevel, err)
    ),
    "ADVANCED OPTIONS" %
    (
        database_storage_options_cli(opt.dbconfig, err)
    ),
    catch_unknown(err)
    );

}



//-------------------------------------------------------------------
build_options
get_modify_options(const cmdline_args& args, modify_options opt)
{
    error_messages err;

    auto cli = modify_mode_cli(opt, err);

    auto result = clipp::parse(args, cli);

    if(!result || err.any()) {
        raise_default_error(err, "modify", modify_mode_usage());
    }

    // use settings from database file as defaults
    auto db = make_database(opt.dbfile, database::scope::metadata_only);

    const auto& ts = db.target_sketcher();
    auto& sk = opt.sketching;
    sk.kmerlen   = ts.kmer_size();
    sk.sketchlen = ts.sketch_size();
    sk.winlen    = ts.window_size();
    sk.winstride = ts.window_stride();

    opt.dbconfig.maxLoadFactor = db.max_load_factor();
    opt.dbconfig.maxLocationsPerFeature = db.max_locations_per_feature();

    // parse again
    clipp::parse(args, cli);

    replace_directories_with_contained_files(opt.infiles);

    if(opt.dbconfig.maxLocationsPerFeature < 0)
        opt.dbconfig.maxLocationsPerFeature = database::max_supported_locations_per_feature();

    return opt;
}



//-------------------------------------------------------------------
string modify_mode_usage() {
    return
    "    rnacache modify <database> <sequence file/directory>... [OPTION]...\n\n"
    "    rnacache modify <database> [OPTION]... <sequence file/directory>...";
}



//-------------------------------------------------------------------
string modify_mode_examples() {
    return
    "    Add reference sequence 'added.fa' to database 'transcripts'\n"
    "        rnacache modify transcripts added.fa\n";
}



//-------------------------------------------------------------------
string modify_mode_docs() {

    build_options opt;
    error_messages err;

    auto cli = modify_mode_cli(opt, err);

    string docs = "SYNOPSIS\n\n";

    docs += modify_mode_usage();

    docs += "\n\n\n"
        "DESCRIPTION\n"
        "\n"
        "    Add reference sequence to an existing database.\n"
        "\n\n";

    docs += clipp::documentation(cli, cli_doc_formatting()).str();

    docs += "\n\n\nEXAMPLES\n";
    docs += modify_mode_examples();

    return docs;
}





/*************************************************************************//**
 *
 *
 *  Q U E R Y   M O D E
 *
 *
 *****************************************************************************/
/// @brief command line interface for classification parameter tuning
clipp::group
classification_params_cli(classification_options& opt, error_messages& err)
{
    using namespace clipp;

    return (
    (   option("-hitmin", "-hit-min", "-hits-min", "-hitsmin") &
        integer("t", opt.hitsMin)
            .if_missing([&]{ err += "Number missing after '-hitmin'!"; })
    )
        %("Sets classification threshhold 't^min' to <t>.\n"
          "All candidates with fewer hits are discarded "
          "from the query's candidate set. Higher values will increase "
          "precision at the expense of recall.\n"
          "default: "s + to_string(opt.hitsMin))
    ,
    (   option("-maxcand", "-max-cand") &
        integer("#", opt.maxNumCandidatesPerQuery)
            .if_missing([&]{ err += "Number missing after '-maxcand'!"; })
    )
        %(std::is_same<classification_candidates, distinct_matches_in_contiguous_window_ranges>() ?
            "Has no effect. (Requires selection of best_distinct_matches_... candidate generator in config.h)." :
            "Maximum number of candidates to consider (before filtering!)."
            "default: "s + to_string(opt.maxNumCandidatesPerQuery))
    ,
    ( 
        option("-hit-cutoff", "-cutoff", "-hits-cutoff", "-hitcutoff", "-hitscutoff") &
        number("t", opt.hitsCutoff)
            .if_missing([&]{ err += "Number missing after '-hit-cutoff'!"; })
    )
        %("Sets classification threshhold 't^cutoff' to <t>.\n"
          "All candidates with fewer hits (relative to the maximal candidate) "
          "are discarded from the query's candidate set. Higher values will increase "
          "precision at the expense of recall.\n"
          "default: "s + to_string(opt.hitsCutoff))
    ,   
        option("-cov-min", "-covmin", "-coverage-min", "-coveragemin", "-coverage") &
        number("p", opt.covMin)
            .if_missing([&]{ err += "Number missing after '-cov-min'!"; })

    ,   
        option("-max-edit", "-max-edit-dist", "-max-edit-distance") &
        integer("t", opt.maxEditDist)
            .if_missing([&]{ err += "Number missing after '-cov-min'!"; })

    ,   
        option("-no-cov-norm", "-no-norm-coverage").set(opt.covNorm, coverage_norm::none)
    ,
        option("-fill-coverage", "-fill-in-coverage").set(opt.covFill, coverage_fill::fill)
    );
}



//-------------------------------------------------------------------
/// @brief build mode command-line options
clipp::group
classification_output_format_cli(classification_output_formatting& opt,
                                 error_messages& err)
{
    using namespace clipp;
    return (
        one_of(
            option("-no-map", "-nomap").set(opt.mapViewMode, map_view_mode::none)
            %("Don't report classification for each individual query "
              "sequence; show summaries only (useful for quick tests).\n"
              "default: "s + (opt.mapViewMode == map_view_mode::none ? "on" : "off"))
            ,
            option("-mapped-only", "-mappedonly").set(opt.mapViewMode, map_view_mode::mapped_only)
            %("Don't list unclassified reads/read pairs.\n"
              "default: "s + (opt.mapViewMode == map_view_mode::mapped_only ? "on" : "off"))
        )
        ,
        option("-taxids", "-taxid").set(opt.taxonStyle.showId)
            %("Print taxon ids in addition to taxon names.\n"
              "default: "s + (opt.taxonStyle.showId ? "on" : "off"))
        ,
        option("-taxids-only", "-taxidsonly")
            .set(opt.taxonStyle.showId).set(opt.taxonStyle.showName,false)
            %("Print taxon ids instead of taxon names.\n"
              "default: "s + (opt.taxonStyle.showId && !opt.taxonStyle.showName ? "on" : "off"))
        ,
        option("-omit-ranks", "-omitranks").set(opt.taxonStyle.showRankName,false)
            %("Do not print taxon rank names.\n"
              "default: "s + (!opt.taxonStyle.showRankName ? "on" : "off"))
        ,
        option("-separate-cols", "-separatecols").set(opt.useSeparateCols)
            %("Prints *all* mapping information (rank, taxon name, taxon ids) "
              "in separate columns (see option '-separator').\n"
              "default: "s + (opt.useSeparateCols ? "on" : "off"))
        ,
        (   option("-separator") &
            value("text", [&](const string& arg) {
                    opt.tokens.column = sanitize_special_chars(arg);
                })
                .if_missing([&]{ err += "Text missing after '-separator'!"; })
        )
            % "Sets string that separates output columns.\n"
              "default: '\\t|\\t'"
        ,
        (   option("-comment") &
            value("text", opt.tokens.comment)
                .if_missing([&]{ err += "Text missing after '-comment'!"; })
        )
            %("Sets string that precedes comment (non-mapping) lines.\n"
              "default: '"s + opt.tokens.comment + "'")
        ,
        option("-queryids", "-query-ids").set(opt.showQueryIds)
            %("Show a unique id for each query.\n"
              "Note that in paired-end mode a query is a pair of two "
              "read sequences. This option will always be activated if "
              "option '-hits-per-seq' is given.\n"
              "default: "s + (opt.showQueryIds ? "on" : "off"))
    );
}



//-------------------------------------------------------------------
/// @brief build mode command-line options
clipp::group
classification_analysis_cli(classification_analysis_options& opt, error_messages& err)
{
    using namespace clipp;
    return (
        "ANALYSIS: ABUNDANCES" %
        (
            (option("-abundances", "-abundance").set(opt.showTaxAbundances) &
             opt_value("file", opt.abundanceFile))
                %("Show absolute and relative abundance of each taxon.\n"
                  "If a valid filename is given, the list will be written to "
                  "this file.\n"
                  "default: "s + (opt.showTaxAbundances ? "on" : "off"))
            ,
            (   option("-abundance-per") &
                value("rank", [&](const string& name) {
                        auto r = rank_from_name(name, err);
                        if(r < taxon_rank::root) opt.showAbundanceEstimatesOnRank = r;
                    })
                    .if_missing([&]{ err += "Taxonomic rank missing after '-abundance-per'!"; })
            )
                %("Show absolute and relative abundances for each "
                  "taxon on one specific rank.\n"
                  "Classifications on higher ranks will be estimated by "
                  "distributing them down according to the relative "
                  "abundances of classifications on or below the given rank. "
                  "(Valid values: "s + taxon_rank_names() + ")\n"s +
                  "If '-abundances <file>' was given, this list will "
                  "be printed to the same file.\n"
                  "default: "s + (false ? "on" : "off"))
        )
        ,
        "ANALYSIS: RAW DATABASE HITS" %
        (
            option("-allhits", "-all-hits").set(opt.showAllHits)
                %("For each query, print all feature hits in database.\n"
                  "default: "s + (opt.showAllHits ? "on" : "off"))
            ,
            option("-locations").set(opt.showLocations).set(opt.showTopHits)
                %("Show locations in candidate reference sequences.\n"
                  "Activates option '-tophits'.\n"
                  "default: "s + (opt.showLocations ? "on" : "off"))
            // ,
            // (   option("-hits-per-ref", "-hits-per-seq",
            //            "-hits-per-tgt", "-hits-per-target")
            //         .set(opt.showHitsPerTargetList) &
            //     opt_value("file", opt.targetMappingsFile)
            // )
            //     %("Unavailable in RNACache.")
        )
        ,
        "ANALYSIS: ALIGNMENTS" %
        group(
            option("-align", "-alignment").set(opt.showAlignment)
                %("Show semi-global alignment to best candidate reference sequence.\n"
                  "Original files of reference sequences must be available.\n"
                  "This feature decreases the querying speed!\n"
                  "default: "s + (opt.showAlignment ? "on" : "off"))
        )
    );
}



//-------------------------------------------------------------------
/// @brief build mode command-line options
clipp::group
classification_evaluation_cli(classification_evaluation_options& opt,
                              error_messages&)
{
    using namespace clipp;
    return (
        option("-ground-truth", "-groundtruth")
            .set(opt.determineGroundTruth).set(opt.showGroundTruth)
            %("Report correct query taxa if known.\n"
              "Queries need to have either a 'taxid|<number>' entry in "
              "their header or a sequence id that is also present in "
              "the database.\n"
              "This feature decreases querying speed!\n"
              "default: "s + (opt.determineGroundTruth ? "on" : "off"))
        ,
        option("-accuracy").set(opt.accuracy).set(opt.determineGroundTruth)
            %("Report accuracy stats "
              "by comparing query origins (ground truth) and mappings.\n"
              "Queries need to have either a 'taxid|<number>' entry in "
              "their header or a sequence id that is also found in "
              "the database.\n"
              "This feature might decrease querying speed!\n"
              "default: "s + (opt.accuracy ? "on" : "off"))
    );
}



//-------------------------------------------------------------------
/// @brief build mode command-line options
clipp::group
performance_options_cli(performance_tuning_options& opt, error_messages& err)
{
    using namespace clipp;
    return (
    (   option("-threads") &
        integer("#", opt.numThreads)
            .if_missing([&]{ err += "Number missing after '-threads'!"; })
    )
        %("Sets the maximum number of parallel threads to use."
          "default (on this machine): "s + to_string(opt.numThreads))
    ,
    #ifdef RC_BAM
    (   option("-bam-threads") &
        integer("#", opt.bamThreads)
            .if_missing([&]{ err += "Number missing after '-bam-threads'!"; })
    ),
    #endif
    (   option("-batch-size", "-batchsize") &
        integer("#", opt.batchSize)
            .if_missing([&]{ err += "Number missing after '-batch-size'!"; })
    )
        %("Process <#> many queries (reads or read pairs) per thread at once.\n"
          "default (on this machine): "s + to_string(opt.batchSize))
    ,
    #ifdef RC_BAM
    (   option("-bam-buffer") &
        integer("#", opt.bamBufSize)
            .if_missing([&]{ err += "Number missing after '-bam-buffer'!"; })
    )
    ,
    #endif
    (   option("-query-limit", "-querylimit") &
        integer("#", opt.queryLimit)
            .if_missing([&]{ err += "Number missing after '-query-limit'!"; })
    )
        %("Classify at max. <#> queries (reads or read pairs) per input file. "
          "and \n"
          "default: "s + (opt.queryLimit < 1 ? "none"s : to_string(opt.queryLimit))
    )
    );
}



//-------------------------------------------------------------------
/// @brief build mode command-line options
clipp::group
query_mode_cli(query_options& opt, error_messages& err)
{
    using namespace clipp;

    return (
    "BASIC PARAMETERS" %
    (
        database_parameter(opt.dbfile, err)
        ,
        opt_values(match::prefix_not{"-"}, "sequence file/directory", opt.infiles)
            % "FASTA or FASTQ files containing transcriptomic sequences "
              "(short reads, long reads, ...) "
              "that shall be classified.\n"
              "* If directory names are given, they will be searched for "
              "sequence files (at most 10 levels deep).\n"
              "* If no input filenames or directories are given, RNACache will "
              "run in interactive query mode. This can be used to load the database into "
              "memory only once and then query it multiple times with different "
              "query options. "
    ),
    "MAPPING RESULTS OUTPUT" %
    one_of(
        (   option("-out") &
            value("file", opt.queryMappingsFile)
                .if_missing([&]{ err += "Output filename missing after '-out'!"; })
        )
            % "Redirect output to file <file>.\n"
              "If not specified, output will be written to stdout. "
              "If more than one input file was given all output "
              "will be concatenated into one file."
        ,
        (   option("-split-out", "-splitout").set(opt.splitOutputPerInput) &
            value("file", opt.queryMappingsFile)
                .if_missing([&]{ err += "Output filename missing after '-split-out'!"; })
        )
            % "Generate output and statistics for each input file "
              "separately. For each input file <in> an output file "
              "with name <file>_<in> will be written."
    ),
    "PAIRED-END READ HANDLING" %
    (   one_of(
            option("-pairfiles", "-pair-files", "-paired-files")
            .set(opt.pairing, pairing_mode::files)
            % "Interleave paired-end reads from two consecutive files, "
              "so that the nth read from file m and the nth read "
              "from file m+1 will be treated as a pair. "
              "If more than two files are provided, their names "
              "will be sorted before processing. Thus, the order "
              "defined by the filenames determines the pairing not "
              "the order in which they were given in the command line."
            ,
            option("-pairseq", "-pair-seq", "-paired-seq")
            .set(opt.pairing, pairing_mode::sequences)
            % "Two consecutive sequences (1+2, 3+4, ...) from each file "
              "will be treated as paired-end reads."
        ),

        (   option("-insertsize", "-insert-size") &
            integer("#", opt.classify.insertSizeMax)
                .if_missing([&]{ err += "Number missing after '-insertsize'!"; })
        )
            % "Maximum insert size to consider.\n"
              "default: sum of lengths of the individual reads"
    )
    ,
    "CLASSIFICATION" %
        classification_params_cli(opt.classify, err)
    ,
    "GENERAL OUTPUT FORMATTING" % (
        option("-no-summary", "-nosummary").set(opt.output.showSummary,false)
            %("Dont't show result summary & mapping statistics at the "
              "end of the mapping output\n"
              "default: "s + (!opt.output.showSummary ? "on" : "off"))
        ,
        option("-no-query-params", "-no-queryparams", "-noqueryparams")
            .set(opt.output.showQueryParams,false)
            %("Don't show query settings at the beginning of the "
              "mapping output\n"
              "default: "s + (!opt.output.showQueryParams ? "on" : "off"))
        ,
        option("-no-err", "-noerr", "-no-errors").set(opt.output.showErrors,false)
            %("Suppress all error messages.\n"
              "default: "s + (!opt.output.showErrors ? "on" : "off"))
    )
    ,
    "CLASSIFICATION RESULT FORMATTING" %
        classification_output_format_cli(opt.output.format, err)
    ,
    classification_analysis_cli(opt.output.analysis, err)
    ,
    "ADVANCED: GROUND TRUTH BASED EVALUATION" %
        classification_evaluation_cli(opt.output.evaluate, err)
    ,
    "ADVANCED: CUSTOM QUERY SKETCHING (SUBSAMPLING)" %
        sketching_options_cli(opt.sketching, err)
    ,
    "ADVANCED: DATABASE MODIFICATION" %
        database_storage_options_cli(opt.dbconfig, err)
    ,
    "ADVANCED: PERFORMANCE TUNING / TESTING" %
        performance_options_cli(opt.performance, err)
    ,
    option("-sam", "-SAM", "-edlib").set(opt.rnaMode, rna_mode::sam)
    ,
    option("-bam").set(opt.rnaMode, rna_mode::bam)
    ,
    option("-search").set(opt.rnaMode, rna_mode::search),
    catch_unknown(err)
    );
}



//-------------------------------------------------------------------
query_options
get_query_options(const cmdline_args& args, query_options opt)
{
    error_messages err;

    auto cli = query_mode_cli(opt, err);

    auto result = clipp::parse(args, cli);

    if(!result || err.any()) {
        raise_default_error(err, "query", query_mode_usage());
    }

    replace_directories_with_contained_files(opt.infiles);

    if(opt.pairing == pairing_mode::files) {
        if(opt.infiles.size() > 1) {
            std::sort(opt.infiles.begin(), opt.infiles.end());
        } else {
            // TODO warning that pairing_mode::files requires at least 2 files
            opt.pairing = pairing_mode::none;
        }
    }

    // interprets numbers > 1 as percentage
    auto& cl = opt.classify;
    if(cl.covMin > 1) cl.covMin *= 0.01;
    if(cl.hitsCutoff > 1) cl.hitsCutoff *= 0.01;

    if(cl.maxNumCandidatesPerQuery < 1) {
        cl.maxNumCandidatesPerQuery = std::numeric_limits<size_t>::max();
    }


    // processing option checks
    auto& perf = opt.performance;
    if(perf.numThreads < 1) perf.numThreads = 1;
    if(perf.batchSize  < 1) perf.batchSize  = 1;
    if(perf.queryLimit < 0) perf.queryLimit = 0;

    #ifdef RC_BAM
    if (opt.rnaMode==rna_mode::bam)
        perf.bamBufSize = 1 << perf.bamBufSize;
    else
        perf.bamBufSize = 1;
    #endif


    //output file consistency checks
    auto& ana = opt.output.analysis;
    if(ana.targetMappingsFile == opt.queryMappingsFile) ana.targetMappingsFile.clear();
    if(ana.abundanceFile == opt.queryMappingsFile) ana.abundanceFile.clear();

    // output option checks and consistency

    //always show query ids if hits per target list requested
    auto& fmt = opt.output.format;

    if(ana.showHitsPerTargetList) fmt.showQueryIds = true;

    // modify output tokens for separate column printig
    if(fmt.useSeparateCols) {
        fmt.collapseUnclassifiedLineages = false;
        fmt.tokens.taxSeparator = fmt.tokens.column;
        fmt.tokens.rankSuffix   = fmt.tokens.column;
        fmt.tokens.taxidPrefix  = fmt.tokens.column;
        fmt.tokens.taxidSuffix  = "";
    }

    if(ana.showAllHits) {
        fmt.mapViewMode = map_view_mode::all;
    }

    // SAM output must be free of other content
    if (opt.rnaMode==rna_mode::sam && opt.queryMappingsFile.empty()) {
        opt.output.showQueryParams = false;
        opt.output.format.mapViewMode = map_view_mode::none;
        opt.output.showSummary = false;
    }

    return opt;
}



//-------------------------------------------------------------------
string query_mode_usage() {
    return
    "    rnacache query <database>\n\n"
    "    rnacache query <database> <sequence file/directory>... [OPTION]...\n\n"
    "    rnacache query <database> [OPTION]... <sequence file/directory>...";
}



//-------------------------------------------------------------------
string query_mode_examples() {
    return
    "    Query all sequences in 'myreads.fna' against pre-built database 'refseq':\n"
    "        rnacache query refseq myreads.fna -out results.txt\n"
    "\n"
    "    Query all sequences in multiple files against database 'refseq':\n"
    "        rnacache query refseq reads1.fna reads2.fna reads3.fna\n"
    "\n"
    "    Query all sequence files in folder 'test' againgst database 'refseq':\n"
    "        rnacache query refseq test\n"
    "\n"
    "    Query multiple files and folder contents against database 'refseq':\n"
    "        rnacache query refseq file1.fna folder1 file2.fna file3.fna folder2\n"
    "\n"
    "    Load database in interactive query mode, then query multiple read batches\n"
    "        rnacache query refseq\n"
    "        reads1.fa reads2.fa -pairfiles -insertsize 400\n"
    "        reads3.fa -pairseq -insertsize 300\n";
}



//-------------------------------------------------------------------
string query_mode_docs() {

    query_options opt;
    error_messages err;
    auto cli = query_mode_cli(opt, err);

    string docs = "SYNOPSIS\n\n";

    docs += query_mode_usage();

    docs += "\n\n\n"
        "DESCRIPTION\n"
        "\n"
        "    Map sequences (short reads, long reads, ...)\n"
        "    to their most likely reference sequences of origin.\n"
        "\n\n";

    docs += clipp::documentation(cli, cli_doc_formatting()).str();

    docs += "\n\n\nEXAMPLES\n\n";
    docs += query_mode_examples();

    docs += "\n\n"
        "OUTPUT FORMAT\n"
        "\n"
        "    RNACache's default read mapping output format is:\n"
        "    read_header | seq_name\n"
        "\n"
        "    Note that the separator '\\t|\\t' can be changed to something else with\n"
        "    the command line option '-separator <text>'.\n";
        // TODO: add more information about cli-parameters.

    return docs;
}





/*************************************************************************//**
 *
 *
 *  I N F O   M O D E
 *
 *
 *****************************************************************************/
/// @brief shared command-line options for taxonomy
clipp::group
info_mode_cli(info_options& opt, error_messages& err)
{
    using namespace clipp;

    return "PARAMETERS" % (
        database_parameter(opt.dbfile, err)
            .required(false).set(opt.mode, info_mode::db_config)
        ,
        one_of(
            command(""), // dummy
            (
                command("reference", "references", "ref",
                        "target", "targets", "tgt",
                        "sequence", "sequences", "seq")
                    .set(opt.mode, info_mode::targets),
                opt_values("sequence_id", opt.targetIds)
            ),
            (
                command("rank").set(opt.mode, info_mode::tax_ranks),
                value("rank_name", [&](const string& name) {
                        opt.rank = rank_from_name(name, err);
                    })
                    .if_missing([&]{ err += "Rank name missing!"; })
                    .doc("Valid values: "s + taxon_rank_names())
            ),
            command("lineages", "lineage", "lin")
                .set(opt.mode, info_mode::tax_lineages)
            ,
            command("statistics", "stat")
                .set(opt.mode, info_mode::db_statistics)
            ,
            command("locations", "loc", "featuremap", "features")
                .set(opt.mode, info_mode::db_feature_map)
            ,
            command("featurecounts")
                .set(opt.mode, info_mode::db_feature_counts)
        )
        ,
        catch_unknown(err)
    );

}



//-------------------------------------------------------------------
info_options
get_info_options(const cmdline_args& args)
{
    info_options opt;
    error_messages err;

    auto cli = info_mode_cli(opt, err);

    auto result = clipp::parse(args, cli);

    if(!result || err.any()) {
        raise_default_error(err, "info", info_mode_usage());
    }

    return opt;
}



//-------------------------------------------------------------------
string info_mode_usage()
{
    info_options opt;
    error_messages err;
    const auto cli = info_mode_cli(opt, err);
    return clipp::usage_lines(cli, "rnacache info", cli_usage_formatting()).str();
}



//-------------------------------------------------------------------
string info_mode_examples() {
    return
    "    List metadata for all reference sequences in database 'refseq':\n"
    "        rnacache info refseq.db ref\n"
    "\n"
    "    List metadata for the sequence with id NC_12345.6 in database 'refseq':\n"
    "        rnacache info refseq.db ref NC_12345.6\n";
}



//-------------------------------------------------------------------
string info_mode_docs() {

    info_options opt;
    error_messages err;
    const auto cli = info_mode_cli(opt, err);

    string docs = "SYNOPSIS\n\n";

    docs += clipp::usage_lines(cli, "rnacache info", cli_usage_formatting()).str();

    docs += "\n\n\n"
        "DESCRIPTION\n"
        "\n"
        "    Display (meta-)information stored in a database.\n"
        "\n\n"
        "SUB-MODES\n"
        "\n"
        "    rnacache info\n"
        "        show basic properties of RNACache executable (data type widths, etc.)\n"
        "\n"
        "    rnacache info <database>\n"
        "        show basic properties of <database>\n"
        "\n"
        "    rnacache info <database> ref[erence]\n"
        "       list meta information for all reference sequences in <database>\n"
        "\n"
        "    rnacache info <database> ref[erence] <sequence_id>...\n"
        "       list meta information for specific reference sequences\n"
        "\n"
        "    rnacache info <database> stat[istics]\n"
        "       print database statistics / hash table properties\n"
        "\n"
        "    rnacache info <database> loc[ations]\n"
        "       print map (feature -> list of reference locations)\n"
        "\n"
        "    rnacache info <database> featurecounts\n"
        "       print map (feature -> number of reference locations)\n"

        "\n\n";

    docs += clipp::documentation{cli, cli_doc_formatting()}.str();

    docs += "\n\n\nEXAMPLES\n\n";
    docs += info_mode_examples();

    return docs;
}



} // namespace mc
