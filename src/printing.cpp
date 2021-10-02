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

#include <ostream>
#include <utility>

#include "database.h"
#include "candidates.h"
#include "classification.h"
#include "classification_statistics.h"
#include "matches_per_target.h"
#include "stat_confusion.h"
#include "options.h"
#include "classify_common.h"
#include "printing.h"


namespace mc {

//-----------------------------------------------------------------------------
void show_query_parameters(std::ostream& os, const database& db, const query_options& opt)
{
    const auto& fmt = opt.output.format;
    const auto& comment = fmt.tokens.comment;

    if(fmt.showMapping) {
        os << comment << "Reporting per-read mappings (non-mapping lines start with '"
           << comment << "').\n";
    }
    else {
        os << comment << "Per-Read mappings will not be shown.\n";
    }

    if(opt.pairing == pairing_mode::files) {
        os << comment << "File based paired-end mode:\n"
           << comment << "  Reads from two consecutive files will be interleaved.\n"
           << comment << "  Max insert size considered " << opt.classify.insertSizeMax << ".\n";
    }
    else if(opt.pairing == pairing_mode::sequences) {
        os << comment << "Per file paired-end mode:\n"
           << comment << "  Reads from two consecutive sequences in each file will be paired up.\n"
           << comment << "  Max insert size considered " << opt.classify.insertSizeMax << ".\n";
    }

    os << comment << "Using " << opt.performance.numThreads << " threads\n";

    auto s = db.query_sketcher();
    os << comment
       << "k: " << int(s.kmer_size()) << '\n';

    os << comment
       << "Window length / stride / size (sketch): "
       << s.window_size() << " / " << s.window_stride() << " / " << s.sketch_size() << '\n';
    
    if (std::is_same<classification_candidates, best_distinct_matches_in_contiguous_window_ranges>()) {
        os << comment
           << "At maximum "
           << opt.classify.maxNumCandidatesPerQuery
           << " classification candidates will be considered per query.\n";
    }

    os << comment
       << "Minimum hit threshold: "
       << opt.classify.hitsMin << " per candidate\n";
    
    os << comment
       << "Hit cutoff: "
       << int(opt.classify.hitsCutoff*100) << "% of the maximal candidate\n";

    os << comment
       << "Coverage cutoff: "
       << int(opt.classify.covMin*100) << '%';

    if (opt.classify.covNorm == coverage_norm::max)
        os << " of the maximal candidate target coverage per query\n";
    else if (opt.classify.covNorm == coverage_norm::none)
        os << " of windows\n";
    
    if (opt.classify.covFill == coverage_fill::fill) {
        os << "Coverage includes caps. (2nd coverage condition waived.)\n";
    }
}



//-------------------------------------------------------------------
void show_target_header(std::ostream& os,
                       const classification_output_formatting& opt,
                       const std::string& prefix)
{
    const auto& style = opt.targetStyle;
    const auto& fmt = opt.tokens;

    if(style.showName) {
        os << prefix << "target_name";
        if(style.showId) {
            os << fmt.tgtidPrefix << prefix << "target_id" << fmt.tgtidSuffix;
        }
    }
    else if(style.showId) {
        os << prefix << "target_id";
    }
}



//-------------------------------------------------------------------
void print_target(std::ostream& os,
                 const std::string& taxName,
                 target_id id,
                 target_print_style style,
                 const formatting_tokens& fmt)
{
    if(style.showName) {
        os << taxName;
        if(style.showId) {
            os << fmt.tgtidPrefix << id << fmt.tgtidSuffix;
        }
    }
    else if(style.showId) {
        os << id;
    }
}



//-------------------------------------------------------------------
void show_target(std::ostream& os,
                const database& db,
                const classification_output_formatting& opt,
                target_id tgt)
{
    if(tgt == database::nulltgt) 
        print_target(os, opt.tokens.none, database::nulltgt, opt.targetStyle, opt.tokens);
    else  {
        print_target(os, db.get_target(tgt).name(), tgt, opt.targetStyle, opt.tokens);
    }
}



//-------------------------------------------------------------------
void show_candidates(std::ostream& os,
                     const database& db,
                     const classification_candidates& cands)
{
    using size_t = classification_candidates::size_type;

    for(size_t i = 0; i < cands.size(); ++i) {
        if(i > 0) os << ',';
        os << db.get_target(cands[i].tgt).name() << ':' << cands[i].hits;
    }

}



//-------------------------------------------------------------------
template<class Locations>
void show_matches(std::ostream& os,
                  const database& db,
                  const Locations& matches)
{
    if(matches.empty()) return;

    auto cur = matches.begin();
    int count = 1;
    for(auto it = matches.begin()+1; it != matches.end(); ++it) {
        if(*cur == *it)
            ++count;
        else {
            const auto tgt = db.get_target(cur->tgt);
            os << tgt.name()
               << '/' << int(cur->win)
               << ':' << count << ',';
            cur = it;
            count = 1;
        }
    }
    const auto tgt = db.get_target(cur->tgt);
    os << tgt.name()
       << '/' << int(cur->win)
       << ':' << count << ',';
}

template void show_matches<match_locations>(
    std::ostream& os,
    const database& db,
    const match_locations& matches);



//-------------------------------------------------------------------
void show_candidate_ranges(std::ostream& os,
                           const database& db,
                           const classification_candidates& cand)
{
    const auto w = db.target_sketcher().window_stride();

    for(const auto& c : cand) {
        os << '[' << (w * c.pos.beg)
           << ',' << (w * c.pos.end + db.target_sketcher().window_size()) << "] ";
    }
}



//-------------------------------------------------------------------
void show_statistics(std::ostream& os,
                     const rna_mapping_statistics& stats,
                     const std::string& prefix)
{
    os << prefix << "Mapping statistics:\n"
       << prefix << "Total Reads (or pairs): " << stats.total() << '\n'
       << prefix << "Total Matches:          " << stats.matches() << '\n'
       << prefix << "Reads Aligned:          " << stats.aligned() << '\n'
       << prefix << "Hits Per Read:          " << stats.hits_per_read() << '\n';
}

//-------------------------------------------------------------------
void show_accuracy(std::ostream& os,
                   const rna_mapping_statistics& stats,
                   const std::string& prefix)
{
    os << prefix << "Accuracy statistics:\n"
       << prefix << "Total reads (or pairs): " << stats.total() << '\n'
       << prefix << "Total matches:          " << stats.matches() << '\n'
       << prefix << "Origin found:           " << stats.correct() << '\n'
       << prefix << "Correctly rejected:     " << stats.true_negatives() << '\n'
       << prefix << "Reads aligned:          " << stats.aligned() << '\n'
       << prefix << "Recall:                 " << stats.recall() << '\n'
       << prefix << "Hits per read:          " << stats.hits_per_read() << '\n'
       << prefix << "Hits per aligned read:  " << stats.hits_per_aligned_read() << '\n'
       << prefix << "Total True Hit Rate:    " << stats.total_true_hit_rate() << '\n'
       << prefix << "Mean True Hit Rate:     " << stats.mean_true_hit_rate() << '\n';
}



/*************************************************************************//**
 *
 * @brief show summary and statistics of classification
 *
 *****************************************************************************/
void show_summary(const query_options& opt,
                  const classification_results& results)
{
    const auto& statistics = results.statistics;
    const auto numQueries = (opt.pairing == pairing_mode::none)
                            ? statistics.total() : 2 * statistics.total();

    const auto speed = numQueries / results.time.minutes();
    const auto& comment = opt.output.format.tokens.comment;
    results.mainOut
        << comment << "queries: " << numQueries << '\n'
        << comment << "time:    " << results.time.milliseconds() << " ms\n"
        << comment << "speed:   " << speed << " queries/min\n";

    if(statistics.total() > 0) {
        if (opt.output.evaluate.statistics) {
            if (opt.output.evaluate.determineGroundTruth)
                show_accuracy(results.mainOut, statistics, comment);
            else
                show_statistics(results.mainOut, statistics, comment);
        }
    } else
        results.mainOut << comment << "No valid query sequences found.\n";
}



//-------------------------------------------------------------------
void print_static_properties(const database& db)
{
    using target_id = database::target_id;
    using window_id = database::window_id;
    using feature_t = database::feature;
    using bkt_sz_t  = database::bucket_size_type;

    std::cout
        << "------------------------------------------------\n"
        << "RNACache version     " << RC_VERSION_STRING << " (" << RC_VERSION << ")\n"
        << "database version     " << RC_DB_VERSION << '\n'
        << "------------------------------------------------\n"
        << "sequence type        " << type_name<database::sequence>() << '\n'
        << "target id type       " << type_name<target_id>() << " " << (sizeof(target_id)*CHAR_BIT) << " bits\n"
        << "target limit         " << std::uint64_t(db.max_target_count()) << '\n'
        << "------------------------------------------------\n"
        << "window id type       " << type_name<window_id>() << " " << (sizeof(window_id)*CHAR_BIT) << " bits\n"
        << "window limit         " << std::uint64_t(db.max_windows_per_target()) << '\n'
        << "window length        " << db.target_sketcher().window_size() << '\n'
        << "window stride        " << db.target_sketcher().window_stride() << '\n'
        << "------------------------------------------------\n"
        << "sketcher type        " << type_name<database::sketcher>() << '\n'
        << "feature type         " << type_name<feature_t>() << " " << (sizeof(feature_t)*CHAR_BIT) << " bits\n"
        << "feature hash         " << type_name<database::feature_hash>() << '\n'
        << "kmer size            " << std::uint64_t(db.target_sketcher().kmer_size()) << '\n'
        << "kmer limit           " << std::uint64_t(db.target_sketcher().max_kmer_size()) << '\n'
        << "sketch size          " << db.target_sketcher().sketch_size() << '\n'
        << "------------------------------------------------\n"
        << "bucket size type     " << type_name<bkt_sz_t>() << " " << (sizeof(bkt_sz_t)*CHAR_BIT) << " bits\n"
        << "max. locations       " << std::uint64_t(db.max_locations_per_feature()) << '\n'
        << "location limit       " << std::uint64_t(db.max_supported_locations_per_feature()) << '\n'
        << "------------------------------------------------"
        << std::endl;
}



//-----------------------------------------------------------------------------
void print_content_properties(const database& db)
{
    if(db.target_count() > 0) {

        std::cout
        << "targets              " << db.target_count() << '\n';
    }

    if(db.feature_count() > 0) {
        auto lss = db.location_list_size_statistics();

        std::cout
        << "------------------------------------------------\n"
        << "buckets              " << db.bucket_count() << '\n'
        << "bucket size          " << "max: " << lss.max()
                                   << " mean: " << lss.mean()
                                   << " +/- " << lss.stddev()
                                   << " <> " << lss.skewness() << '\n'
        << "features             " << db.feature_count() << '\n'
        << "dead features        " << db.dead_feature_count() << '\n'
        << "locations            " << db.location_count() << '\n';
    }
    std::cout
        << "------------------------------------------------\n";
}



//-----------------------------------------------------------------------------
void print_coverages(const database& db, matches_per_target_light& mpt)
{
    for (size_t i = 0; i < db.target_count(); ++i) {
        auto& tgt = db.get_target(i);
        std::cout << tgt.name() << " " << tgt.source().windows << " " << mpt.num_hits(i) << " " << double(mpt.num_hits(i))/tgt.source().windows << "\n";
    }
}



//-----------------------------------------------------------------------------
void print_coverage_profiles(const database& db, matches_per_target_light& mpt)
{
    for (size_t i = 0; i < db.target_count(); ++i) {
        auto& tgt = db.get_target(i);
        std::cout << tgt.name() << ' ' << tgt.source().windows << ' ' << mpt.num_hits(i) << " " << double(mpt.num_hits(i))/tgt.source().windows;
        if (mpt.num_hits(i) > 0) {
            std::cout << " :";
            for (auto& win: mpt.find(i)->second)
                std::cout << ' ' << win;
        }
        std::cout << '\n';
    }
} 




} // namespace mc
