/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller   (muellan@uni-mainz.de)
 *                       & Robin Kobus    (rkobus@uni-mainz.de)
 *                       & Julian Cascitt (jcascitt@students.uni-mainz.de)
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

/******************************************************************************
 *
 * @file types and functions used by both DNA and RNA mapping functions
 *
 *****************************************************************************/
#ifndef MC_CLASSIFY_COMMON_H_
#define MC_CLASSIFY_COMMON_H_


#include "config.h"
#include "candidates.h"
#include "classification_statistics.h"
#include "alignment.h"
#include "sequence_view.h"
#include "timer.h"

// SEE NOTES
#include "querying.h"
// !!!!!


namespace mc {


/// @brief forward declarations
// struct classification_options;
// struct classification_output_options;
// struct evaluation_options;
// struct query_options;


/// @brief typedefs
using taxon_list         = std::vector<const taxon*>;
using query_matches      = std::vector<taxon_list>;




/*************************************************************************//**
 *
 * @brief classification candidates + derived best classification
 *
 *****************************************************************************/
struct classification
{
    classification(classification_candidates cand):
        candidates{std::move(cand)}, best{nullptr}, groundTruth{nullptr}
    {}

    classification_candidates candidates;
    const taxon* best;
    const taxon* groundTruth;
};



/*************************************************************************//**
 *
 * @brief classification result target
 *
 *****************************************************************************/
struct classification_results
{
    explicit
    classification_results(std::ostream& readOutputTgt   = std::cout,
                           std::ostream& targetOutputTgt = std::cout,
                           std::ostream& taxonOutputTgt  = std::cout,
                           std::ostream& alignmentOutputTgt = std::cout,
                           std::ostream& statusTarget    = std::cerr)
    :
        perReadOut(readOutputTgt),
        perTargetOut(targetOutputTgt),
        perTaxonOut(taxonOutputTgt),
        alignmentOut(alignmentOutputTgt),
        status(statusTarget)
    {}

    void flush_all_streams() {
        perReadOut.flush();
        perTargetOut.flush();
        perTaxonOut.flush();
        alignmentOut.flush();
        status.flush();
    }

    std::ostream& perReadOut;
    std::ostream& perTargetOut; // unused in RNA
    std::ostream& perTaxonOut;  // unused in RNA
    std::ostream& alignmentOut; // RNA SAM mode
    std::ostream& status;
    timer time;

    rna_mapping_statistics statistics;
};



/*************************************************************************//**
 *
 * @brief Compare taxa by rank in descending order; root > ... > species.
 *        If ranks are equal, compare using sequence ids.
 *
 *****************************************************************************/
struct rank_higher {
    bool operator() (const taxon* lhs, const taxon* rhs) const noexcept {
        if(lhs->rank() > rhs->rank()) return true;
        if(lhs->rank() < rhs->rank()) return false;
        return lhs->id() < rhs->id();
    }
};

using taxon_count_map = std::map<const taxon*, float, rank_higher>;
// using taxon_count_map = std::unordered_map<const taxon*, query_id>;



/*************************************************************************//**
 *
 * @brief returns the ground truth taxon based on a sequence header info
 *        (needed for precision tests)
 *
 *****************************************************************************/
const taxon*
ground_truth_taxon(const database&, const std::string& header);


/*************************************************************************//**
 *
 * @brief returns the next main-rank ancestor taxon of
 *        a sequence's ground truth taxon (extracted from sequence header info)
 *
 *****************************************************************************/
const taxon*
ground_truth_ranked_taxon(const database&, const std::string& header);



/*************************************************************************//**
 *
 * @brief print header line for mapping table
 *
 *****************************************************************************/
void show_query_mapping_header(std::ostream&,
                               const classification_output_options&);



/*************************************************************************//**
 *
 * @brief generate classification candidates
 *
 *****************************************************************************/
template<class Locations>
classification_candidates
make_classification_candidates(const database& db,
                               const classification_options& opt,
                               const sequence_query& query,
                               const Locations& allhits);



/*************************************************************************//**
 *
 * @brief shows one query mapping line
 *        [query id], query_header, classification [, [top|all]hits list]
 *
 *****************************************************************************/
template<class Locations>
void show_query_mapping(
    std::ostream& os,
    const database& db,
    const classification_output_options& opt,
    const sequence_query& query,
    const classification& cls,
    const Locations& allhits);


/*************************************************************************//**
 *
 * @brief makes non-owning view to (sub)sequence
 *
 *****************************************************************************/
template<class Sequence>
inline auto
make_view_from_window_range(const Sequence& s, const window_range& range,
                            int size, int stride)
{
    auto end = s.begin() + (stride * range.end) + size;
    if(end > s.end()) end = s.end();

    return make_view(s.begin() + (stride * range.beg), end);
}



/*************************************************************************//**
 *
 * @brief performs a semi-global alignment
 *
 *****************************************************************************/
template<class Subject>
alignment<default_alignment_scheme::score_type,typename sequence::value_type>
make_semi_global_alignment(const sequence_query& query,
                           const Subject& subject)
{
    std::size_t score  = 0;
    std::size_t scorer = 0;

    const auto scheme = default_alignment_scheme{};

    //compute alignment
    auto align = align_semi_global(query.seq1, subject, scheme);
    score = align.score;
    //reverse complement
    auto query1r = make_reverse_complement(query.seq1);
    auto alignr = align_semi_global(query1r, subject, scheme);
    scorer = alignr.score;

    //align paired read as well
    if(!query.seq2.empty()) {
        score += align_semi_global_score(query.seq2, subject, scheme);
        auto query2r = make_reverse_complement(query.seq2);
        scorer += align_semi_global_score(query2r, subject, scheme);
    }

    return (score > scorer) ? align : alignr;
}



/*************************************************************************//**
 *
 * @brief compute alignment of top hits and optionally show it
 *
 *****************************************************************************/
void show_alignment(std::ostream&,
                    const database&,
                    const classification_output_options&,
                    const sequence_query&,
                    const classification_candidates&);



} // namespace mc


#endif
