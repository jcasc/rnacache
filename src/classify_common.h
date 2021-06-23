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
#include "sequence_view.h"
#include "timer.h"
#include "querying.h"

namespace mc {



/*************************************************************************//**
 *
 * @brief candidates + groundTruth
 *
 *****************************************************************************/
struct classification
{
    classification(classification_candidates cand):
        candidates{std::move(cand)}, groundTruth{database::nulltgt}
    {}

    classification_candidates candidates;
    target_id groundTruth;
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
                           std::ostream& alignmentOutputTgt = std::cout,
                           std::ostream& statusTarget    = std::cerr)
    :
        perReadOut(readOutputTgt),
        alignmentOut(alignmentOutputTgt),
        status(statusTarget)
    {}

    void flush_all_streams() {
        perReadOut.flush();
        alignmentOut.flush();
        status.flush();
    }

    std::ostream& perReadOut;
    std::ostream& alignmentOut; // BAM / SAM
    std::ostream& status;
    timer time;

    rna_mapping_statistics statistics;
};



/*************************************************************************//**
 *
 * @brief returns the ground truth target based on a sequence header info
 *        (needed for accuracy tests)
 *
 *****************************************************************************/
target_id
ground_truth_target(const database&, const std::string& header);



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

} // namespace mc


#endif
