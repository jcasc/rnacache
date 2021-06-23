/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller    (muellan@uni-mainz.de)
 *                       & Robin Kobus     (rkobus@uni-mainz.de)
 *                       & Julian Cascitti (jcascitt@students.uni-mainz.de)
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

#include <thread>

#include "classify_common.h"
#include "sequence_io.h"
#include "candidates.h"
#include "printing.h"



namespace mc {

using std::string;


/*************************************************************************//**
 *
 * @brief print header line for mapping table
 *
 *****************************************************************************/
void show_query_mapping_header(std::ostream& os,
                               const classification_output_options& opt)
{
    if(opt.format.mapViewMode == map_view_mode::none) return;

    const auto& colsep = opt.format.tokens.column;

    os << opt.format.tokens.comment << "TABLE_LAYOUT: ";

    if(opt.format.showQueryIds) os << "query_id" << colsep;

    os << "query_header" << colsep;

    if(opt.evaluate.showGroundTruth) {
        show_target_header(os, opt.format, "truth_");
        os << colsep;
    }

    if(opt.analysis.showAllHits) os << "all_hits" << colsep;
    if(opt.analysis.showTopHits) os << "top_hits" << colsep;
    if(opt.analysis.showLocations) os << "candidate_locations" << colsep;

    os << '\n';
}



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
                               const Locations& allhits)
{
    candidate_generation_rules rules;

    rules.maxWindowsInRange = window_id( 2 + (
        std::max(query.seq1.size() + query.seq2.size(), opt.insertSizeMax) /
        db.target_sketcher().window_stride() ));

    rules.maxCandidates = opt.maxNumCandidatesPerQuery;

    return classification_candidates{allhits, rules};

}

template classification_candidates
make_classification_candidates<match_locations>(const database& db,
                               const classification_options& opt,
                               const sequence_query& query,
                               const match_locations& allhits);


//-----------------------------------------------------------------------------
target_id
ground_truth_target(const database& db, const string& header)
{
    //try to extract query id and find the corresponding target in database
    target_id tgt = database::nulltgt;
    tgt = db.target_with_name(extract_accession_string(header, sequence_id_type::acc_ver));
    if(tgt) return tgt;

    tgt = db.target_with_similar_name(extract_accession_string(header, sequence_id_type::acc));
    if(tgt) return tgt;

    //try to extract id from header
    tgt = extract_target_id(header);
    if(tgt) return tgt;

    //try to find entire header as sequence identifier
    tgt = db.target_with_name(header);

    return tgt;
}


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
    const Locations& allhits)
{
    const auto& fmt = opt.format;

    if(fmt.mapViewMode == map_view_mode::none ||
        (fmt.mapViewMode == map_view_mode::mapped_only &&
        cls.candidates.empty()))
    {
        return;
    }

    const auto& colsep = fmt.tokens.column;

    if(fmt.showQueryIds) os << query.id << colsep;

    //print query header (first contiguous string only)
    auto l = query.header.find(' ');
    if(l != string::npos) {
        auto oit = std::ostream_iterator<char>{os, ""};
        std::copy(query.header.begin(), query.header.begin() + l, oit);
    }
    else {
        os << query.header;
    }
    os << colsep;

    if(opt.evaluate.showGroundTruth) {
        show_target(os, db, opt.format, cls.groundTruth);
        os << colsep;
    }

    if(opt.analysis.showAllHits) {
        show_matches(os, db, allhits);
        os << colsep;
    }
    
    if(opt.analysis.showTopHits)
    {
        show_candidates(os, db, cls.candidates);
        os << colsep;
    }
    if(opt.analysis.showLocations) {
        show_candidate_ranges(os, db, cls.candidates);
        os << colsep;
    }

    os << '\n';
}


template void show_query_mapping<match_locations>(
    std::ostream& os,
    const database& db,
    const classification_output_options& opt,
    const sequence_query& query,
    const classification& cls,
    const match_locations& allhits);


} // namespace mc
