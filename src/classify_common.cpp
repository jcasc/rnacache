/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller   (muellan@uni-mainz.de)
 *                       & Robin Kobus    (rkobus@uni-mainz.de)
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

#include "sequence_io.h"
#include "candidates.h"
#include "printing.h"

// SEE NOTES
#include "querying.h"
// !!!!!

#include "classify_common.h"


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
        show_taxon_header(os, opt.format, "truth_");
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

    rules.mergeBelow    = opt.lowestRank;
    rules.maxCandidates = opt.maxNumCandidatesPerQuery;

    return classification_candidates{db, allhits, rules};

}


// WARNING !!!!!!!!!!!!!!!!!
template classification_candidates
make_classification_candidates<match_locations>(const database& db,
                               const classification_options& opt,
                               const sequence_query& query,
                               const match_locations& allhits);


//-----------------------------------------------------------------------------
const taxon*
ground_truth_taxon(const database& db, const string& header)
{
    //try to extract query id and find the corresponding target in database
    const taxon* tax = nullptr;
    tax = db.taxon_with_name(extract_accession_string(header, sequence_id_type::acc_ver));
    if(tax) return tax;

    tax = db.taxon_with_similar_name(extract_accession_string(header, sequence_id_type::acc));
    if(tax) return tax;

    //try to extract id from header
    tax = db.taxon_with_id(extract_taxon_id(header));
    if(tax) return tax;

    //try to find entire header as sequence identifier
    tax = db.taxon_with_name(header);

    return tax;
}


//-----------------------------------------------------------------------------
const taxon*
ground_truth_ranked_taxon(const database& db, const string& header)
{
    const taxon* tax = ground_truth_taxon(db, header);
    if(tax) return db.next_ranked_ancestor(tax);

    return nullptr;
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
        show_taxon(os, db, opt.format, cls.groundTruth);
        os << colsep;
    }

    if(opt.analysis.showAllHits) {
        show_matches(os, db, allhits, fmt.lowestRank);
        os << colsep;
    }
    
    if(opt.analysis.showTopHits)
    {
        show_candidates(os, db, cls.candidates, fmt.lowestRank);
        os << colsep;
    }
    if(opt.analysis.showLocations) {
        show_candidate_ranges(os, db, cls.candidates);
        os << colsep;
    }

    os << '\n';
}


// TEMPLATE WARNING TEMPLATE
template void show_query_mapping<match_locations>(
    std::ostream& os,
    const database& db,
    const classification_output_options& opt,
    const sequence_query& query,
    const classification& cls,
    const match_locations& allhits);



/*************************************************************************//**
 *
 * @brief compute alignment of top hits and optionally show it
 *
 *****************************************************************************/
void show_alignment(std::ostream& os,
                    const database& db,
                    const classification_output_options& opt,
                    const sequence_query& query,
                    const classification_candidates& tophits)
{
    //try to align to top target
    const taxon* tgtTax = tophits[0].tax;
    if(tgtTax && tgtTax->rank() == taxon_rank::Sequence) {
        const auto& src = tgtTax->source();
        try {
            //load candidate file and forward to sequence
            auto reader = make_sequence_reader(src.filename);
            reader->skip(src.index-1);

            if(reader->has_next()) {
                auto tgtSequ = reader->next().data;
                auto subject = make_view_from_window_range(
                               tgtSequ, tophits[0].pos,
                               db.target_sketcher().window_size(),
                               db.target_sketcher().window_stride());

                auto align = make_semi_global_alignment(query, subject);

                //print alignment to top candidate
                const auto w = db.target_sketcher().window_stride();
                const auto& comment = opt.format.tokens.comment;
                os  << '\n'
                    << comment << "  score  " << align.score
                    << "  aligned to "
                    << src.filename << " #" << src.index
                    << " in range [" << (w * tophits[0].pos.beg)
                    << ',' << (w * tophits[0].pos.end + w) << "]\n"
                    << comment << "  query  " << align.query << '\n'
                    << comment << "  target " << align.subject;
            }
        }
        catch(std::exception& e) {
            if(opt.showErrors) std::cerr << e.what() << '\n';
        }
    }
}


} // namespace mc
