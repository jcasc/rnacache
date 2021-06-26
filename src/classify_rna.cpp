/******************************************************************************
 *
 * RNACache - Transcriptomic Mapping Tool
 *
 * Copyright (C) 2016-2021 Julian Cascitti (jcascitt@students.uni-mainz.de)
 *                       & André Müller    (muellan@uni-mainz.de)
 *                       & Robin Kobus     (rkobus@uni-mainz.de)
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
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <memory>
#include <iomanip>

#include <chrono>

#include "querying.h"
#include "candidates.h"
#include "matches_per_target.h"
#include "printing.h"
#include "querying.h"
#include "options.h"

#include "classification.h"
#include "classify_common.h"

#include "../dep/edlib.h"

#ifdef RC_BAM
#include <sam.h>
#endif


namespace mc {

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;

#ifdef RC_BAM
using bam_vector = std::vector<bam1_t>;
#endif

/*************************************************************************//**
 *
 * @brief print header line for mapping table
 *
 *****************************************************************************/
void show_query_mapping_header(std::ostream& os,
                               const classification_output_options& opt)
{
    if(!opt.format.showMapping) return;

    const auto& colsep = opt.format.tokens.column;

    os << opt.format.tokens.comment << "TABLE_LAYOUT: ";

    if(opt.format.showQueryIds) os << "query_id" << colsep;

    os << "query_header" << colsep;

    if(opt.evaluate.showGroundTruth) {
        show_target_header(os, opt.format, "truth_");
        os << colsep;
    }

    if(opt.analysis.showAllHits) os << "all_hits" << colsep;
    os << "top_hits" << colsep;
    if(opt.analysis.showLocations) os << "candidate_locations" << colsep;

    os << '\n';
}


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
 * @brief
 *
 *****************************************************************************/
template<class Consumer>
void do_in_parallel(int numThreads, std::size_t iterations,
                    Consumer&& consume)
{
    const auto prt = int(iterations / numThreads);
    const auto res = int(iterations % numThreads);

    auto job = [&](int tid) {
        const auto beg = tid * prt + std::min(tid, res);
        const auto end = (tid+1) * prt + std::min(tid+1, res);

        for(int i = beg; i < end; ++i) {
            consume(i);
        }
    };

    std::vector<std::thread> threads;
    for(int tid = 0; tid < numThreads; ++tid) {
        threads.emplace_back(job, tid);
    }
    for(auto& t: threads) {
        t.join();
    }
}





/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
template<class Action>
void for_each_eligible_mapping(const classification_options& opt,
                      const classification_candidates& cands,
                      Action&& action)
{
    match_candidate::count_type maxHits = 0;
    for(const auto& cand: cands) {
        maxHits = std::max(maxHits, cand.hits);
    }

    for(const auto& cand: cands) {
        if(cand.hits >= opt.hitsMin
            && double(cand.hits)/maxHits >= opt.hitsCutoff)
        {
            action(cand);
        }
    }
}

/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
template<class Action>
void for_each_eligible_mapping(const database& db,
                      const classification_options& opt,
                      const sequence_query& query,
                      const match_locations& allhits,
                      Action&& action)
{
    auto cands = make_classification_candidates(db, opt, query, allhits);

    for_each_eligible_mapping(opt, cands, std::forward<Action>(action));
}



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
void hits_cutoff_filter(const classification_options& opt,
                        classification_candidates& cands)
{
    match_candidate::count_type maxHits = 0;
    for (const auto& cand: cands)
        maxHits = std::max(maxHits, cand.hits);
        
    cands.erase(
        std::remove_if(cands.begin(), cands.end(),
            [&](match_candidate& cand) {
                // I don't remember why I checked for nullptr in for_each_eligible, but I'll do it to be safe.
                return !(cand.hits >= opt.hitsMin && double(cand.hits)/maxHits >= opt.hitsCutoff);
            }),
    cands.end());
}


/*************************************************************************//**
 *
 * @brief applies coverage filter to a list of taxa
 *
 *****************************************************************************/
void coverage_filter(const database& db,
                     const classification_options& opt,
                     classification_candidates& cands,
                     const matches_per_target_light& mpt)
{
    double norm = 0.0;

    switch(opt.covNorm) {
        default:
            norm = 1.0;
            break;
        case coverage_norm::max:
            for(const auto& cand: cands) {
                double cov = double(mpt.num_hits(cand.tgt))
                             / db.get_target(cand.tgt).source().windows;
                norm = std::max(norm, cov);
            }
            if(std::abs(norm) <= std::numeric_limits<double>::min()) return;
            norm = 1.0 / norm;
            break;
    }

    cands.erase(
        std::remove_if(cands.begin(), cands.end(),
        [&](match_candidate& cand) {
            double cov = double(mpt.num_hits(cand.tgt))
                         / db.get_target(cand.tgt).source().windows;
            return cov * norm < opt.covMin;
        }), cands.end());
}





/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
template<class Locations>
classification
classify(const database& db,
         const classification_options& opt,
         const sequence_query& query,
         const Locations& allhits,
         const matches_per_target_light& cov)
{
    classification cls { make_classification_candidates(db, opt, query, allhits) };

    hits_cutoff_filter(opt, cls.candidates);

    coverage_filter(db, opt, cls.candidates, cov);

    return cls;
}



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
void evaluate_classification(
    const classification_evaluation_options& opt,
    classification& cls,
    rna_mapping_statistics& statistics)
{
    if (opt.determineGroundTruth) {
        if (cls.groundTruth == database::nulltgt) {
            if (cls.candidates.empty()) {
                return statistics.mapped_noise();
            }
        } else {
            for (const auto& cand: cls.candidates) {
                if (cand.tgt == cls.groundTruth) {
                    return statistics.mapped_correct(cls.candidates.size());
                }
            }
        }
    }
    return statistics.mapped(cls.candidates.size());
}

#ifdef RC_BAM
struct bam_buffer
{
    bam_buffer() = default;
    bam_buffer(size_t n) : buf(n) {}

    void add_bam(size_t l_qname, const char *qname,
        uint16_t flag, int32_t tid, hts_pos_t pos, uint8_t mapq,
        size_t n_cigar, const uint32_t *cigar,
        int32_t mtid, hts_pos_t mpos, hts_pos_t isize,
        size_t l_seq, const char *seq, const char *qual,
        size_t l_aux)
    {
        vec.emplace_back();
        bam1_t& b = vec.back();
        bam_set_mempolicy(&b, BAM_USER_OWNS_STRUCT);
        if (!full) {
            bam_set_mempolicy(&b, BAM_USER_OWNS_STRUCT | BAM_USER_OWNS_DATA);
            b.data = &buf[used];
            b.m_data = buf.size()-used;
        }

        int ret = bam_set1(&b, l_qname, qname, flag, tid, pos, mapq, n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual, l_aux);
        
        if (ret < 0) {
            std::cerr << "!!! bam_set1 ERROR: " << ret << " !!!!" << std::endl;
        } else {
        }

        if ((bam_get_mempolicy(&b) & BAM_USER_OWNS_DATA) != 0) {
            b.m_data = std::min(b.m_data, (uint32_t(b.l_data) + 7) & (~7U));
            used += b.m_data;
        } else {
            full = true;
        }

    }

    bam_vector vec;

private:
    std::vector<uint8_t> buf;
    bool full = 0;
    size_t used = 0;
};
#endif


/*************************************************************************//**
 *
 * @brief per-batch buffer for output and (target -> hits) lists
 *
 *****************************************************************************/
struct mappings_buffer
{
    std::ostringstream out;
    std::ostringstream align_out;

    matches_per_target_light coverage;

    #ifdef RC_BAM
    bam_buffer bam_buf;
    mappings_buffer() = default;
    mappings_buffer(size_t n) : bam_buf(n) {}
    #endif
};


/*************************************************************************//**
 *
 * @brief edlib alignment container
 *
 *****************************************************************************/
struct edlib_alignment {

    enum struct status {FORWARD, REVERSE, UNALIGNED};

    edlib_alignment(const std::string& query, target_id tgt, const database& db, int max_edit_distance):
        tgt_(tgt), status_(status::UNALIGNED), score_(query.size()), cigar_(nullptr)
    {
        const std::string& target = db.get_target(tgt).seq();

        auto edlib_config = edlibNewAlignConfig(max_edit_distance, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities.data(), additionalEqualities.size());
        
        auto regular = edlibAlign(query.c_str(), query.size(), target.c_str(), target.size(), edlib_config);
        std::string reverse_query = make_reverse_complement(query);
        auto reverse_complement = edlibAlign(reverse_query.c_str(), reverse_query.size(), target.c_str(), target.size(), edlib_config);
        
        if (regular.status != EDLIB_STATUS_OK || regular.status != EDLIB_STATUS_OK) {
            throw std::runtime_error{"edlib failed!"};
        }

        // keep in mind ../dep/edlip.cpp:212
        // start loc is always on reference, but end can be negative
        if (regular.editDistance >= 0 && regular.endLocations[0] >= 0 &&
           (reverse_complement.editDistance < 0 || reverse_complement.endLocations[0] < 0 ||
            reverse_complement.editDistance >= regular.editDistance))
        {
            status_ = status::FORWARD;
            start_ = regular.startLocations[0];
            end_ = regular.endLocations[0];
            score_ = regular.editDistance;
            cigar_ = edlibAlignmentToCigar(regular.alignment, regular.alignmentLength, EDLIB_CIGAR_STANDARD);
        } 
        else if (reverse_complement.editDistance >= 0 && reverse_complement.endLocations[0] >= 0)
        {
            status_ = status::REVERSE;
            start_ = reverse_complement.startLocations[0];
            end_ = reverse_complement.endLocations[0];
            score_ = reverse_complement.editDistance;
            cigar_ = edlibAlignmentToCigar(reverse_complement.alignment, reverse_complement.alignmentLength, EDLIB_CIGAR_STANDARD);
        }
        edlibFreeAlignResult(regular);
        edlibFreeAlignResult(reverse_complement);
    }

    edlib_alignment(const edlib_alignment&) = delete;
    edlib_alignment& operator=(const edlib_alignment&) = delete;
    edlib_alignment& operator=(edlib_alignment&&) = delete;
    
    edlib_alignment(edlib_alignment&& other):
        tgt_(other.tgt_), status_(other.status_), score_(other.score_), start_(other.start_), end_(other.end_), cigar_(other.cigar_)
    {
        other.cigar_ = nullptr;
    }

    ~edlib_alignment() {free(cigar_);}

    bool aligned() const noexcept {return status_!=status::UNALIGNED;}
    int score() const noexcept {return score_;}
    status orientation() const noexcept {return status_;}
    target_id tgt() const noexcept {return tgt_;}
    int start() const noexcept {return start_;}
    int end() const noexcept {return end_;}
    const char * cigar() const noexcept {return cigar_;}

private:
    static std::vector<EdlibEqualityPair> additionalEqualities;
    target_id tgt_;
    status status_;
    int score_;
    int start_, end_;
    char * cigar_;
};

std::vector<EdlibEqualityPair> edlib_alignment::additionalEqualities({{'a', 'A'}, {'t', 'T'}, {'c', 'C'}, {'g', 'G'}});

struct edlib_alignment_pair {
    edlib_alignment_pair(const sequence_query& query, target_id tgt, const database& db, int max_edit_distance):
        first(query.seq1, tgt, db, max_edit_distance),
        second(query.seq2, tgt, db, max_edit_distance)
    {}

    bool aligned() const noexcept {return first.aligned() || second.aligned();}
    int score() const noexcept {return first.score() + second.score();}
    target_id tgt() const noexcept {return first.tgt();};

    edlib_alignment first, second;
};

using alns_vector = std::vector<edlib_alignment_pair>;



void show_sam_minimal(std::ostream& os, const target& tgt, const sequence_query& query, bool primary) {

    // function only applicable for mapped reads

    std::string qname = query.header.substr(0, query.header.size() - 2);
    size_t tgtlen = tgt.seq().size();
    size_t readlen = query.seq1.size();
    int64_t tlen = std::min(tgtlen, readlen);

    uint16_t flag1 = 0, flag2 = 0;
    
    // function only applicable for paired reads atm
    flag1 |= 0x1;
    flag2 |= 0x1;

    flag1 |= 0x2;
    flag2 |= 0x2;

    flag1 |= 0x40;
    flag2 |= 0x80;

    flag1 |= 0x20;
    flag2 |= 0x10;

    if (!primary) {
        flag1 |= 0x100;
        flag2 |= 0x100;
    }
    
    //----------------------------------- 
    // mate 1

    // QNAME
    os << qname << '\t';

    // FLAG
    os << flag1 << '\t';

    // RNAME
    os << tgt.header() << '\t';

    // POS
    os << 1 << '\t';

    // MAPQ
    os << 255 << '\t';

    // CIGAR
    if (tgtlen < readlen) {
        os << tgtlen << "M" << (readlen-tgtlen) << "I\t";
    } else {
        os << readlen << "M" << "\t";
    }

    // RNEXT
    os << "=" << '\t';

    // PNEXT
    os << 1 << '\t';

    // TLEN
    os << tlen << '\t';

    // SEQ
    os << query.seq1 << '\t';

    // QUAL
    os << "*";

    os << '\n';

    //----------------------------------- 
    // mate 2

    // QNAME
    os << qname << '\t';

    // FLAG
    os << flag2 << '\t';

    // RNAME
    os << tgt.header() << '\t';

    // POS
    os << 1  << '\t';

    // MAPQ
    os << 255 << '\t';

    // CIGAR
    if (tgtlen < readlen) {
        os << tgtlen << "M" << (readlen-tgtlen) << "I\t";
    } else {
        os << readlen << "M" << "\t";
    }

    // RNEXT
    os << "=" << '\t';

    // PNEXT
    os << 1 << '\t';

    // TLEN
    os << -tlen << '\t';

    // SEQ
    os << make_reverse_complement(query.seq2) << '\t';

    // QUAL
    os << "*";

    os << '\n';
}



void show_sam_alignment(std::ostream& os,
                        const database& db,
                        const sequence_query& query,
                        const edlib_alignment_pair& alignment,
                        bool primary)
{
    const target& tgt = db.get_target(alignment.tgt());

    uint16_t flag1 = 0, flag2 = 0;
    flag1 |= 0x1;
    flag2 |= 0x1;
    if (alignment.first.aligned() && alignment.second.aligned()) {
        flag1 |= 0x2;
        flag2 |= 0x2;
    }
    if (!alignment.first.aligned()) {
        flag1 |= 0x4;
        flag2 |= 0x8;
    }
    if (!alignment.second.aligned()) {
        flag1 |= 0x8;
        flag2 |= 0x4;
    }
    if (alignment.first.orientation() == edlib_alignment::status::REVERSE) {
        flag1 |= 0x10;
        flag2 |= 0x20;
    }
    if (alignment.second.orientation() == edlib_alignment::status::REVERSE) {
        flag1 |= 0x20;
        flag2 |= 0x10;
    }
    flag1 |= 0x40;
    flag2 |= 0x80;
    if (!primary) {
        flag1 |= 0x100;
        flag2 |= 0x100;
    }

    int64_t tlen = 0;
    if (alignment.first.aligned() && alignment.second.aligned()) 
    {
        tlen = std::max(alignment.first.end(), alignment.second.end())
             - std::min(alignment.first.start(), alignment.second.start()) + 1;
    }
    
    //----------------------------------- 
    // mate 1

    // QNAME
    os << query.header.substr(0, query.header.size() - 2) << '\t';

    // FLAG
    os << flag1 << '\t';

    // RNAME
    os << tgt.header()  << '\t';

    // POS
    os << (alignment.first.aligned() ? alignment.first.start() : alignment.second.start()) + 1  << '\t';

    // MAPQ
    os << 255 << '\t';

    // CIGAR
    os << (alignment.first.aligned() ? alignment.first.cigar() : "*")  << '\t';

    // RNEXT
    os << "=" << '\t';

    // PNEXT
    os << (alignment.second.aligned() ? alignment.second.start() : alignment.first.start())  + 1 << '\t';

    // TLEN
    os << (alignment.first.start() <= alignment.second.start() ? tlen: -tlen) << '\t';

    // SEQ
    if (alignment.first.orientation() == edlib_alignment::status::REVERSE)
        os << make_reverse_complement(query.seq1) << '\t';
    else
        os << query.seq1 << '\t';

    // QUAL
    os << "*" << '\t';

    os << '\n';

    //----------------------------------- 
    // mate 2

    // QNAME
    os << query.header.substr(0, query.header.size()-2) << '\t';

    // FLAG
    os << flag2 << '\t';

    // RNAME
    os << tgt.header()  << '\t';

    // POS
    os << (alignment.second.aligned() ? alignment.second.start() : alignment.first.start()) + 1  << '\t';

    // MAPQ
    os << 255 << '\t';

    // CIGAR
    os << (alignment.second.aligned() ? alignment.second.cigar() : "*")  << '\t';

    // RNEXT
    os << "=" << '\t';

    // PNEXT
    os << (alignment.first.aligned() ? alignment.first.start() : alignment.second.start())  + 1 << '\t';

    // TLEN
    os << (alignment.second.start() < alignment.first.start() ? tlen: -tlen) << '\t';

    // SEQ
    if (alignment.second.orientation() == edlib_alignment::status::REVERSE)
        os << make_reverse_complement(query.seq2) << '\t';
    else
        os << query.seq2 << '\t';

    // QUAL
    os << "*" << '\t';

    os << '\n';
}

#ifdef RC_BAM
void show_bam_alignment(bam_buffer& bam_buf, const sequence_query& query, const edlib_alignment_pair& alignment, bool primary) {

    // function only applicable for mapped reads atm
    // function only applicable for paired reads atm

    target_id tgt_id = alignment.tgt();

    uint16_t flag1 = BAM_FPAIRED | BAM_FREAD1;
    uint16_t flag2 = BAM_FPAIRED | BAM_FREAD2;

    if (!primary) {
        flag1 |= BAM_FSECONDARY;
        flag2 |= BAM_FSECONDARY;
    }
    if (!alignment.first.aligned()) {
        flag1 |= BAM_FUNMAP;
        flag2 |= BAM_FMUNMAP;
    }
    if (!alignment.second.aligned()) {
        flag1 |= BAM_FMUNMAP;
        flag2 |= BAM_FUNMAP;
    }
    if (alignment.first.aligned() && alignment.second.aligned()) {
        flag1 |= BAM_FPROPER_PAIR;
        flag2 |= BAM_FPROPER_PAIR;
    }
    if (alignment.first.orientation() == edlib_alignment::status::REVERSE) {
        flag1 |= BAM_FREVERSE;
        flag2 |= BAM_FMREVERSE;
    }
    if (alignment.second.orientation() == edlib_alignment::status::REVERSE) {
        flag1 |= BAM_FMREVERSE;
        flag2 |= BAM_FREVERSE;
    }

    int64_t tlen1 = 0;
    if (alignment.first.aligned() && alignment.second.aligned())
    {
        tlen1 = std::max(alignment.first.end(), alignment.second.end())
             - std::min(alignment.first.start(), alignment.second.start()) + 1;
    }
    int64_t tlen2 = tlen1;
    if (alignment.first.start() <= alignment.second.start())
        tlen2 = -tlen2;
    else
        tlen1 = -tlen1;
    
    size_t pos1 = (alignment.first.aligned() ? alignment.first.start() : alignment.second.start());
    size_t pos2 = (alignment.second.aligned() ? alignment.second.start() : alignment.first.start());

    size_t l_cigar = 0;
    uint32_t *cigar = nullptr;
    const char* seq = query.seq1.data();
    string rcseq;
    
    if (alignment.first.aligned())
        sam_parse_cigar(alignment.first.cigar(), nullptr, &cigar, &l_cigar);
    
    if (alignment.first.orientation() == edlib_alignment::status::REVERSE)
        seq = (rcseq = make_reverse_complement(query.seq1)).data();

    // mate 1
    bam_buf.add_bam(query.header.size()-2, query.header.data(), flag1, tgt_id, pos1, 255, l_cigar, cigar,
                    tgt_id, pos2, tlen1, query.seq1.size(), seq, nullptr, 0);

    if (alignment.second.aligned())
        sam_parse_cigar(alignment.second.cigar(), nullptr, &cigar, &l_cigar);
    
    seq = query.seq2.data();
    if (alignment.second.orientation() == edlib_alignment::status::REVERSE)
        seq = (rcseq = make_reverse_complement(query.seq2)).data();
    
    // mate 2
    bam_buf.add_bam(query.header.size()-2, query.header.data(), flag2, tgt_id, pos2, 255, l_cigar, cigar,
                    tgt_id, pos1, tlen2, query.seq2.size(), seq, nullptr, 0);

    free(cigar);
}
#endif



#ifdef RC_BAM
void show_bam_minimal(bam_buffer& bam_buf, const database& db, const sequence_query& query, target_id tgt, bool primary) {

    // function only applicable for mapped reads atm
    // function only applicable for paired reads atm

    size_t l_tgt = db.get_target(tgt).seq().size();
    size_t l_read = query.seq1.size();
    int64_t l_template = std::min(l_tgt, l_read);

    uint16_t flag1 = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FMREVERSE | BAM_FREAD1 | ((!primary)*BAM_FSECONDARY);
    uint16_t flag2 = BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FREVERSE | BAM_FREAD2 | ((!primary)*BAM_FSECONDARY);

    size_t n_cigar = 1;
    uint32_t cigar[2];

    cigar[0] = (l_template << 4);
    if (l_tgt < l_read) {
        cigar[1] = ((l_read-l_tgt) << 4) | 1;
        ++n_cigar;
    }

    // mate 1
    bam_buf.add_bam(query.header.size()-2, query.header.data(), flag1, tgt, 0, 255, n_cigar, cigar,
             tgt, 0, l_template, l_read, query.seq1.data(), nullptr, 0);
    
    // mate 2
    std::string recv2 = make_reverse_complement(query.seq2);

    bam_buf.add_bam(query.header.size()-2, query.header.data(), flag2, tgt, 0, 255, n_cigar, cigar,
             tgt, 0, -l_template, l_read, recv2.data(), nullptr, 0);
}
#endif

#ifdef RC_BAM
void prepare_bam(const database& db, const query_options& opt, classification_results& results) {
    std::string sam_header_text = db.get_sam_header();
    results.bamOut = sam_open(opt.samFile.data(), "wb1");
    hts_set_threads(results.bamOut, opt.performance.bamThreads);
    results.bamHdr = sam_hdr_parse(sam_header_text.size(), sam_header_text.data());
    sam_hdr_write(results.bamOut, results.bamHdr);
    // TODO handle errors
}
#endif


void show_as_alignment(mappings_buffer& buf, const database& db, 
    const query_options& opt, const sequence_query& query, 
    const classification_candidates& cands)
{
    if (opt.output.samMode == sam_mode::none || cands.empty()) // TODO allow showing of unmapped in SAM / BAM
        return;

    size_t primary = 0;
    for (size_t i = 0; i < cands.size(); ++i)
        if (cands[primary].hits < cands[i].hits)
            primary = i;

    if (opt.output.samMode == sam_mode::sam)
        for (size_t i = 0; i < cands.size(); ++i)
            show_sam_minimal(buf.align_out, db.get_target(cands[i].tgt), query, i == primary);
    
    #ifdef RC_BAM
    else if (opt.output.samMode == sam_mode::bam)
        for (size_t i = 0; i < cands.size(); ++i)
            show_bam_minimal(buf.bam_buf, db, query, cands[i].tgt, i == primary);
    #endif
}


void align_candidates(mappings_buffer& buf, const database& db, 
    const query_options& opt, const sequence_query& query, 
    classification_candidates& cands)
{
    if (cands.empty()) return; //TODO allow showing of unmapped in SAM / BAM
    
    alns_vector alns;
    size_t primary = 0;

    const auto align_candidate = [&](const auto& cand) {
        alns.emplace_back(query, cand.tgt, db, opt.classify.maxEditDist);
        if (!alns.back().aligned()) {
            alns.pop_back();
            return true;
        } else if (alns.back().score() < alns[primary].score()) {
            primary = alns.size()-1;
        }
        return false;
    };

    cands.erase(std::remove_if(cands.begin(), cands.end(), align_candidate), cands.end());

    if (opt.output.samMode == sam_mode::sam)
        for (size_t i = 0; i < alns.size(); ++i) //TODO allow showing of unmapped in SAM / BAM
            show_sam_alignment(buf.align_out, db, query, alns[i], i == primary);
    
    #ifdef RC_BAM
    else if (opt.output.samMode == sam_mode::bam)
        for (size_t i = 0; i < alns.size(); ++i) //TODO allow showing of unmapped in SAM / BAM
            show_bam_alignment(buf.bam_buf, query, alns[i], i == primary);
    #endif
}


/*************************************************************************//**
 *
 * @brief transcriptome classification scheme 2-pass variant;
 *        saves memory at expense of speed
 *
 *****************************************************************************/
void map_queries_to_targets_2pass(
    const vector<string>& infiles,
    const database& db, const query_options& opt,
    classification_results& results)
{
    matches_per_target_light coverage_;

    const auto makeCovBuffer = [] { return matches_per_target_light(); };

    const auto processCoverage = [&] (matches_per_target_light& buf,
        const sequence_query& query, const auto& allhits)
    {
        if (query.empty()) return;

        for_each_eligible_mapping(db, opt.classify, query, allhits, [&](const auto& cand) {
            buf.insert(allhits, cand, opt.classify.covFill);
        });
    };

    const auto mergeCoverage = [&] (matches_per_target_light&& buf) {
        coverage_.merge(std::move(buf));
    };

    const auto appendToOutput = [] (const std::string&) {
        // nothing because we might need an intact sam file
    };

    // 1st pass: generate coverage
    query_database(infiles, db, opt.pairing, opt.performance,
                   makeCovBuffer, processCoverage, mergeCoverage,
                   appendToOutput);
    
    if (opt.output.samMode == sam_mode::sam)
        db.show_sam_header(results.samOut);
    #ifdef RC_BAM
    else if (opt.output.samMode == sam_mode::bam)
        prepare_bam(db, opt, results);
    #endif

    const auto makeBatchBuffer = [&] {
        #ifdef RC_BAM
        if (opt.output.samMode == sam_mode::bam)
            return mappings_buffer(opt.performance.bamBufSize);
        else
        #endif
            return mappings_buffer();
    };

    const auto processQuery = [&] (mappings_buffer& buf,
        const sequence_query& query, const auto& allhits)
    {
        if(query.empty()) return;

        classification cls = classify(db, opt.classify, query, allhits, coverage_);
       
        if (opt.output.evaluate.determineGroundTruth)
            cls.groundTruth = ground_truth_target(db, query.header);

        if (opt.classify.align)
            align_candidates(buf, db, opt, query, cls.candidates); // removes unalignable candidates
        else
            show_as_alignment(buf, db, opt, query, cls.candidates);

        show_query_mapping(buf.out, db, opt.output, query, cls, allhits);
            
        evaluate_classification(opt.output.evaluate, cls, results.statistics);
    };

    const auto finalizeBatch = [&] (mappings_buffer&& buf) {
        results.mainOut << buf.out.str();
        results.samOut << buf.align_out.str();

        #ifdef RC_BAM
        if (opt.output.samMode == sam_mode::bam) {
            for (bam1_t& aln: buf.bam_buf.vec) {
                sam_write1(results.bamOut, results.bamHdr, &aln); //TODO: handle errors
                bam_destroy1(&aln);
            }
        }
        #endif
    };

    // 2nd pass: process queries
    query_database(infiles, db, opt.pairing, opt.performance,
                   makeBatchBuffer, processQuery, finalizeBatch,
                   appendToOutput);

    #ifdef RC_BAM
    if (results.bamOut) sam_close(results.bamOut);
    if (results.bamHdr) sam_hdr_destroy(results.bamHdr);
    #endif
}

/*************************************************************************//**
 *
 * @brief transcriptome classification scheme
 *
 *****************************************************************************/
void map_queries_to_targets(const vector<string>& infiles,
                            const database& db,
                            const query_options& opt,
                            classification_results& results)
{
    std::cerr << db.query_sketcher().window_stride() << " " << db.query_sketcher().window_size() << " " << db.query_sketcher().sketch_size() << '\n'
              << db.target_sketcher().window_stride() << " " << db.target_sketcher().window_size() << " " << db.target_sketcher().sketch_size() << std::endl;
    std::cerr << opt.classify.hitsMin << " " << opt.classify.hitsCutoff << " " << opt.classify.covMin << " [" << (opt.classify.covNorm==coverage_norm::max?"NORM":"NO NORM") << "]" << std::endl;

    if(opt.output.format.showMapping)
        show_query_mapping_header(results.mainOut, opt.output);
    map_queries_to_targets_2pass(infiles, db, opt, results);
}



} // namespace mc
