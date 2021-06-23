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
    const database& db,
    const classification_evaluation_options& opt,
    const sequence_query& query,
    classification& cls,
    rna_mapping_statistics& statistics)
{
    if (opt.determineGroundTruth)
        cls.groundTruth = ground_truth_target(db, query.header);
    
    if (opt.accuracy) {
        if (!cls.groundTruth) {
            if (cls.candidates.empty())
                return statistics.mapped_noise();
        } else
            for (const auto& cand: cls.candidates)
                if (cand.tgt == cls.groundTruth)
                    return statistics.mapped_correct(cls.candidates.size());
        statistics.mapped(cls.candidates.size());
    } else {
        statistics.mapped();
    }
}

#ifdef RC_BAM
struct bam_buffer
{
    bam_vector vec;
    std::unique_ptr<uint8_t[]> buf;
    bool full = 0;
    size_t used = 0;
    size_t N;

    bam_buffer(size_t n) : buf(new uint8_t[n]), N(n) {}

    template<typename ...Arg>
    void add_bam(Arg&& ...arg) {
        vec.emplace_back();
        bam1_t& b = vec.back();
        bam_set_mempolicy(&b, BAM_USER_OWNS_STRUCT);
        if (!full) {
            bam_set_mempolicy(&b, BAM_USER_OWNS_STRUCT | BAM_USER_OWNS_DATA);
            b.data = &buf[used];
            b.m_data = N-used;
        }

        int ret = bam_set1(&b, std::forward<Arg>(arg)...);
        
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
    mappings_buffer(size_t n) : bam_buf(n) {}
    #endif
};

struct alignment_targets {
    
    struct target {
        std::string header;
        std::string seq;
    };

    alignment_targets() = default;

    alignment_targets(const database& db) {
        load(db);
    }

    void load(const database& db) {
        using indexed_targets = std::unordered_map<size_t, target_id>;
        using catalogue = std::pair<std::unique_ptr<sequence_reader>, indexed_targets>;
        
        std::unordered_map<std::string, catalogue> catalogues;

        for (target_id tgt = 0; tgt < db.target_count(); ++tgt) {
            const auto& src = db.get_target(tgt).source();
            if (!catalogues.count(src.filename))
                catalogues.emplace(src.filename, catalogue{make_sequence_reader(src.filename), indexed_targets()});
            catalogues[src.filename].second.emplace(src.index, tgt);
        }

        targets_.resize(db.target_count());
        for (auto& i: catalogues) {
            auto& cat = i.second;
            auto& reader = cat.first;
            auto& targets = cat.second;
            while (reader->has_next()) {
                auto idx = reader->index();
                if(targets.count(idx)) {
                    auto seq = reader->next();
                    targets_[targets[idx]] = {std::move(seq.header), std::move(seq.data)};
                } else {
                    reader->skip(1);
                }
            }
        }
        // debug();
    }

    const target& operator[](target_id tgt) const {
        return targets_[tgt];
    }

    auto begin() const {return targets_.cbegin();};
    auto end() const {return targets_.cend();}
    size_t size() const {return targets_.size();};
    
    alignment_targets(const alignment_targets&) = delete;
    alignment_targets& operator=(const alignment_targets&) = delete;

    void show_sam_header(std::ostream& os) const {
        os << "@HD\tVN:1.0 SO:unsorted\n";
        for(const auto& tgt: targets_)
            os << "@SQ\tSN:" << tgt.header << "\tLN:" << tgt.seq.size() << '\n';
        os << "@PG\tID:rnaache\tPN:rnaache\tVN:" << RC_VERSION_STRING << '\n';
    }

    std::string get_sam_header() const {
        std::ostringstream tmp;
        show_sam_header(tmp);
        return tmp.str();
    }

private:
    std::vector<target> targets_;
};

/*************************************************************************//**
 *
 * @brief edlib alignment container
 *
 *****************************************************************************/
struct edlib_alignment {

    enum struct status {FORWARD, REVERSE, UNMAPPED};

    edlib_alignment(const std::string& query, target_id tgt, const alignment_targets& refs, int max_edit_distance):
        tgt_(tgt), status_(status::UNMAPPED), cigar_(nullptr)
    {
        const std::string& target = refs[tgt].seq;

        auto edlib_config = edlibNewAlignConfig(max_edit_distance, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities.data(), additionalEqualities.size());
        
        auto regular = edlibAlign(query.c_str(), query.size(), target.c_str(), target.size(), edlib_config);
        std::string reverse_query = make_reverse_complement(query);
        auto reverse_complement = edlibAlign(reverse_query.c_str(), reverse_query.size(), target.c_str(), target.size(), edlib_config);
        
        if (regular.status != EDLIB_STATUS_OK || regular.status != EDLIB_STATUS_OK) {
            throw std::runtime_error{"edlib failed!"};
        }

        // keep in mind ../dev/edlip.cpp:212
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
    
    edlib_alignment(edlib_alignment&& other):
        tgt_(other.tgt_), status_(other.status_), score_(other.score_), start_(other.start_), end_(other.end_), cigar_(other.cigar_)
    {
        other.cigar_ = nullptr;
    }

    ~edlib_alignment() {free(cigar_);}

    bool mapped() const noexcept {return status_!=status::UNMAPPED;}
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
    edlib_alignment_pair(const sequence_query& query, target_id tgt, const alignment_targets& refs, int max_edit_distance):
        first(query.seq1, tgt, refs, max_edit_distance),
        second(query.seq2, tgt, refs, max_edit_distance)
    {}

    int score() const {
        return first.score() + second.score();
    }

    edlib_alignment first, second;
};

void show_sam_alignment(std::ostream& os,
                        const alignment_targets& refs,
                        const sequence_query& query,
                        const edlib_alignment_pair& alignment,
                        bool primary = true)
{
    if (!alignment.first.mapped() && !alignment.second.mapped()) return;

    uint_least16_t flag1 = 0, flag2 = 0;
    flag1 |= 0x1;
    flag2 |= 0x1;
    if (alignment.first.mapped() && alignment.second.mapped()) {
        flag1 |= 0x2;
        flag2 |= 0x2;
    }
    if (!alignment.first.mapped()) {
        flag1 |= 0x4;
        flag2 |= 0x8;
    }
    if (!alignment.second.mapped()) {
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

    int tlen = 0;
    if (alignment.first.mapped() && alignment.second.mapped()) 
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
    os << (alignment.first.mapped() ? refs[alignment.first.tgt()].header : refs[alignment.second.tgt()].header)  << '\t';

    // POS
    os << (alignment.first.mapped() ? alignment.first.start() : alignment.second.start()) + 1  << '\t';

    // MAPQ
    os << 255 << '\t';

    // CIGAR
    os << (alignment.first.mapped() ? alignment.first.cigar() : "*")  << '\t';

    // RNEXT
    os << "=" << '\t';

    // PNEXT
    os << (alignment.second.mapped() ? alignment.second.start() : alignment.first.start())  + 1 << '\t';

    // TLEN
    os << (alignment.first.start() <= alignment.second.start() ? tlen: -tlen) << '\t';

    // SEQ
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
    os << (alignment.second.mapped() ? refs[alignment.second.tgt()].header : refs[alignment.first.tgt()].header)  << '\t';

    // POS
    os << (alignment.second.mapped() ? alignment.second.start() : alignment.first.start()) + 1  << '\t';

    // MAPQ
    os << 255 << '\t';

    // CIGAR
    os << (alignment.second.mapped() ? alignment.second.cigar() : "*")  << '\t';

    // RNEXT
    os << "=" << '\t';

    // PNEXT
    os << (alignment.first.mapped() ? alignment.first.start() : alignment.second.start())  + 1 << '\t';

    // TLEN
    os << (alignment.second.start() <= alignment.first.start() ? tlen: -tlen) << '\t';

    // SEQ
    os << query.seq2 << '\t';

    // QUAL
    os << "*" << '\t';

    os << '\n';
}

void show_sam_minimal(std::ostream& os, const alignment_targets& refs, const sequence_query& query, target_id tgt, bool primary) {

    // function only applicable for mapped reads

    std::string qname = query.header.substr(0, query.header.size() - 2);
    size_t tgtlen = refs[tgt].seq.size();
    size_t readlen = query.seq1.size();
    size_t tlen = std::min(tgtlen, readlen);

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
    os << refs[tgt].header << '\t';

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
    os << refs[tgt].header << '\t';

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
    os << tlen << '\t';

    // SEQ
    auto recv = query.seq2;
    reverse_complement(recv);
    os << recv << '\t';

    // QUAL
    os << "*";

    os << '\n';
}

#ifdef RC_BAM
void add_bam_minimal(bam_buffer& bam_buf, const alignment_targets& refs, const sequence_query& query, target_id tgt, bool primary) {

    // function only applicable for mapped reads atm
    // function only applicable for paired reads atm

    std::string qname = query.header.substr(0, query.header.size() - 2);
    size_t l_tgt = refs[tgt].seq.size();
    size_t l_read = query.seq1.size();
    size_t l_template = std::min(l_tgt, l_read);

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
    bam_buf.add_bam(qname.size(), qname.data(), flag1, tgt, 0, 255, n_cigar, cigar,
             tgt, 0, l_template, l_read, query.seq1.data(), nullptr, 0);
    
    // mate 2
    std::string recv2 = query.seq2;
    reverse_complement(recv2);

    bam_buf.add_bam(qname.size(), qname.data(), flag2, tgt, 0, 255, n_cigar, cigar,
             tgt, 0, l_template, l_read, recv2.data(), nullptr, 0);
}
#endif

/*************************************************************************//**
 *
 * @brief show alignment(s) in SAM format (no header)
 *
 *****************************************************************************/
void show_sam(std::ostream& os,
              const query_options& opt,
              const sequence_query& query,
              const classification_candidates& cands,
              const alignment_targets& refs)
{
    if (opt.rnaMode == rna_mode::sam_align) {

        // find all valid alignments for this query
        std::vector<edlib_alignment_pair> alignments;
        size_t primary = 0;
        for (const auto& cand: cands) {

            alignments.emplace_back(query, cand.tgt, refs, opt.classify.maxEditDist);
            if (alignments[primary].score() < alignments.back().score())
                primary = alignments.size()-1;
        }

        if (alignments.empty()) return;

        // show primary alignment
        show_sam_alignment(os, refs, query, alignments[primary]);

        // show secondary alignments
        for (size_t i = 0; i < alignments.size(); ++i) {
            if (i == primary) continue;
            show_sam_alignment(os, refs, query, alignments[i], false);
        }
    }
    
    else if (opt.rnaMode == rna_mode::sam) {
        size_t max_cand = 0;
        for (size_t i = 0; i < cands.size(); ++i) {
            const auto& cand = cands[i];

            if (cand.hits > cands[max_cand].hits) {
                max_cand = i;
            }
        }
        for (size_t i = 0; i < cands.size(); ++i) {
            const auto& cand = cands[i];

            show_sam_minimal(os, refs, query, cand.tgt, i == max_cand);
        }
    }
}


#ifdef RC_BAM
void add_bam_alignments(bam_buffer& bam_buf,
                        const query_options& opt,
                        const sequence_query& query,
                        const classification_candidates& cands,
                        const alignment_targets& refs)
{
    if (opt.rnaMode != rna_mode::bam) return;
    
    size_t max_cand = 0;
    for (size_t i = 0; i < cands.size(); ++i) {
        const auto& cand = cands[i];
        if (cand.hits > cands[max_cand].hits) {
            max_cand = i;
        }
    }
    for (size_t i = 0; i < cands.size(); ++i) {
        const auto& cand = cands[i];
        add_bam_minimal(bam_buf, refs, query, cand.tgt, i == max_cand);
    }


}
#endif




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

    const auto finalizeCoverage = [&] (matches_per_target_light&& buf) {
        coverage_.merge(std::move(buf));
    };

    //runs if something needs to be appended to the output
    const auto appendToOutput = [] (const std::string&) {
        // nothing because we might need an intact sam file
    };

    // 1st pass: generate coverage
    query_database(infiles, db, opt.pairing, opt.performance,
                   makeCovBuffer, processCoverage, finalizeCoverage,
                   appendToOutput);
    
    // reread all references
    
    // show sam header
    alignment_targets refs;
    
    #ifdef RC_BAM
    samFile   *bam_file   = nullptr;
    sam_hdr_t *bam_header = nullptr;

    #endif

    if (opt.rnaMode == rna_mode::sam || opt.rnaMode == rna_mode::sam_align) {
        refs.load(db);
        refs.show_sam_header(results.alignmentOut);
    }
    
    #ifdef RC_BAM
    else if (opt.rnaMode == rna_mode::bam) {

        refs.load(db);
        std::string sam_header_text = refs.get_sam_header();
        bam_file = sam_open("mc_debug.bam", "wb1");
        hts_set_threads(bam_file, opt.performance.bamThreads);
        bam_header = sam_hdr_parse(sam_header_text.size(), sam_header_text.data());
        if (!bam_header)
            std::cerr << "sam_hdr_parse failed!" << std::endl;
        if (sam_hdr_write(bam_file, bam_header) < 0)
            std::cerr << "sam_hdr_write failed!" << std::endl;

    }
    #endif

    const auto makeBatchBuffer = [&] {
        #ifdef RC_BAM
        return mappings_buffer(opt.performance.bamBufSize);
        #else
        return mappings_buffer();
        #endif
    };

    // updates buffer with the database answer of a single query
    const auto processQuery = [&] (mappings_buffer& buf,
        const sequence_query& query, const auto& allhits)
    {
        if(query.empty()) return;

        auto cls = classify(db, opt.classify, query, allhits, coverage_);

        evaluate_classification(db, opt.output.evaluate, query, cls, results.statistics);

        show_query_mapping(buf.out, db, opt.output, query, cls, allhits);
        
        show_sam(buf.align_out, opt, query, cls.candidates, refs);
        
        #ifdef RC_BAM
        add_bam_alignments(buf.bam_buf, opt, query, cls.candidates, refs);
        #endif

    };

    // size_t asdf = 0;
    const auto finalizeBatch = [&] (mappings_buffer&& buf) {
        results.perReadOut << buf.out.str();
        results.alignmentOut << buf.align_out.str();

        #ifdef RC_BAM
        if (opt.rnaMode == rna_mode::bam) {
            // if (buf.bam_buf.full)
                // std::cerr << "BAM buf status = " << (buf.bam_buf.full?"full":"free") << std::endl;
            for (bam1_t& aln: buf.bam_buf.vec) {
                // std::cerr << "ASDFFFFF FINALIZE" << std::endl;
                if (sam_write1(bam_file, bam_header, &aln) < 0)
                    std::cerr << "SAM_WRITE1 ERROR" << std::endl;
                bam_destroy1(&aln);
            }
        }
        #endif
        // asdf+=opt.performance.batchSize;
        // std::cerr << asdf << std::endl;

    };

    // 2nd pass: process queries
    query_database(infiles, db, opt.pairing, opt.performance,
                   makeBatchBuffer, processQuery, finalizeBatch,
                   appendToOutput);

    #ifdef RC_BAM
    if (bam_header) sam_hdr_destroy(bam_header);
    if (bam_file) sam_close(bam_file);
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
    
    if(opt.rnaMode == rna_mode::search) {
        // pass
    } else {
        if(opt.output.format.mapViewMode != map_view_mode::none)
            show_query_mapping_header(results.perReadOut, opt.output);
        map_queries_to_targets_2pass(infiles, db, opt, results);
    }
}



} // namespace mc
