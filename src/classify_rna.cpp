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
#include "sam.h"
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
        if(cand.tax && cand.hits >= opt.hitsMin
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
                return !(cand.tax && cand.hits >= opt.hitsMin && double(cand.hits)/maxHits >= opt.hitsCutoff);
            }),
    cands.end());
}


/*************************************************************************//**
 *
 * @brief applies coverage filter to a list of taxa
 *
 *****************************************************************************/
void coverage_filter(const classification_options& opt,
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
                             / cand.tax->source().windows;
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
                         / cand.tax->source().windows;
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

    coverage_filter(opt, cls.candidates, cov);

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
    if(opt.determineGroundTruth) {
        cls.groundTruth = ground_truth_taxon(db, query.header);
        
        if (!cls.groundTruth) {
            if (cls.candidates.empty())
                return statistics.mapped_noise();
        } else
            for (const auto& cand: cls.candidates)
                if (cand.tax == cls.groundTruth)
                    return statistics.mapped_correct(cls.candidates.size());
    }
    statistics.mapped(cls.candidates.size());
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

    // query_matches mappings;
    // taxon_list truth;

    // search
    matches_per_target_param pcoverage;
    // std::vector<std::vector<uint16_t>> hitcounts;

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
        using indexed_taxons = std::unordered_map<size_t, const taxon *>;
        using catalogue = std::pair<std::unique_ptr<sequence_reader>, indexed_taxons>;
        
        std::unordered_map<std::string, catalogue> catalogues;

        for (const auto& tax: db.target_taxa()) {
            const auto& src = tax.source();
            if (!catalogues.count(src.filename))
                catalogues.emplace(src.filename, catalogue{make_sequence_reader(src.filename), indexed_taxons()});
            catalogues[src.filename].second.emplace(src.index, &tax);
            // std::cerr << "#!@ "<<src.index << " " << &tax << std::endl;
        }

        for (const auto& i: catalogues) {
            const auto& cat = i.second;
            const auto& reader = cat.first;
            const auto& taxons = cat.second;
            while (reader->has_next()) {
                auto idx = reader->index();
                if(taxons.count(idx)) {
                    auto seq = reader->next();
                    // std::cerr << "2#!@ " << idx << " " << taxons.at(idx) << std::endl;
                    tax_to_refid_.emplace(taxons.at(idx), targets_.size());
                    targets_.push_back({std::move(seq.header), std::move(seq.data)});
                } else {
                    // std::cerr << "SKIPPEDI " << idx << std::endl;
                    reader->skip(1);
                }
            }
        }
        // debug();
    }

    const target& operator()(const taxon* tax) const {
        return targets_.at(tax_to_refid_.at(tax));
    }

    uint32_t id(const taxon* tax) const {
        return tax_to_refid_.at(tax);
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
        os << "@PG\tID:metacache\tPN:metacache\tVN:" << MC_VERSION_STRING << '\n';
    }

    std::string get_sam_header() const {
        std::ostringstream tmp;
        show_sam_header(tmp);
        return tmp.str();
    }

private:
    std::unordered_map<const taxon*, uint32_t> tax_to_refid_;
    std::vector<target> targets_;
};

/*************************************************************************//**
 *
 * @brief edlib alignment container
 *
 *****************************************************************************/
struct edlib_alignment {

    enum struct status {FORWARD, REVERSE, UNMAPPED};

    edlib_alignment(const std::string& query, const taxon* tax, const alignment_targets& refs, int max_edit_distance):
        tax_(tax), status_(status::UNMAPPED), cigar_(nullptr)
    {
        const std::string& target = refs(tax).seq;

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
        tax_(other.tax_), status_(other.status_), score_(other.score_), start_(other.start_), end_(other.end_), cigar_(other.cigar_)
    {
        other.cigar_ = nullptr;
    }

    ~edlib_alignment() {free(cigar_);}

    bool mapped() const noexcept {return status_!=status::UNMAPPED;}
    int score() const noexcept {return score_;}
    status orientation() const noexcept {return status_;}
    const taxon * tax() const noexcept {return tax_;}
    int start() const noexcept {return start_;}
    int end() const noexcept {return end_;}
    const char * cigar() const noexcept {return cigar_;}


private:
    static std::vector<EdlibEqualityPair> additionalEqualities;
    const taxon* tax_;
    status status_;
    int score_;
    int start_, end_;
    char * cigar_;
};

std::vector<EdlibEqualityPair> edlib_alignment::additionalEqualities({{'a', 'A'}, {'t', 'T'}, {'c', 'C'}, {'g', 'G'}});

struct edlib_alignment_pair {
    edlib_alignment_pair(const sequence_query& query, const taxon* tax, const alignment_targets& refs, int max_edit_distance):
        first(query.seq1, tax, refs, max_edit_distance),
        second(query.seq2, tax, refs, max_edit_distance)
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
    os << (alignment.first.mapped() ? refs(alignment.first.tax()).header : refs(alignment.second.tax()).header)  << '\t';

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
    os << (alignment.second.mapped() ? refs(alignment.second.tax()).header : refs(alignment.first.tax()).header)  << '\t';

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

void show_sam_minimal(std::ostream& os, const alignment_targets& refs, const sequence_query& query, const taxon* tax, bool primary) {

    // function only applicable for mapped reads

    std::string qname = query.header.substr(0, query.header.size() - 2);
    size_t tgt = refs(tax).seq.size();
    size_t read = query.seq1.size();
    size_t tlen = std::min(tgt, read);

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
    os << refs(tax).header << '\t';

    // POS
    os << 1 << '\t';

    // MAPQ
    os << 255 << '\t';

    // CIGAR
    if (tgt < read) {
        os << tgt << "M" << (read-tgt) << "I\t";
    } else {
        os << read << "M" << "\t";
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
    os << refs(tax).header << '\t';

    // POS
    os << 1  << '\t';

    // MAPQ
    os << 255 << '\t';

    // CIGAR
    if (tgt < read) {
        os << tgt << "M" << (read-tgt) << "I\t";
    } else {
        os << read << "M" << "\t";
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
void add_bam_minimal(bam_buffer& bam_buf, const alignment_targets& refs, const sequence_query& query, const taxon* tax, bool primary) {

    // function only applicable for mapped reads atm
    // function only applicable for paired reads atm

    std::string qname = query.header.substr(0, query.header.size() - 2);
    size_t l_tgt = refs(tax).seq.size();
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
    bam_buf.add_bam(qname.size(), qname.data(), flag1, refs.id(tax), 0, 255, n_cigar, cigar,
             refs.id(tax), 0, l_template, l_read, query.seq1.data(), nullptr, 0);
    
    // mate 2
    std::string recv2 = query.seq2;
    reverse_complement(recv2);

    bam_buf.add_bam(qname.size(), qname.data(), flag2, refs.id(tax), 0, 255, n_cigar, cigar,
             refs.id(tax), 0, l_template, l_read, recv2.data(), nullptr, 0);
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
            // still neccessary??? //TODO: check
            if(!cand.tax || cand.tax->rank() != taxon_rank::Sequence)
                continue;
                
            alignments.emplace_back(query, cand.tax, refs, opt.classify.maxEditDist);
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
            if(!cand.tax || cand.tax->rank() != taxon_rank::Sequence)
                continue;
            if (cand.hits > cands[max_cand].hits) {
                max_cand = i;
            }
        }
        for (size_t i = 0; i < cands.size(); ++i) {
            const auto& cand = cands[i];
            if(!cand.tax || cand.tax->rank() != taxon_rank::Sequence)
                continue;
            show_sam_minimal(os, refs, query, cand.tax, i == max_cand);
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
        if(!cand.tax || cand.tax->rank() != taxon_rank::Sequence)
            continue;
        if (cand.hits > cands[max_cand].hits) {
            max_cand = i;
        }
    }
    for (size_t i = 0; i < cands.size(); ++i) {
        const auto& cand = cands[i];
        if(!cand.tax || cand.tax->rank() != taxon_rank::Sequence)
            continue;
        add_bam_minimal(bam_buf, refs, query, cand.tax, i == max_cand);
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


#ifndef RC_BAM
// parameter search !
//
// they don't call it high *readability* computing
//
void map_queries_to_targets_rna_search_2pass(
    const vector<string>& infiles,
    const database& db, const query_options& opt)
{    
    struct eval_result {
        // parameters
        size_t hit_thresh = 0;
        double cutoff = 0;
        double cov_stat = 0;
        bool cov_maxnorm = false;

        // results
        size_t reads_total = 0;

        size_t matches_total = 0; 
        size_t reads_aligned = 0;
        size_t origin_mapped = 0;
        size_t correctly_rejected = 0;
    };  
    struct cmp {
        bool operator() (const eval_result& lhs, const eval_result& rhs) const noexcept{
            if (lhs.matches_total == rhs.matches_total)
                return lhs.origin_mapped > rhs.origin_mapped;
            else
                return lhs.matches_total < rhs.matches_total;
        }
    };
    std::multiset<eval_result, cmp> global_results;
    std::atomic<size_t> query_counter{0};
    matches_per_target_param pcoverage_;

    std::vector<std::pair<uint16_t, double>> settings;
    
    for (int rel = 0; rel <=100; rel+=1)
        for (int abs = 1; abs <=64; abs*=2)
            settings.emplace_back(abs, double(rel)/100);
    // settings.emplace_back(opt.classify.hitsMin, opt.classify.hitsCutoff);

    
    constexpr size_t N = 100;

    std::unordered_map<size_t, std::atomic<size_t>> total_matches, origin_mapped, aligned, correctly_rejected;
    std::unordered_map<size_t, std::atomic<size_t>> total_matches_maxnorm, origin_mapped_maxnorm, aligned_maxnorm, correctly_rejected_maxnorm;
    for (size_t i = 0; i < settings.size()*(N+1); ++i) {
        total_matches_maxnorm.emplace(i,0);
        origin_mapped_maxnorm.emplace(i,0);

        total_matches.emplace(i,0);
        origin_mapped.emplace(i,0);
        aligned.emplace(i,0);
        correctly_rejected.emplace(i,0);
    }

    for (size_t i = 0; i < settings.size(); ++i) { 
        correctly_rejected_maxnorm.emplace(i,0);
        aligned_maxnorm.emplace(i,0);
    }



    std::vector<std::unordered_map<const taxon*, std::atomic<double>>> cov_cache;
    for (size_t i = 0; i < settings.size(); ++i) {
        cov_cache.emplace_back();
        for (const taxon& tax: db.target_taxa())
            cov_cache.back().emplace(&tax, -1.0);
    }

    const auto makeBatchBuffer = [] { return mappings_buffer(); };

    const auto processCoverage = [&] (mappings_buffer& buf,
        const sequence_query& query, const auto& allhits)
    {
        if(query.empty()) return;
        ++query_counter;

        buf.pcoverage.insert(allhits, make_classification_candidates(db, opt.classify, query, allhits), opt.classify.covFill);
    };

    //runs before a batch buffer is discarded
    const auto finalizeCoverage = [&] (mappings_buffer&& buf) {
        pcoverage_.merge(std::move(buf.pcoverage));
        std::cerr << query_counter << std::endl;
    };

    //runs if something needs to be appended to the output
    const auto appendToOutput = [&] (const std::string&) {
    };

    std::cerr << "STARTING FIRST PASS" << std::endl;
    //run (parallel) database queries according to processing options
    query_database(infiles, db, opt.pairing, opt.performance,
                   makeBatchBuffer, processCoverage, finalizeCoverage,
                   appendToOutput);
    
    
    query_counter = 0;

    const auto processQuery = [&] (mappings_buffer&,
        const sequence_query& query, const auto& allhits)
    {
        if(query.empty()) return;
        ++query_counter;
                
        classification cls { make_classification_candidates(db, opt.classify, query, allhits) } ;
        cls.groundTruth = ground_truth_taxon(db, query.header);

        for (size_t queue_ = 0; queue_ < settings.size(); ++queue_) {
            
            uint16_t abs = settings[queue_].first;
            double rel = settings[queue_].second;

            auto cands = cls.candidates;
            auto filter_ops = opt.classify;
            filter_ops.hitsMin = abs;
            filter_ops.hitsCutoff = rel;
            hits_cutoff_filter(filter_ops, cands);

            size_t idx = (N+1) * queue_;
            
            if (cands.empty()) {
                if (cls.groundTruth == nullptr) {
                    correctly_rejected_maxnorm[queue_]++;
                    correctly_rejected[idx]++;
                }
                continue;
            }
            
            uint_least64_t max_hits = 0;
            for (const auto& cand: cands) max_hits = std::max(cand.hits, max_hits);

            double max_cov = 0;
            for (const auto& cand: cands) {
                if (cov_cache[queue_][cand.tax] == -1.0)
                    cov_cache[queue_][cand.tax] = double(pcoverage_.num_hits(cand.tgt, abs, rel))/cand.tax->source().windows;
                max_cov = std::max(max_cov, cov_cache[queue_][cand.tax].load());
            }

            for (const auto& cand: cands) {
                
                double cov_stat = cov_cache[queue_][cand.tax]/max_cov;

                total_matches_maxnorm[idx+N* cov_stat]++;
                total_matches[idx+N* cov_cache[queue_][cand.tax]]++; 

                if (cls.groundTruth == cand.tax) {
                    origin_mapped_maxnorm[idx+N* cov_stat]++;
                    origin_mapped[idx+N* cov_cache[queue_][cand.tax]]++;
                }
            }
            
            aligned_maxnorm[queue_]++;
            aligned[idx+N* max_cov]++;

            if (cls.groundTruth == nullptr) {
                size_t reject_idx = idx+N*max_cov +1;
                if (reject_idx <= idx+N)
                    correctly_rejected[reject_idx]++;
            }
        }
    };

    //runs before a batch buffer is discarded
    const auto finalizeBatch = [&] (mappings_buffer&&) {
        std::cerr << query_counter << std::endl;
    };

    query_database(infiles, db, opt.pairing, opt.performance,
                   makeBatchBuffer, processQuery, finalizeBatch,
                   appendToOutput);

    // prefix_sums

    do_in_parallel(opt.performance.numThreads, settings.size(), [&](size_t q) {
        size_t idx = (N+1) * q;
        for (size_t i = N-1; i < N; --i) {
            total_matches_maxnorm.at(idx+ i) += total_matches_maxnorm.at(idx+ i+1);
            origin_mapped_maxnorm[idx+ i] += origin_mapped_maxnorm[idx+ i+1];

            total_matches[idx+ i] += total_matches[idx+ i+1];
            origin_mapped[idx+ i] += origin_mapped[idx+ i+1];
            aligned[idx+ i] += aligned[idx+ i+1];
            correctly_rejected[idx+ N-i] += correctly_rejected[idx+ N-i-1];
        }
    });

    for (size_t q = 0; q < settings.size(); ++q) {
        size_t idx = (N+1) * q;
        for (size_t i = 0; i <= N; ++i) {
            // filter here if neccessary

            uint16_t abs = settings[q].first;
            double rel = settings[q].second;
            global_results.insert({abs, rel, double(i)/N, true, query_counter, total_matches_maxnorm[idx+ i], aligned_maxnorm[q], origin_mapped_maxnorm[idx+ i], correctly_rejected_maxnorm[q]});
            global_results.insert({abs, rel, double(i)/N, false, query_counter, total_matches[idx+ i], aligned[idx+ i], origin_mapped[idx+ i], correctly_rejected[idx+ i]});
        }
    }

    std::cerr << global_results.size() << " results" << std::endl;
    std::cout << "hit_threshold,hit_cutoff,relative_coverage,reads,matches,aligned,TP,TN" << std::endl;
    for (const auto& i: global_results) {
        std::cout << i.hit_thresh << "," << i.cutoff << "," << i.cov_stat << "," << i.cov_maxnorm << ","
                << i.reads_total << "," << i.matches_total << "," << i.reads_aligned << "," << i.origin_mapped << "," << i.correctly_rejected << std::endl;
    }
    std::cout << std::endl;
}
#endif

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
        #ifndef RC_BAM
        map_queries_to_targets_rna_search_2pass(infiles, db, opt);
        #endif
    } else {
        if(opt.output.format.mapViewMode != map_view_mode::none)
            show_query_mapping_header(results.perReadOut, opt.output);
        map_queries_to_targets_2pass(infiles, db, opt, results);
    }
}



/*************************************************************************//**
 *
 * @brief needed for 'merge' mode: default classification scheme & output
 *        try to map candidates to a taxon with the lowest possible rank
 *
 *****************************************************************************/
void map_candidates_to_targets(
    const vector<string>&,
    const vector<classification_candidates>&,
    const database&, const query_options&,
    classification_results&)
{
    throw std::runtime_error{"merge mode not available for RNA mapping"};
}

} // namespace mc
