/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2020 André Müller (muellan@uni-mainz.de)
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

#ifndef MC_CLASSIFICATION_STATISTICS_H_
#define MC_CLASSIFICATION_STATISTICS_H_

#include <atomic>

#include "taxonomy.h"
#include "stat_confusion.h"


namespace mc {

/*************************************************************************//**
 *
 * @brief for tracking RNA classification properties
 *
 *****************************************************************************/
class rna_mapping_statistics
{
public:
    using count_t = std::uint_least64_t;

    double recall() const noexcept {
        return double(originMapped_) / totalReads_;
    }

    double hits_per_reads() const noexcept {
        return double(totalMatches_) / totalReads_;
    }

    double total_true_hit_rate() const noexcept {
        return double(originMapped_) / totalMatches_;
    }

    double mean_true_hit_rate() const noexcept {
        return double(originMappedWeighted_) / alignedReads_;
    }

    count_t total() const noexcept {
        return totalReads_;
    }

    count_t matches() const noexcept {
        return totalMatches_;
    }

    count_t correct() const noexcept {
        return originMapped_;
    }

    count_t true_negatives() const noexcept {
        return correctlyRejected_;
    }

    count_t aligned() const noexcept {
        return alignedReads_;
    }

    // count mapping
    void mapped(size_t numMatches) noexcept {
        ++totalReads_;
        if(numMatches > 0) {
            ++alignedReads_;
            totalMatches_ += numMatches;
        }
    }

    // eval: identified origin
    void mapped_correct(size_t numMatches) noexcept {
        ++totalReads_;
        totalMatches_ += numMatches;
        ++originMapped_;
        ++alignedReads_;
        double tmp = originMappedWeighted_, inc = 1.0/numMatches;
        while(!originMappedWeighted_.compare_exchange_weak(tmp, tmp+inc));
    }

    // eval: identified noise
    void mapped_noise() noexcept {
        ++totalReads_;
        ++correctlyRejected_;
    }

private:
    using atom_t = std::atomic<count_t>;
    atom_t totalReads_{0};   // total #reads
    atom_t totalMatches_{0}; // total #candidates
    atom_t alignedReads_{0}; // #reads where #cands >= 1
    
    // requires ground truth
    atom_t originMapped_{0}; // #mappings including ground truth
    std::atomic<double> originMappedWeighted_{0}; // #mappings including ground truth
    atom_t correctlyRejected_{0}; // #noise reads left unmapped

};


} // namespace mc


#endif
