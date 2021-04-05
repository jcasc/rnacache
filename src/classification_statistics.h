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
 * @brief useful for tracking classification properties (precision, ...)
 *
 *****************************************************************************/
class taxon_mapping_statistics
{
public:
    //---------------------------------------------------------------
    using rank = taxonomy::rank;
    using count_t = std::uint_least64_t;


    //---------------------------------------------------------------
    taxon_mapping_statistics() noexcept :
        assigned_{}, known_{}, correct_{}, wrong_{}, coverage_{}
    {
        for(auto& x : assigned_) x = 0;
        for(auto& x : known_)    x = 0;
        for(auto& x : correct_)  x = 0;
        for(auto& x : wrong_)    x = 0;
    }

    taxon_mapping_statistics(const taxon_mapping_statistics&) = delete;
    taxon_mapping_statistics(taxon_mapping_statistics&&) = delete;

    taxon_mapping_statistics& operator = (const taxon_mapping_statistics&) = delete;
    taxon_mapping_statistics& operator = (taxon_mapping_statistics&&) = delete;


    //---------------------------------------------------------------
    /**
     * @brief    counts one assignment
     * @details  concurrency-safe
     */
    void assign(rank assigned) noexcept
    {
        if(assigned == rank::none) {
            ++assigned_[int(rank::none)];
        }
        else {
            for(rank r = assigned; r <= rank::root; ++r) ++assigned_[int(r)];
        }
    }

    //---------------------------------------------------------------
    /**
     * @brief   counts one assignment including ground truth and correctness
     *          assessment
     *
     * @details concurrency-safe
     *
     * @param assigned : lowest rank of current assignment
     * @param known    : lowest rank for which ground truth was known
     * @param correct  : lowest rank for which current assignment is correct
     */
    void assign_known_correct(rank assigned, rank known, rank correct) noexcept
    {
        assign(assigned);

        //plausibility check
        if(correct < assigned) correct = assigned;
        if(correct < known)    correct = known;

        //if ground truth known -> count correct and wrong assignments
        if(known == rank::none) {
            ++known_[int(rank::none)];
        }
        else {
            for(rank r = known; r <= rank::root; ++r) ++known_[int(r)];

            if(correct == rank::none) {
                ++correct_[int(rank::none)];
            }
            else {
                for(rank r = correct; r <= rank::root; ++r) ++correct_[int(r)];
            }
            //if ranks below the correct rank are known and assigned,
            //then all ranks below the correct rank are wrong
            if(correct > known && correct > assigned) {
                for(rank r = rank::Sequence; r < correct; ++r) {
                    ++wrong_[int(r)];
                }
            }
        }
    }


    //---------------------------------------------------------------
    /// @details concurrency-safe
    void count_coverage_true_pos(rank r) noexcept {
        coverage_[int(r)].count_true_pos();
    }
    /// @details concurrency-safe
    void count_coverage_false_pos(rank r) noexcept {
        coverage_[int(r)].count_false_pos();
    }
    /// @details concurrency-safe
    void count_coverage_true_neg(rank r) noexcept {
        coverage_[int(r)].count_true_neg();
    }
    /// @details concurrency-safe
    void count_coverage_false_neg(rank r) noexcept {
        coverage_[int(r)].count_false_neg();
    }

    //-----------------------------------------------------
    const confusion_statistics&
    coverage(rank r) const noexcept {
        return coverage_[int(r)];
    }


    count_t assigned() const noexcept {
        return assigned_[int(rank::root)];
    }
    /**
     * @brief number of assignments on a taxonomix rank (and above)
     */
    count_t assigned(rank r) const noexcept {
        return assigned_[int(r)];
    }
    count_t unassigned() const noexcept {
        return assigned_[int(rank::none)];
    }

    //---------------------------------------------------------------
    count_t total() const noexcept {
        return assigned() + unassigned();
    }

    count_t known() const noexcept {
        return known_[int(rank::root)];
    }
    /**
     * @brief number of cases with ground truth known on ranks >= r
     */
    count_t known(rank r) const noexcept {
        return known_[int(r)];
    }
    count_t unknown() const noexcept {
        return known_[int(rank::none)];
    }

    /**
     * @brief number of known correct assignments on ranks >= r
     */
    count_t correct(rank r) const noexcept {
        return correct_[int(r)];
    }
    count_t correct() const noexcept {
        return correct_[int(rank::root)];
    }

    /**
     * @brief number of known wrong assignments on ranks <= r
     */
    count_t wrong(rank r) const noexcept {
        return wrong_[int(r)];
    }
    count_t wrong() const noexcept {
        return wrong_[int(rank::root)];
    }


    //---------------------------------------------------------------
    double known_rate(rank r) const noexcept {
        return total() > 0 ?  known(r) / double(total()) : 0;
    }
    double known_rate() const noexcept {
        return total() > 0 ?  known() / double(total()) : 0;
    }
    double unknown_rate() const noexcept {
        return total() > 0 ?  unknown() / double(total()) : 0;
    }
    double classification_rate(rank r) const noexcept {
        return total() > 0 ? assigned(r) / double(total()) : 0;
    }
    double unclassified_rate() const noexcept {
        return total() > 0 ? unassigned() / double(total()) : 0;
    }

    double sensitivity(rank r) const noexcept {
        return known(r) > 0 ? correct(r) / double(known(r)) : 0;
    }
    double precision(rank r) const noexcept {
        //note that in general tot != assigned(r) and tot != known(r)
        double tot = correct(r) + wrong(r);
        return tot > 0 ? correct(r) / tot : 0;
    }


private:
    //---------------------------------------------------------------
    std::atomic<count_t> assigned_[taxonomy::num_ranks+1];
    std::atomic<count_t> known_[taxonomy::num_ranks+1];
    std::atomic<count_t> correct_[taxonomy::num_ranks+1];
    std::atomic<count_t> wrong_[taxonomy::num_ranks+1];
    confusion_statistics coverage_[taxonomy::num_ranks+1];
};






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
