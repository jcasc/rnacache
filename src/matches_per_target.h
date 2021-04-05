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

#ifndef MC_MATCHES_PER_TARGET_H_
#define MC_MATCHES_PER_TARGET_H_


#include <vector>
#include <unordered_map>
#include <algorithm> // move, sort, lower_bound
#include <iterator>  // make_move_iterator

#include "candidates.h"
#include "options.h"

namespace mc {


/*************************************************************************//**
 *
 * @brief records matches (and their query origin) per classification target
 *
 *****************************************************************************/
class matches_per_target
{

public:
    //---------------------------------------------------------------
    using window_id        = database::window_id;
    using match_count_type = database::match_count_type;
    using location         = database::match_locations::value_type;

    struct window_matches {
        window_matches() = default;

        constexpr
        window_matches(window_id w, match_count_type c) noexcept :
            win{w}, hits{c}
        {}

        friend bool
        operator < (const window_matches& a, const window_matches& b) noexcept {
            return (a.win < b.win);
        }
        friend bool
        operator > (const window_matches& a, const window_matches& b) noexcept {
            return (a.win > b.win);
        }

        window_id win = 0;
        match_count_type hits = 0;
    };

    using matches_per_window = std::vector<window_matches>;

    struct candidate {
        candidate(query_id qid, matches_per_window mpw):
            qeryid{qid}, matches{std::move(mpw)}
        {}
        query_id qeryid = 0;
        matches_per_window matches;
    };

private:
    using hits_per_target =
        std::unordered_map<target_id, std::vector<candidate>>;

public:
    using iterator        = hits_per_target::iterator;
    using const_iterator  = hits_per_target::const_iterator;


    //-------------------------------------------------------------------
    bool empty() const noexcept {
        return hitsPerTarget_.empty();
    }

    std::size_t size() const noexcept {
        return hitsPerTarget_.size();
    }

    const_iterator find(target_id tgt) const {
        return hitsPerTarget_.find(tgt);
    }

    const_iterator begin() const noexcept { return hitsPerTarget_.begin(); }
    const_iterator end()   const noexcept { return hitsPerTarget_.end(); }

    iterator erase(const_iterator pos) { return hitsPerTarget_.erase(pos); }

    hits_per_target::size_type erase(const hits_per_target::key_type& key) {
        return hitsPerTarget_.erase(key);
    }

    //---------------------------------------------------------------
    template<class Locations>
    void insert(query_id qid,
                const Locations& matches,
                const classification_candidates& candidates,
                std::uint_least32_t minHitsPerCandidate = 0)
    {
        for(const auto& cand : candidates) {

            if(cand.tax && cand.hits >= minHitsPerCandidate) {

                auto tgt = cand.tgt;

                // find candidate in matches
                location lm{cand.pos.beg, tgt};
                auto it = std::lower_bound(matches.begin(), matches.end(), lm);
                // fill window vector
                if(it == matches.end()) return;

                // create new window vector
                matches_per_window mpw;
                mpw.reserve(cand.pos.end - cand.pos.beg + 1);

                while(it != matches.end() &&
                        it->tgt == tgt &&
                        it->win <= cand.pos.end)
                {
                    if(mpw.size() > 0 && mpw.back().win == it->win)
                        mpw.back().hits++;
                    else
                        mpw.emplace_back(it->win, 1);
                    ++it;
                }
                // insert into map
                hitsPerTarget_[tgt].emplace_back(qid, std::move(mpw));
            }
        }
    }

    //---------------------------------------------------------------
    void merge(matches_per_target&& other)
    {
        for(auto& mapping : other) {

            auto& source = mapping.second;
            auto& target = hitsPerTarget_[mapping.first];

            target.insert(target.end(),
                          std::make_move_iterator(source.begin()),
                          std::make_move_iterator(source.end()) );
        }
    }

    //---------------------------------------------------------------
    void sort_match_lists()
    {
        for(auto& mapping : hitsPerTarget_) {
            std::sort(mapping.second.begin(), mapping.second.end(),
                [] (const candidate& a, const candidate& b) {
                    if(a.matches.front() < b.matches.front()) return true;
                    if(a.matches.front() > b.matches.front()) return false;
                    if(a.matches.back() < b.matches.back()) return true;
                    if(a.matches.back() > b.matches.back()) return false;
                    return (a.qeryid < b.qeryid);
                });
        }
    }

private:
    hits_per_target hitsPerTarget_;
};


/*************************************************************************//**
 *
 * @brief records matches (and their query origin) per classification target
 *        Lightweight version of the above, no query- or candidate-level
 *        information is retained
 *
 *****************************************************************************/
class matches_per_target_light
{

public:
    //---------------------------------------------------------------
    using window_id        = database::window_id;
    using location         = database::match_locations::value_type;

private:
    using window_hits = std::unordered_set<window_id>;
    using hits_per_target =
        std::unordered_map<target_id, window_hits>;

public:
    using const_iterator  = hits_per_target::const_iterator;


    //-------------------------------------------------------------------
    bool empty() const noexcept {
        return hitsPerTarget_.empty();
    }

    std::size_t size() const noexcept {
        return hitsPerTarget_.size();
    }

    const_iterator find(target_id tgt) const {
        return hitsPerTarget_.find(tgt);
    }

    const_iterator begin() const noexcept { return hitsPerTarget_.begin(); }
    const_iterator end()   const noexcept { return hitsPerTarget_.end(); }



    //---------------------------------------------------------------
    void insert(const match_locations& matches,
                const match_candidate& cand,
                const coverage_fill fill)
    {   
        auto& hpt = hitsPerTarget_[cand.tgt];
       
        if (fill == coverage_fill::matches) {
            // find candidate in matches
            location lm{cand.pos.beg, cand.tgt};
            auto it = std::lower_bound(matches.begin(), matches.end(), lm);

            window_id cur_win = -1;
            while(it != matches.end() && it->tgt == cand.tgt && it->win <= cand.pos.end) {
                if(it->win != cur_win) {
                    hpt.emplace(it->win);
                    cur_win = it->win;
                }
                ++it;
            }
        } else { // fill == coverage:fill
            for (window_id win = cand.pos.beg; win<=cand.pos.end; ++win)
                hitsPerTarget_[cand.tgt].emplace(win);
        }
    }

    //---------------------------------------------------------------
    void merge(matches_per_target_light&& other)
    {
        for(auto& mapping : other) {
            const auto it = hitsPerTarget_.emplace(std::move(mapping));
            if (!it.second) {
                auto& source = mapping.second;
                auto& target = it.first->second;
                target.insert(source.begin(), source.end());
            }
        }
    }

    //---------------------------------------------------------------
    size_t num_hits(target_id tgt) const {
        if (!hitsPerTarget_.count(tgt))
            return 0;
        return hitsPerTarget_.at(tgt).size();
    }

    // void print(const database& db) const {
    //     std::ofstream os("covdist.tsv");
    //     for (const auto& tgt: *this)
    //         os << tgt.first << "\t" << double(num_hits(tgt.first))/db.taxon_of_target(tgt.first)->source().windows << '\n';
    // }

private:
    hits_per_target hitsPerTarget_;
};

/*************************************************************************//**
 *
 * @brief records matches (and their query origin) per classification target
 *        Experimental version of the above, used for parameter search.
 *        Saves highest thresholds with which a position
 *        would have been covered.
 *        
 *
 *****************************************************************************/
class matches_per_target_param
{

public:
    //---------------------------------------------------------------
    using window_id        = database::window_id;
    using location         = database::match_locations::value_type;

private:

    struct hit_type {
        uint16_t abs;
        double rel;

        friend bool operator<(const hit_type& lhs, const hit_type& rhs) {
            if (lhs.abs == rhs.abs)
                return lhs.rel < rhs.rel;
            return lhs.abs < rhs.abs;
        }
    };

    struct parametrized_coverage_hit {
        
        void hit(uint16_t abs, double rel) {
            hit_type h{abs, rel};
            auto it = hits_.begin();
            while (it != hits_.end() && !(*it < h)) {
                if (it->rel >= rel)
                    return;
                ++it;
            }
            if (it != hits_.end() && it->rel <= rel) {
                *it++ = h;
                auto end_del = it;
                while (end_del != hits_.end() && end_del->rel <= rel)
                    ++end_del;
                hits_.erase(it, end_del);
            } else {
                hits_.emplace(it, h);
            }
        }

        bool get(uint16_t abs, double rel) const {
            auto it = hits_.crbegin();
            while (it != hits_.crend() && it->abs < abs) {
                if (it->rel <= rel)
                    return false;
                ++it;
            }
            return it != hits_.crend() && it->rel >= rel;
        }

        void merge(parametrized_coverage_hit other) {
            for (const auto h: other.hits_)
                hit(h.abs, h.rel);
        }
        
        size_t size() const {
            return hits_.size();
        }

        private:
            std::vector<hit_type> hits_;
    };

    using window_hits = std::unordered_map<window_id, parametrized_coverage_hit>;
    using hits_per_target =
        std::unordered_map<target_id, window_hits>;
    
public:

    using const_iterator  = hits_per_target::const_iterator;

    //-------------------------------------------------------------------
    bool empty() const noexcept {
        return hitsPerTarget_.empty();
    }

    std::size_t size() const noexcept {
        return hitsPerTarget_.size();
    }

    const_iterator find(target_id tgt) const {
        return hitsPerTarget_.find(tgt);
    }

    const_iterator begin() const noexcept { return hitsPerTarget_.begin(); }
    const_iterator end()   const noexcept { return hitsPerTarget_.end(); }

    //---------------------------------------------------------------
    void insert(const match_locations& matches,
                const classification_candidates& candidates,
                const coverage_fill fill)
    {
        uint_least32_t max_hits = 0;
        for(const auto& cand: candidates) if (cand.hits > max_hits) max_hits = cand.hits;

        for(const auto& cand: candidates) {
            auto& hpt = hitsPerTarget_[cand.tgt];

            if(fill == coverage_fill::matches) {
                // find candidate in matches
                location lm{cand.pos.beg, cand.tgt};
                auto it = std::lower_bound(matches.begin(), matches.end(), lm);

                window_id cur_win = -1;
                while(it != matches.end() &&
                        it->tgt == cand.tgt &&
                        it->win <= cand.pos.end)
                {
                    if (it->win != cur_win) {
                        hpt[it->win].hit(cand.hits, double(cand.hits)/max_hits);
                        cur_win = it->win;
                    }
                    ++it;
                }
            } else { // fill == coverage_fill::fill
                for (window_id win = cand.pos.beg; win<=cand.pos.end; ++win)
                    hpt[win].hit(cand.hits, double(cand.hits)/max_hits);
            }
        }
    }

    //---------------------------------------------------------------
    void merge(matches_per_target_param&& other)
    {
        for (auto& other_seq : other) {
            const auto it_seq = hitsPerTarget_.emplace(std::move(other_seq));
            if (!it_seq.second) {
                for (auto& other_win : other_seq.second) {
                    const auto it_win = it_seq.first->second.emplace(std::move(other_win));
                    if (!it_win.second)
                        it_win.first->second.merge(other_win.second);
                }
            }
        }
    }

    //---------------------------------------------------------------
    size_t num_hits(const target_id tgt, size_t abs, double rel) const {
        if (!hitsPerTarget_.count(tgt)) return 0;
        size_t covered = 0;
        for (auto const win: hitsPerTarget_.at(tgt)) {
            covered += win.second.get(abs, rel);
        }
        return covered;
    }

private:
    hits_per_target hitsPerTarget_;

};

} // namespace mc

#endif
