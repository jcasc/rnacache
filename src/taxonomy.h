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

#ifndef MC_TAXONOMY_H_
#define MC_TAXONOMY_H_

#include <array>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <cassert>

#include "io_serialize.h"


namespace mc {


/*************************************************************************//**
 *
 * @brief id-based taxonomy
 *
 * @details the directed graph is stored implicitly
 *
 *****************************************************************************/
class taxonomy
{
public:
    //---------------------------------------------------------------
    using taxon_id   = std::int_least64_t;
    using taxon_name = std::string;

    static constexpr taxon_id none_id() noexcept { return 0; }

    /************************************************************
     *
     * @brief taxonomic node
     *
     ************************************************************/
    class taxon {
        friend class taxonomy;
    public:
        //-----------------------------------------------------
        struct file_source {
            using index_t   = std::uint_least64_t;
            using window_id = std::uint_least64_t;

            explicit
            file_source(std::string filename = "", index_t index = 0,
                        window_id numWindows = 0)
            :
                filename{std::move(filename)}, windows{numWindows}, index{index}
            {}

            std::string filename;
            window_id windows;
            index_t index;
        };

        //default: empty taxon
        explicit
        taxon(taxon_id taxonId = none_id(),
              std::string taxonName = "--",
              file_source source = file_source{})
        :
            id_{taxonId},
            name_{std::move(taxonName)},
            source_{std::move(source)}
        {}

        //makes taxa sortable
        inline friend bool
        operator < (const taxon& a, const taxon& b) noexcept {
            return a.id_ < b.id_;
        }

        taxon_id id() const noexcept { return id_; }

        const taxon_name& name() const noexcept { return name_; }

        const file_source& source() const noexcept { return source_; }

        //-----------------------------------------------------
        friend
        void read_binary(std::istream& is, taxon& t) {
            read_binary(is, t.id_);
            read_binary(is, t.name_);
            read_binary(is, t.source_.filename);
            read_binary(is, t.source_.index);
            read_binary(is, t.source_.windows);
        }

        //-----------------------------------------------------
        friend
        void write_binary(std::ostream& os, const taxon& t) {
            write_binary(os, t.id_);
            write_binary(os, t.name_);
            write_binary(os, t.source_.filename);
            write_binary(os, t.source_.index);
            write_binary(os, t.source_.windows);
        }

        //-----------------------------------------------------
    private:
        taxon_id id_;
        taxon_name name_;
        file_source source_;
    };


private:
    using taxon_store = std::set<taxon>;

public:
    //-----------------------------------------------------
    using file_source    = taxon::file_source;
    using size_type      = taxon_store::size_type;
    using const_iterator = taxon_store::const_iterator;
    //-----------------------------------------------------
    class const_range {
    public:
        explicit
        const_range(const_iterator beg, const_iterator end):
            beg_{beg}, end_{end}
        {}
        const_iterator begin() const { return beg_; }
        const_iterator end() const   { return end_; }
    private:
        const_iterator beg_;
        const_iterator end_;
    };


    //-------------------------------------------------------------------
    void clear() {
        taxa_.clear();
    }


    //-----------------------------------------------------
    const_iterator
    emplace(taxon_id taxonId,
            std::string taxonName = "",
            file_source source = file_source{})
    {
        if(taxonId == none_id()) return taxa_.end();

        return taxa_.emplace(taxonId,
                             std::move(taxonName),
                             std::move(source)).first;
    }


    //---------------------------------------------------------------
    const_iterator
    insert(const taxon& t) {
        return taxa_.insert(t).first;
    }
    //-----------------------------------------------------
    const_iterator
    insert(taxon&& t) {
        return taxa_.insert(std::move(t)).first;
    }


    //---------------------------------------------------------------
    const_iterator
    insert_or_replace(const taxon& t)
    {
        auto i = taxa_.find(taxon{t.id()});
        if(i != taxa_.end()) {
            taxa_.erase(i);
            return taxa_.insert(t).first;
        } else {
            return taxa_.insert(t).first;
        }
    }
    //-----------------------------------------------------
    const_iterator
    insert_or_replace(taxon&& t)
    {
        auto i = taxa_.find(taxon{t.id()});
        if(i != taxa_.end()) {
            taxa_.erase(i);
            return taxa_.insert(std::move(t)).first;
        } else {
            return taxa_.insert(t).first;
        }
    }

    //---------------------------------------------------------------
    const_iterator
    find(taxon_id id) const {
        return taxa_.find(taxon{id});
    }
    /**
     * @return taxon to 'id' or nullptr;
     */
    const taxon*
    operator [] (taxon_id id) const {
        auto i = find(id);
        if(i == taxa_.end()) return nullptr;
        return &(*i);
    }


    //---------------------------------------------------------------
    bool
    contains(taxon_id id) const {
        return find(id) != taxa_.end();
    }



    //---------------------------------------------------------------
    size_type size() const noexcept {
        return taxa_.size();
    }
    //-----------------------------------------------------
    bool empty() const noexcept {
        return taxa_.empty();
    }


    //---------------------------------------------------------------
    const_iterator begin() const { return taxa_.begin(); }
    const_iterator end()   const { return taxa_.end(); }
    //-----------------------------------------------------
    const_range
    full_range() const {
        return const_range{taxa_.begin(), taxa_.end()};
    }
    //-----------------------------------------------------
    const_range
    subrange(taxon_id first, taxon_id last) const {
        return const_range{taxa_.lower_bound(taxon{first}),
                           taxa_.upper_bound(taxon{last})};
    }
    //-----------------------------------------------------
    const_range
    subrange_from(taxon_id first) const {
        return const_range{taxa_.lower_bound(taxon{first}), taxa_.end()};
    }
    //-----------------------------------------------------
    const_range
    subrange_until(taxon_id last) const {
        return const_range{taxa_.begin(), taxa_.upper_bound(taxon{last})};
    }


    //---------------------------------------------------------------
    friend void
    read_binary(std::istream& is, taxonomy& tax)
    {
        tax.taxa_.clear();
        std::uint64_t nt = 0;
        read_binary(is, nt);
        for(std::uint64_t i = 0; i < nt; ++i) {
            taxon t;
            read_binary(is, t);
            tax.taxa_.insert(std::move(t));
        }
    }

    //-----------------------------------------------------
    friend void
    write_binary(std::ostream& os, const taxonomy& tax)
    {
        write_binary(os, std::uint64_t(tax.taxa_.size()));
        for(const auto& t : tax.taxa_) {
            write_binary(os, t);
        }
    }


private:
    //---------------------------------------------------------------
    taxon_store taxa_;
};




/*************************************************************************//**
 *
 * @brief pull some types from taxonomy into global namespace
 *
 *****************************************************************************/
using taxon          = taxonomy::taxon;
using taxon_id       = taxonomy::taxon_id;


} // namespace mc


#endif
