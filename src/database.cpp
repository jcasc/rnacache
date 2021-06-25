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

#include "database.h"


namespace mc {

// ----------------------------------------------------------------------------
bool database::add_target(const sequence& seq, target_name sid,
                          file_source source)
{
    using std::begin;
    using std::end;

    //reached hard limit for number of targets
    if(targets_.size() >= max_target_count()) {
        throw target_limit_exceeded_error{};
    }

    if(seq.empty()) return false;

    //don't allow non-unique sequence ids
    if(name2tax_.find(sid) != name2tax_.end()) return false;

    const auto targetCount = target_id(targets_.size());

    //sketch sequence -> insert features
    source.windows = add_all_window_sketches(seq, targetCount);

    //store sequence metadata
    targets_.emplace_back(sid, std::move(source));

    //allows lookup via sequence id
    name2tax_.insert({std::move(sid), targetCount});

    return true;
}



// ----------------------------------------------------------------------------
void database::read(const std::string& filename, scope what)

{
    std::ifstream is{filename, std::ios::in | std::ios::binary};

    if(!is.good()) {
        throw file_access_error{"can't open file " + filename};
    }

    //database version info
    using std::uint64_t;
    using std::uint8_t;
    uint64_t dbVer = 0;
    read_binary(is, dbVer);

    if(uint64_t( RC_DB_VERSION ) != dbVer) {
        throw file_read_error{
            "Database " + filename + " (version " + std::to_string(dbVer) + ")"
            + " is incompatible\nwith this version of RNACache"
            + " (uses version " + std::to_string(RC_DB_VERSION) + ")" };
    }

    //data type info
    {
        //data type widths
        uint8_t featureSize = 0; read_binary(is, featureSize);
        uint8_t targetSize = 0;  read_binary(is, targetSize);
        uint8_t windowSize = 0;  read_binary(is, windowSize);
        uint8_t bucketSize = 0;  read_binary(is, bucketSize);
        uint8_t tgtidSize = 0;   read_binary(is, tgtidSize);

        if( (sizeof(feature) != featureSize) ||
            (sizeof(target_id) != targetSize) ||
            (sizeof(bucket_size_type) != bucketSize) ||
            (sizeof(window_id) != windowSize) ||
            (sizeof(target_id) != tgtidSize) )
        {
            throw file_read_error{
                "Database " + filename +
                " is incompatible with this variant of RNACache" +
                " due to different data type sizes"};
        }
    }

    clear();

    //sketching parameters
    read_binary(is, targetSketcher_);
    read_binary(is, querySketcher_);

    //target insertion parameters
    read_binary(is, maxLocsPerFeature_);

    //target metadata
    read_binary(is, targets_);

    //sequence id lookup
    name2tax_.clear();
    for(target_id t = 0; t < targets_.size(); ++t) {
        name2tax_.insert({targets_[t].name(), t});
    }

    if(what == scope::metadata_only) return;

    //hash table
    read_binary(is, features_);

    if (what == scope::everything)
        reread_targets();
}



// ----------------------------------------------------------------------------
void database::write(const std::string& filename) const
{
    using std::uint64_t;
    using std::uint8_t;

    std::ofstream os{filename, std::ios::out | std::ios::binary};

    if(!os.good()) {
        throw file_access_error{"can't open file " + filename};
    }

    //database version info
    write_binary(os, uint64_t( RC_DB_VERSION ));

    //data type widths
    write_binary(os, uint8_t(sizeof(feature)));
    write_binary(os, uint8_t(sizeof(target_id)));
    write_binary(os, uint8_t(sizeof(window_id)));
    write_binary(os, uint8_t(sizeof(bucket_size_type)));
    write_binary(os, uint8_t(sizeof(target_id)));

    //sketching parameters
    write_binary(os, targetSketcher_);
    write_binary(os, querySketcher_);

    //target insertion parameters
    write_binary(os, maxLocsPerFeature_);

    //target metadata
    write_binary(os, targets_);

    //hash table
    write_binary(os, features_);
}



// ----------------------------------------------------------------------------
void database::max_locations_per_feature(bucket_size_type n)
{
    if(n < 1) n = 1;
    if(n >= max_supported_locations_per_feature()) {
        n = max_supported_locations_per_feature();
    }
    else if(n < maxLocsPerFeature_) {
        for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
            if(i->size() > n) features_.shrink(i, n);
        }
    }
    maxLocsPerFeature_ = n;
}



// ----------------------------------------------------------------------------
database::feature_count_type
database::remove_features_with_more_locations_than(bucket_size_type n)
{
    //note that features are not really removed, because the hashmap
    //does not support erasing keys; instead all values belonging to
    //the key are cleared and the key is kept without values
    feature_count_type rem = 0;
    for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
        if(i->size() > n) {
            features_.clear(i);
            ++rem;
        }
    }
    return rem;
}



// ----------------------------------------------------------------------------
database::feature_count_type
database::remove_ambiguous_features(bucket_size_type maxambig)
{
    feature_count_type rem = 0;

    if(maxambig == 0) maxambig = 1;

    for(auto i = features_.begin(), e = features_.end(); i != e; ++i) {
        if(!i->empty()) {
            std::set<target_id> targets;
            for(auto loc : *i) {
                targets.insert(loc.tgt);
                if(targets.size() > maxambig) {
                    features_.clear(i);
                    ++rem;
                    break;
                }
            }
        }
    }
    return rem;
}



// ----------------------------------------------------------------------------
void database::clear() {
    targets_.clear();
    name2tax_.clear();
    features_.clear();
}


// ----------------------------------------------------------------------------
/**
 * @brief very dangerous! clears feature map without memory deallocation
 */
void database::clear_without_deallocation() {
    targets_.clear();
    name2tax_.clear();
    features_.clear_without_deallocation();
}




// ----------------------------------------------------------------------------
database
make_database(const std::string& filename, database::scope what, info_level info)
{
    if(filename.empty()) throw file_access_error{"No database name given"};

    database db;

    const bool showInfo = info != info_level::silent;

    if(showInfo) {
        std::cerr << "Reading database from file '"
                  << filename << "' ... " << std::flush;
    }
    try {
        db.read(filename, what);
        if(showInfo) std::cerr << "done." << std::endl;
    }
    catch(const file_access_error& e) {
        std::cerr << "FAIL" << std::endl;
        throw file_access_error{"Could not read database file '" + filename + "'"};
    }

    return db;
}


} // namespace mc
