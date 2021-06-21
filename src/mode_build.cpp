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

#include <cstdint>
#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <thread>
#include <chrono>

#include "timer.h"
#include "options.h"
#include "cmdline_utility.h"
#include "filesys_utility.h"
#include "io_error.h"
#include "io_options.h"
#include "database.h"
#include "printing.h"
#include "sequence_io.h"

#include "batch_processing.h"


namespace mc {

using std::string;
using std::cout;
using std::cerr;
using std::flush;
using std::endl;



// ---------------------------------------------------------------------------
struct input_sequence {
    sequence_reader::header_type header;
    sequence_reader::data_type data;
    database::file_source fileSource;
};

using input_batch = std::vector<input_sequence>;



/*************************************************************************//**
 *
 * @brief add batch of reference targets to database
 *
 * @return false, if database not ready / insertion error
 *
 *****************************************************************************/
void add_targets_to_database(
    database& db,
    const input_batch& batch,
    info_level infoLvl = info_level::moderate)
{
    for(const auto& seq : batch) {
        if(!seq.data.empty()) {
            auto seqId = extract_accession_string(
                             seq.header, sequence_id_type::any);

            // make sure sequence id is not empty,
            // use entire header if neccessary
            if(seqId.empty()) seqId = seq.header;

            if(infoLvl == info_level::verbose) {
                cout << "[" << seqId;
                cout << "] ";
            }

            // try to add to database
            bool added = db.add_target(
                seq.data, seqId, seq.fileSource);

            if(infoLvl == info_level::verbose && !added) {
                cout << seqId << " not added to database" << endl;
            }
        }
        if(db.add_target_failed()) break;
    }
}



/*************************************************************************//**
 *
 * @brief adds reference sequences from *several* files to database
 *
 *****************************************************************************/
void add_targets_to_database(database& db,
    const std::vector<string>& infiles,
    info_level infoLvl = info_level::moderate)
{
    int n = infiles.size();
    int i = 0;

    // make executor that runs database insertion (concurrently) in batches
    // IMPORTANT: do not use more than one worker thread!
    batch_processing_options execOpt;
    execOpt.batch_size(8);
    execOpt.queue_size(4);
    execOpt.concurrency(1);

    execOpt.abort_if([&] { return db.add_target_failed(); });

    execOpt.on_error([&] (std::exception& e) {
        if(dynamic_cast<database::target_limit_exceeded_error*>(&e)) {
            cout << endl;
            cerr << "Reached maximum number of targets per database ("
                 << db.max_target_count() << ").\n"
                 << "See 'README.md' on how to compile RNACache with "
                 << "support for databases with more reference targets.\n";
        }
        else if(infoLvl == info_level::verbose) {
            cout << "FAIL: " << e.what() << endl;
        }
    });

    execOpt.on_work_done([&] { db.wait_until_add_target_complete(); });

    batch_executor<input_sequence> executor { execOpt,
        [&] (int, const auto& batch) {
            add_targets_to_database(db, batch, infoLvl);
        }};

    // read sequences in main thread
    for(const auto& filename : infiles) {
        if(infoLvl == info_level::verbose) {
            cout << "  " << filename << " ... " << flush;
        } else if(infoLvl != info_level::silent) {
            show_progress_indicator(cout, i/float(n));
        }

        try {
            const auto fileId = extract_accession_string(
                                    filename, sequence_id_type::acc_ver);

            auto reader = make_sequence_reader(filename);

            while(reader->has_next() && executor.valid()) {
                // get (ref to) next input sequence storage and fill it
                auto& seq = executor.next_item();
                seq.fileSource.filename = filename;
                seq.fileSource.index = reader->index();
                reader->next_header_and_data(seq.header, seq.data);
            }

            if(infoLvl == info_level::verbose) {
                cout << "done." << endl;
            }
        }
        catch(std::exception& e) {
            if(infoLvl == info_level::verbose) {
                cout << "FAIL: " << e.what() << endl;
            }
        }
        ++i;
    }
}



/*************************************************************************//**
 *
 * @brief prepares datbase for build
 *
 *****************************************************************************/
void prepare_database(database& db, const build_options& opt)
{
    const auto dbconf = opt.dbconfig;
    if(dbconf.maxLocationsPerFeature > 0) {
        db.max_locations_per_feature(dbconf.maxLocationsPerFeature);
        cerr << "Max locations per feature set to "
             << dbconf.maxLocationsPerFeature << '\n';
    }

    if(dbconf.maxLoadFactor > 0.4 && dbconf.maxLoadFactor < 0.99) {
        db.max_load_factor(dbconf.maxLoadFactor);
        cerr << "Using custom hash table load factor of "
             << dbconf.maxLoadFactor << '\n';
    }

    if(dbconf.removeAmbigFeatures &&
       opt.infoLevel != info_level::silent)
    {
        cerr << "Ambiguous features will be removed afterwards.\n";
    }
}



/*************************************************************************//**
 *
 * @brief database features post-processing
 *
 *****************************************************************************/
void post_process_features(database& db, const build_options& opt)
{
    const bool notSilent = opt.infoLevel != info_level::silent;

    const auto& dbconf = opt.dbconfig;

    if(dbconf.removeOverpopulatedFeatures) {
        auto old = db.feature_count();
        auto maxlpf = db.max_locations_per_feature() - 1;
        if(maxlpf > 0) { //always keep buckets with size 1
            if(notSilent) {
                cout << "\nRemoving features with more than "
                     << maxlpf << " locations... " << flush;
            }
            auto rem = db.remove_features_with_more_locations_than(maxlpf);

            if(notSilent) {
                cout << rem << " of " << old << " removed." << endl;
                if(rem != old) print_content_properties(db);
            }
        }
    }

    if(dbconf.removeAmbigFeatures)
    {
        if(notSilent) {
            cout << "\nRemoving ambiguous features... " << flush;
        }

        auto old = db.feature_count();
        auto rem = db.remove_ambiguous_features(dbconf.maxTaxaPerFeature);

        if(notSilent) {
            cout << rem << " of " << old << "." << endl;
            if(rem != old) print_content_properties(db);
        }

    }
}



/*************************************************************************//**
 *
 * @brief prepares datbase for build, adds targets and writes database to disk
 *
 *****************************************************************************/
void add_to_database(database& db, const build_options& opt)
{
    prepare_database(db, opt);

    const bool notSilent = opt.infoLevel != info_level::silent;
    if(notSilent) print_static_properties(db);

    timer time;
    time.start();

    if(!opt.infiles.empty()) {
        const auto initNumTargets = db.target_count();

        if(notSilent) cout << "Processing reference sequences." << endl;

        add_targets_to_database(db, opt.infiles, opt.infoLevel);

        if(notSilent) {
            clear_current_line(cout);
            cout << "Added "
                 << (db.target_count() - initNumTargets) << " reference sequences "
                 << "in " << time.seconds() << " s" << endl;

            print_content_properties(db);
        }
    }

    post_process_features(db, opt);

    if(notSilent) {
        cout << "Writing database to file '" << opt.dbfile << "' ... " << flush;
    }
    try {
        db.write(opt.dbfile);
        if(notSilent) cout << "done." << endl;
    }
    catch(const file_access_error&) {
        if(notSilent) cout << "FAIL" << endl;
        cerr << "Could not write database file!\n";
    }

    time.stop();

    if(notSilent) {
        cout << "Total build time: " << time.seconds() << " s" << endl;
    }

    //prevents slow deallocation
    db.clear_without_deallocation();
}



/*************************************************************************//**
 *
 * @brief adds reference sequences to an existing database
 *
 *****************************************************************************/
void main_mode_modify(const cmdline_args& args)
{
    auto opt = get_modify_options(args);

    cout << "Modify database " << opt.dbfile << endl;

    auto db = make_database(opt.dbfile);

    if(opt.infoLevel != info_level::silent && !opt.infiles.empty()) {
        cout << "Adding reference sequences to database..." << endl;
    }

    add_to_database(db, opt);
}



/*************************************************************************//**
 *
 * @brief builds a database from reference input sequences
 *
 *****************************************************************************/
void main_mode_build(const cmdline_args& args)
{
    auto opt = get_build_options(args);

    if(opt.infoLevel != info_level::silent) {
        cout << "Building new database '" << opt.dbfile
             << "' from reference sequences." << endl;
    }

    //configure sketching scheme
    auto sketcher = database::sketcher{};
    sketcher.kmer_size(opt.sketching.kmerlen);
    sketcher.sketch_size(opt.sketching.sketchlen);
    sketcher.window_size(opt.sketching.winlen);
    sketcher.window_stride(opt.sketching.winstride);

    auto db = database{sketcher};

    add_to_database(db, opt);
}


} // namespace mc
