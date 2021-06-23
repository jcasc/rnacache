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
#include <map>
#include <array>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <typeinfo>

#include "options.h"
#include "candidates.h"
#include "database.h"
#include "printing.h"
#include "typename.h"


namespace mc {

using std::cout;
using std::cerr;
using std::endl;
using std::string;


/*************************************************************************//**
 *
 *
 *****************************************************************************/
void show_database_config(const string& dbfile)
{
    auto db = make_database(dbfile, database::scope::metadata_only);
    print_static_properties(db);
    print_content_properties(db);
}



/*************************************************************************//**
 *
 *
 *****************************************************************************/
void print_query_config()
{
    cout << "hit classifier       "
         << type_name<classification_candidates>() << '\n'
         << "------------------------------------------------\n";
}



/*************************************************************************//**
 *
 *
 *****************************************************************************/
void show_database_statistics(const string& dbfile)
{
    auto db = make_database(dbfile);
    print_static_properties(db);
    print_content_properties(db);
}




/*************************************************************************//**
 *
 *
 *****************************************************************************/
void show_feature_map(const string& dbfile)
{
    auto db = make_database(dbfile);
    print_static_properties(db);
    print_content_properties(db);
    cout << "===================================================\n";
    db.print_feature_map(cout);
    cout << "===================================================\n";
}



/*************************************************************************//**
 *
 *
 *****************************************************************************/
void show_feature_counts(const string& dbfile)
{
    auto db = make_database(dbfile);
    print_static_properties(db);
    print_content_properties(db);
    cout << "===================================================\n";
    db.print_feature_counts(cout);
    cout << "===================================================\n";
}



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
void show_target_info(std::ostream& os, const database& db, target_id tgt)
{
    const auto& meta = db.get_target(tgt);
    os  << "Target " << meta.name() << "):\n"
        << "    source:     "
        << meta.source().filename << " / " << meta.source().index
        << "\n    length:     " << meta.source().windows << " windows"
        << "\n    (ID) Name:  " << "(" << tgt << ") " << meta.name()
        << '\n';
}



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
void show_target_info(const info_options& opt)
{
    auto db = make_database(opt.dbfile, database::scope::metadata_only);

    if(!opt.targetIds.empty()) {
        for(const auto& tid : opt.targetIds) {
            const auto tgt = db.target_with_name(tid);
            if(tgt != database::nulltgt) {
                show_target_info(cout, db, tgt);
            }
            else {
                cout << "Target (reference sequence) '" << tid
                     << "' not found in database.\n";
            }
        }
    }
    else {
        cout << "Targets (reference sequences) in database:\n";
        for(target_id tgt = 0; tgt < db.target_count(); ++tgt) {
            show_target_info(cout, db, tgt);
        }
    }
}



/*************************************************************************//**
 *
 * @brief
 *
 *****************************************************************************/
void show_basic_exec_info()
{
    database db;
    print_static_properties(db);
    print_query_config();
    cout << endl;
}



/*************************************************************************//**
 *
 * @brief shows database properties
 *
 *****************************************************************************/
void main_mode_info(const cmdline_args& args)
{
    auto opt = get_info_options(args);

    switch(opt.mode) {
        default:
        case info_mode::basic:
            show_basic_exec_info();
            break;
        case info_mode::targets:
            show_target_info(opt);
            break;
        case info_mode::db_config:
            show_database_config(opt.dbfile);
            print_query_config();
            break;
        case info_mode::db_statistics:
            show_database_statistics(opt.dbfile);
            break;
        case info_mode::db_feature_map:
            show_feature_map(opt.dbfile);
            break;
        case info_mode::db_feature_counts:
            show_feature_counts(opt.dbfile);
            break;
    }
}


} // namespace mc
