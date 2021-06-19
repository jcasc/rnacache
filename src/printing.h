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

#ifndef MC_PRINT_RESULTS_H_
#define MC_PRINT_RESULTS_H_

#include <string>
#include <iosfwd>

#include "config.h"
#include "taxonomy.h"
#include "classify_common.h"
#include "classification_statistics.h"
#include "matches_per_target.h"


namespace mc {

// forward declaration
class database;
class classification_output_formatting;


/*************************************************************************//**
 *
 * @brief prints classification parameters
 *
 *****************************************************************************/
void show_query_parameters(std::ostream&, const query_options&);


/*************************************************************************//**
 *
 * @brief prints taxon information according to output options
 *
 *****************************************************************************/
void show_taxon(std::ostream&,
                const database&,
                const classification_output_formatting&,
                const taxon* classified);


/*************************************************************************//**
 *
 * @brief prints header for taxon information
 *
 *****************************************************************************/
void show_taxon_header(std::ostream&,
                       const classification_output_formatting&,
                       const std::string& prefix = "");


/*************************************************************************//**
 *
 * @brief prints top classification candidates
 *
 *****************************************************************************/
void show_candidates(std::ostream&,
                     const classification_candidates&);


/*************************************************************************//**
 *
 * @brief prints target.window hit matches
 *
 *****************************************************************************/
template<class Locations>
void show_matches(std::ostream&,
                  const database&,
                  const Locations&);


/*************************************************************************//**
 *
 * @brief prints target.window hit statistics from database
 *
 *****************************************************************************/
void show_candidate_ranges(std::ostream&,
                           const database&,
                           const classification_candidates&);


/*************************************************************************//**
 *
 * @brief prints a list of all match locations for each classification target
 *
 *****************************************************************************/
void show_matches_per_targets(std::ostream&,
                              const database&,
                              const matches_per_target&,
                              const classification_output_formatting&);

/*************************************************************************//**
 *
 * @brief prints RNA mapping statistics
 *
 *****************************************************************************/
void show_statistics(std::ostream&,
                     const rna_mapping_statistics&,
                     const std::string& prefix = "");


/*************************************************************************//**
 *
 * @brief show summary and statistics of classification
 *
 *****************************************************************************/
void show_summary(const query_options& opt,
                  const classification_results& results);


/*************************************************************************//**
 *
 * @brief prints database properties to stdout
 *
 *****************************************************************************/
void print_static_properties(const database&);


/*************************************************************************//**
 *
 * @brief prints database properties to stdout
 *
 *****************************************************************************/
void print_content_properties(const database&);


} // namespace mc


#endif
