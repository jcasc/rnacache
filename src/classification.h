/******************************************************************************
 *
 * MetaCache - Meta-Genomic Classification Tool
 *
 * Copyright (C) 2016-2019 André Müller (muellan@uni-mainz.de)
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

/******************************************************************************
 *
 * @file contains the declaration of to the two main
 *       entry points for querying (read mapping)
 *
 *****************************************************************************/
#ifndef MC_CLASSIFICATION_H_
#define MC_CLASSIFICATION_H_

#include <vector>
#include <string>

#include "database.h"


namespace mc {


/// @brief forward declarations
struct query_options;
struct classification_results;


/*************************************************************************//**
 *
 * @brief try to map each read from the input files to a target
 *        according to the query options
 *
 *****************************************************************************/
void map_queries_to_targets(
    const std::vector<std::string>& inputFilenames,
    const database&, const query_options&,
    classification_results&);

} // namespace mc

#endif
