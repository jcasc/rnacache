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

#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "options.h"
#include "filesys_utility.h"


namespace mc {


//-------------------------------------------------------------------
void main_mode_help(const cmdline_args& args)
{
    if(args.size() < 3 || args[1] != "help" || args[2] == "help") {

        if(args.size() > 1 && args[1] != "help") {
            std::cerr << "ERROR: Invalid command line arguments!\n\n";
        }
        else {
            std::cout <<
                "RNACache Copyright (C) 2016-2021 Julian Cascitti & André Müller & Robin Kobus\n"
                "This program comes with ABSOLUTELY NO WARRANTY.\n"
                "This is free software, and you are welcome to redistribute it\n"
                "under certain conditions. See the file 'LICENSE' for details.\n\n";
        }

        std::cout <<
            "USAGE:\n"
            "\n"
            "    rnacache <MODE> [OPTION...]\n"
            "\n"
            "    Available modes:\n"
            "\n"
            "    help        shows documentation \n"
            "    query       classify read sequences using pre-built database\n"
            "    build       build new database from reference sequences\n"
            "    modify      add reference sequences\n"
            "    info        show database and reference sequence properties\n"
            "\n"
            "\n"
            "EXAMPLES:\n"
            "\n"
            "    Query single FASTA file 'myreads.fna' against pre-built database 'refseq':\n"
            "        rnacache query refseq myreads.fna -out results.txt\n"
            "    same with output to the console:\n"
            "        rnacache query refseq myreads.fna\n"
            "\n"
            "    Query all sequence files in folder 'test' againgst database 'refseq':\n"
            "        rnacache query refseq test -out results.txt\n"
            "\n"
            "    Query paired-end reads in separate files:\n"
            "        rnacache query refseq reads1.fa reads2.fa -pairfiles -out results.txt\n"
            "\n"
            "    Query paired-end reads in one file (a1,a2,b1,b2,...):\n"
            "        rnacache query refseq paired_reads.fa -pairseq -out results.txt\n"
            "    \n"
            "    View documentation for query mode:\n"
            "        rnacache help query\n"
            "\n"
            "    View documentation on how to build databases:\n"
            "        rnacache help build\n";
    }
    else if(args[2] == "build") {
        std::cout << build_mode_docs() << '\n';
    }
    else if(args[2] == "modify") {
        std::cout << modify_mode_docs() << '\n';
    }
    else if(args[2] == "query") {
        std::cout << query_mode_docs() << '\n';
    }
    else if(args[2] == "info") {
        std::cout << info_mode_docs() << '\n';
    }
    else {
        std::cerr
            << "You need to specify a mode for which to show help :\n"
            << "    " << args[0] << " help <mode>\n\n"
            << "Unknown mode '" << args[2] << "'\n\n"
            << "Available modes are:\n"
            << "    build\n"
            << "    modify\n"
            << "    query\n"
            << "    info\n";
    }
}


} // namespace mc
