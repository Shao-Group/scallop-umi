/*
Part of Berth
(c) 2024 by Xiaofei Carl Zang, Mingfu Shao, and Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __ALIGNER_H__
#define __ALIGNER_H__

#include "util.h"

using namespace std;

int subseq_pos(const string& seq1, const string& seq2, int nm, int indel_penalty, int mis_penalty);

string revcomp(const string& s);
char   revcomp_char(const char&   c);


/**
 * @param sequences Variable number of sequences to validate
 * @return true if all sequences contain only ATCG, false otherwise
 */
template<typename... Strings>
bool validate_dna_seq(const Strings&... sequences) 
{
    auto check_sequence = [](const string& seq) -> bool 
    {
        static const set<char> dnachars {'A', 'T', 'C', 'G'};
        for (char c : seq) 
        {
            if (dnachars.find(c) == dnachars.end()) return false;        
        }
        return true;
    };
    return (... && check_sequence(sequences));
}

#endif
