///////////////////////////////////////////////////////////////////////////////
// maxprotein.hh
//
// Compute the set of foods that maximizes protein, within a calorie budget,
// with the greedy method or exhaustive search.
//
///////////////////////////////////////////////////////////////////////////////


#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

// Simple structure for a single protein
struct Protein {
	Protein() {
		description = "";
		sequence = "";
	}
	Protein(std::string desc, std::string seq) {
		description = desc;
		sequence = seq;
	}
	std::string		description;
	std::string 	sequence;
};

// Alias for a vector of shared pointers to Protein objects.
typedef std::vector<std::shared_ptr<Protein>> ProteinVector;


// -------------------------------------------------------------------------
// Load all the proteins from a standard FASTA format file with one line
// per sequence (multi-line sequences are not allowed).
// Returns false on I/O error.
// -------------------------------------------------------------------------
bool load_proteins(ProteinVector & proteins, const std::string& path)
{
	//std::cout << "Loading proteins from [" << path << "]" << std::endl;
	proteins.clear();
	std::ifstream ifs(path.c_str());
	if (!ifs.is_open() || !ifs.good()) {
		std::cout << "Failed to open [" << path << "]" << std::endl;
		return false;
	}
	int proteinsLoaded = 0;
	bool have_description = false;
	std::shared_ptr<Protein> newProtein = nullptr;
	while (!ifs.eof()) {
		std::string lineBuffer;
		std::getline(ifs, lineBuffer);
		if (ifs.eof()) {
			break;
		}
		if (lineBuffer.size() == 0) {
			continue;
		}
		if (lineBuffer[0] == '>') {
			newProtein = std::shared_ptr<Protein>(new Protein);
			newProtein->description = lineBuffer.substr(1);
			have_description = true;
		}
		else if (have_description) {
			newProtein->sequence = lineBuffer;
			proteins.push_back(newProtein);
			proteinsLoaded++;
			have_description = false;
		}
	}

	ifs.close();
	//std::cout << "Loaded " << proteinsLoaded << " proteins from [" << path << "]" << std::endl;

	return true;
}


// -------------------------------------------------------------------------
int dynamicprogramming_longest_common_subsequence(const std::string & string1,
	const std::string & string2)
{
	const int n = string1.size();

	const int m = string2.size();

	int** d = new int*[n + 1];

	for (int i = 0; i < n + 1; i++) {
		d[i] = new int[m + 1];
	}

	int up;
	int left;
	int diag;

	for (int i = 0; i <= n; i++) {
		d[i][0] = 0;
	}
	for (int j = 0; j <= m; j++) {
		d[0][j] = 0;
	}
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= m; j++) {
			up = d[i - 1][j];
			left = d[i][j - 1];
			diag = d[i - 1][j - 1];
			if (string1[i - 1] == string2[j - 1]) {
				diag = diag + 1;

			}
			d[i][j] = std::max(std::max(up, left), diag);
		}

	}
	int result = d[n][m];
	for (int i = 0; i < n + 1; i++) {
		delete[] d[i];
	}
	delete[] d;
	return result;
}

// -------------------------------------------------------------------------
std::unique_ptr<std::vector<std::string>> generate_all_subsequences(const std::string & sequence)
{

	std::unique_ptr<std::vector<std::string>> subsequences(new std::vector<std::string>);
	double n = pow(2, sequence.size());
	unsigned int i, j;
	std::string subsequence;

	for (i = 0; i < n; i++) {
		subsequence = "";
		for (j = 0; j < sequence.size(); j++) {
			if (((i >> j) & 1) == 1) {
				subsequence += sequence[j];
			}
		}
		subsequences->push_back(subsequence);
	}





	return subsequences;
}


// -------------------------------------------------------------------------
int exhaustive_longest_common_subsequence(const std::string & string1,
	const std::string & string2)
{
	std::unique_ptr<std::vector<std::string>> seq1;
	std::unique_ptr<std::vector<std::string>> seq2;
	seq1 = generate_all_subsequences(string1);
	seq2 = generate_all_subsequences(string2);
	unsigned int best_score = 0;
	for (std::string s1 : *seq1) {
		for (std::string s2 : *seq2) {
			if (s1 == s2 && s1.size() > best_score) {
				best_score = s1.size();


			}
		}

	}
	return best_score;
}


// -------------------------------------------------------------------------
std::shared_ptr<Protein> exhaustive_best_match(ProteinVector & proteins, const std::string & string1)
{
	int best_i = 0;
	int best_score = 0;
	int score;
	for (unsigned int i = 0; i < proteins.size(); i++) {
		score = exhaustive_longest_common_subsequence(proteins[i]->sequence, string1);
		if (score > best_score) {
			best_score = score;
			best_i = i;

		}
	}
	return proteins[best_i];
}

// -------------------------------------------------------------------------
std::shared_ptr<Protein> dynamicprogramming_best_match(ProteinVector & proteins, const std::string & string1)
{
	int best_i = 0;
	int best_score = 0;
	int score;
	for (unsigned int i = 0; i < proteins.size(); i++) {
		score = dynamicprogramming_longest_common_subsequence(proteins[i]->sequence, string1);
		if (score > best_score) {
			best_score = score;
			best_i = i;

		}
	}
	return proteins[best_i];
	
}


