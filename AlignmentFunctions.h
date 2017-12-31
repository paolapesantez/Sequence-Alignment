#ifndef ALIGNMENTFUNCTIONS_H_INCLUDED
#define ALIGNMENTFUNCTIONS_H_INCLUDED

#include "AlignmentComponents.h"

class AlignmentFunctions
{
	public:
		AlignmentFunctions(int ma, int mi, int h, int g, int option);
		~AlignmentFunctions();
		string* correctSequences(string s1, string s2);
		DP_cell** initializeMatricesGlobalAlignment(int m, int n, int option_matrix);
		DP_cell** initializeMatricesLocalAlignment(int m, int n, int option_matrix);
		int substitutionEdit(char a, char b);	
		bool findSubstring(string sequence1, string sequence2);
		int findMaximumVector(int*previos_scores, int number_items);
		int* findMaximumMatrix(DP_cell** matrix, int m, int n);
		void printMatrix(DP_cell** matrix, int m, int n);
		DP_cell** copyMatrix(DP_cell** matrix, int m, int n);
		vector<DP_cell**> globalAlignAffineGapPenalty(string sequence1, string sequence2);
		vector<DP_cell**> localAlignAffineGapPenalty(string sequence1, string sequence2);		
		void traceBackGlobalAlignment(string sequence1, string sequence2, DP_cell** S,DP_cell** D,DP_cell** I);
		void traceBackLocalAlignment(string sequence1, string sequence2, DP_cell** S,DP_cell** D,DP_cell** I,  int position_option);
		void printOutputGlobalAlignment(string s1,string sm, string s2);
		void printOutputLocalAlignment(string s1,string sm, string s2, int initial_pos_s1, int initial_pos_s2, int position_option);
		void generateOutputFile(string text);

		int position_best_local_alignments[4];

	private:		
		string name_sequence1;
		string name_sequence2;
		int match;
		int mismatch;
		int h_gap;
		int g_gap;
		int option_alignment;
		int maximum_score;		
		
};
#endif
