#ifndef ALIGNMENTCOMPONENTS_H_INCLUDED
#define ALIGNMENTCOMPONENTS_H_INCLUDED

#include "Header.h"

#define NUMBER_PARAMETERS 4 //number the parameters defined to be used when executing the program
#define ALIGNING_POSITIONS 50
#define NUMBER_COLOCAL_ALIGNMENTS 3

struct DP_cell
{
	int score;	
};


class AlignmentComponents
{
	public:
		AlignmentComponents(string input, int option, string parameters);
		~AlignmentComponents();
		int readFileParametersConfig();
		int readFileInputSequences();
		void printDataFromFiles();
		
		vector<string> sequences;
		int match;
		int mismatch;
		int h;
		int g;
		int option_alignment;

	private:		
		void generateOutputFile(string text);

		string input_file_name;		
		string parameters_file_name;
		
		

};
#endif