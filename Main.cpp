// ConsoleApplication3.cpp : Defines the entry point for the console application.
//

#include "AlignmentComponents.h"
#include "AlignmentFunctions.h"


//Function taht allows to check whether the option of alignment pass by the user is a 0 or a 1. No more digits can be used.
int readInt(string input)
{
	int result = 0;
	try
		{			
			result = stoi(input);
			if((result != 0) && (result != 1))
			{
				result = -1;
			}
		}
		catch(...)
		{		
			result = -1;
			return result;
		}
	return result;
}
 
/* ******************************************************************************************************************* */


int main(int argc, char * argv[])
{	
	try
	{				
		if((argc == NUMBER_PARAMETERS)||(argc == NUMBER_PARAMETERS-1)) //check the number of parameters that the user should have inserted
		{
			int option = readInt(argv[2]); //Validate that the inserted option 0 or 1
			if(option!=-1)	//if the option of alignment is correct
			{
				if(argc == NUMBER_PARAMETERS-1) argv[argc]="parameters.config"; //if the user didn't specify a file for the parameters, assign the default file
				vector<DP_cell**> matrices1;	//to store matrices S,D,I returned by the particular alignment
				vector<DP_cell**> matrices2;	//to store matrices S,D,I returned by local alignment when we are computing the second best one
				string previous1 = ""; string previous2 = ""; //used to store the portions of sequence that are before the ones align in local alignment
				string posterior1 = ""; string posterior2 = ""; //used to store the portions of sequence that are after the ones align in local alignment
				int max_score1 = 0; int max_score2 =0; //used to store the scores obtaines once we apply local alignment algorith to previous and posterior parts, if the parts exist.				
				int pass = 0; //to keep track if the files were read correctly, otherwise end the program
				AlignmentComponents aligment_components(argv[1],option,argv[3]); //create an instance of the class alignment components
				pass = aligment_components.readFileParametersConfig();//read parameters
				if(pass == 0)//correct
				{
					pass = aligment_components.readFileInputSequences();//read input
					if(pass==0)//correct
					{
						aligment_components.printDataFromFiles();				
						AlignmentFunctions aligment_functions(aligment_components.match,aligment_components.mismatch,aligment_components.h,aligment_components.g,aligment_components.option_alignment);
						string* sequences_align = aligment_functions.correctSequences(aligment_components.sequences[0],aligment_components.sequences[1]);				
						if((sequences_align[0].length()>0)&&(sequences_align[1].length()>0))//the sequences can't be empty
						{
							if(option==0)//global alignment
							{
								matrices1 = aligment_functions.globalAlignAffineGapPenalty(sequences_align[0],sequences_align[1]);
								aligment_functions.traceBackGlobalAlignment(sequences_align[0],sequences_align[1],matrices1[0],matrices1[1],matrices1[2]);						
							}
							else //local alignment
							{	
								//first best alignment
								matrices1 = aligment_functions.localAlignAffineGapPenalty(sequences_align[0],sequences_align[1]);						
								aligment_functions.traceBackLocalAlignment(sequences_align[0],sequences_align[1],matrices1[0],matrices1[1],matrices1[2],1);													
								//second best alignment
								if(aligment_functions.position_best_local_alignments[0]>1)//if not all the characters were consumed from the beggining of sequence1
									previous1 = sequences_align[0].substr(0,aligment_functions.position_best_local_alignments[0]-1);//store the subsequence
								if(aligment_functions.position_best_local_alignments[2]>1)//if not all the characters were consumed from the beggining of sequence2
									previous2 = sequences_align[1].substr(0,aligment_functions.position_best_local_alignments[2]-1);//store the subsequence
								if(aligment_functions.position_best_local_alignments[1]<sequences_align[0].length())//if not all the characters were consumed until the end of sequence1
									posterior1 = sequences_align[0].substr(aligment_functions.position_best_local_alignments[1],sequences_align[0].length());//store the subsequence
								if(aligment_functions.position_best_local_alignments[3]<sequences_align[1].length())//if not all the characters were consumed until the end of sequence1
									posterior2 = sequences_align[1].substr(aligment_functions.position_best_local_alignments[3],sequences_align[1].length());//store the subsequence
								if((previous1.length()>0)&&(previous2.length()>0))//if we have secuences to compare whose characters come from before
								{
									matrices1 = aligment_functions.localAlignAffineGapPenalty(previous1,previous2);
									int* values = aligment_functions.findMaximumMatrix(matrices1[0],previous1.length(),previous2.length());
									max_score1 = values[0];							
								}
								if((posterior1.length()>0)&&(posterior2.length()>0))//if we have secuences to compare whose characters come from after
								{
									matrices2 = aligment_functions.localAlignAffineGapPenalty(posterior1, posterior2);
									int* values = aligment_functions.findMaximumMatrix(matrices2[0],posterior1.length(),posterior2.length());
									max_score2 = values[0];							
								}
								if((max_score1!=0)||(max_score2!=0))//if there is oen or two second optimal scores
								{
									stringstream line;
									line << "::: Second Optimal Local Alignment :::" << endl  << string(50,'-') << endl;
									aligment_functions.generateOutputFile(line.str());
									if(max_score1>=max_score2)//comparison of which one has a higher value. If they are equal the traceback the one obtained from before substring
									{								
										aligment_functions.traceBackLocalAlignment(previous1,previous2,matrices1[0],matrices1[1],matrices1[2],2);	
									}
									else
									{								
										aligment_functions.traceBackLocalAlignment(posterior1,posterior2,matrices2[0],matrices2[1],matrices2[2],3);	
									}
								}
								else
								{
									cout << "::: There are no more disjoint optimal local alignments :::" << endl;
								}
							}									
						}
						else
						{
							cerr << "::: Incorrect input sequences to align. No empty string allowed :::" << endl;
							system("pause");		
							return -1;
						}
					}
					else
					{
						cerr << "::: Invalid input sequences file name :::" << endl;
						system("pause");		
						return -1;						
					}
				}
				else
				{
					cerr << "::: Invalid parameters file name :::" << endl;
					system("pause");		
					return -1;
				}
			}
			else
			{			
				cerr << "Usage: <executable name> <input file containing both s1 and s2> <0: global, 1: local> <optional: path to parameters config file>" << endl;				
				system("pause");		
				return -1;
			}
		}
		else
		{						
			cerr << "Usage: <executable name> <input file containing both s1 and s2> <0: global, 1: local> <optional: path to parameters config file>" << endl;				
			system("pause");		
			return -1;
		}
	}
	catch (std::runtime_error &e)
	{
		std::cout<<e.what()<<std::endl;
	}
	catch(...)
	{
		std::cout<<"Unknown error Main"<<std::endl;
	}
    
		cout << "Alignment Algorithm has finished."<< endl;
		system("pause");
		return 0;	
}	
	