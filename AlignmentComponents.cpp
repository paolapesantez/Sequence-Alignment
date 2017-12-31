#include "AlignmentComponents.h"
#include "AlignmentFunctions.h"
#include "StringSplitter.h"


//Function that deletes space from a whole string
string trim(const string &s)
{
  int last = s.size() - 1;
  while (last >= 0 && s[last] == ' ')
    --last;
  return s.substr(0, last + 1);
}


/* ****************************************************************************************************************** */

//Class constructor that receives as parameters the information input by the user while trying to execute the program
AlignmentComponents::AlignmentComponents(string input, int option, string parameters)
{
	input_file_name = input;				//name of the file that contains the input sequences to align
	option_alignment =  option;				//specifies if the wanted alignment is 0:global or 1:local
	parameters_file_name = parameters;		//name of the file that contains the requiered parameters to perform the alignment 
	match=0;								//award value for every match
	mismatch=0;								//penalty value for each mismatch
	h=0;									//penalty for an opening gap
	g=0;									//penalty for a continuation gap
}


/* ****************************************************************************************************************** */

//Class destructot
AlignmentComponents::~AlignmentComponents(){}


/* ****************************************************************************************************************** */

//Procedure that opens and read the file that contains the requiered parameters to perform the alignment 
int AlignmentComponents::readFileParametersConfig()
{		
	
	ifstream input_file(parameters_file_name);			
	bool is_good = true;
	try
	{
		if(input_file.is_open() == false)				//If the file hasn't been found
		{		
			return -1;									//no file has been found. End the program.
			is_good = false;
		}	
		if(is_good)									//If the file has been found
		{
			while (input_file.good())				//read line by line
			{
				int items_found = 0;
				string str_line = "";
				getline(input_file,str_line);
				string *pieces = StringSplitter::split(str_line," ",items_found); //StringSplitter is a class that parses a line according a delimiter, dividing the line in its different parts.
				if (pieces[0].compare("match") == 0)	//store the parameters
				{
					match = stoi(pieces[1]);
				}
				else if (pieces[0].compare("mismatch") == 0)
				{
					mismatch = stoi(pieces[1]);
				}
				else if (pieces[0].compare("h") == 0)
				{
					h = stoi(pieces[1]);
				}
				else if (pieces[0].compare("g") == 0)
				{
					g = stoi(pieces[1]);
				}				
				delete[] pieces;		// pieces is a pointer so we need to deallocate the memory being occupied
			}			
			input_file.close();			//close the file from which we were reading
			return 0;					//reading completed succesfully
		}
	}
	catch (runtime_error &e)
	{
		cout<<e.what()<< " read file parameters " <<endl;
	}
	catch(...)
	{
		cout<<"Unknown error in readFileParamentes"<<endl;
	}
}

/* ****************************************************************************************************************** */

//Procedure that opens and read the file that contains the input sequences over which we will perform the alignment 
int AlignmentComponents::readFileInputSequences()
{	
	ifstream input_file(input_file_name);
	bool is_good = true;
	try
	{
		if(input_file.is_open() == false)		//If the file hasn't been found
		{			
			return -1;							//no file has been found. End the program.
			is_good = false;
		}	
		if(is_good)								//If the file has been found
		{
			stringstream sequence;				//used to store the complete sequence that in the file occupies several lines			
			while (input_file.good())			//read line by line
			{
				int items_found = 0;				
				string str_line = "";
				getline(input_file,str_line);
				string *pieces = StringSplitter::split(str_line," ",items_found); //StringSplitter is a class that parses a line according a delimiter, dividing the line in its different parts.
				if ((items_found > 0) && (pieces[0].length()>0)) // When we find an empty line we know the next one will have the start of the next sequence
				{
					if (pieces[0][0] == '>')	//start of a sequence. Parsing the name of the sequence
					{
						sequence.str("");
						if(isalpha(pieces[0][1]))
						{
							sequence << pieces[0].substr(1,pieces[0].length()) << "/";//extract the name of the sequence and append / to it so we can separate the name from the sequence itself later
						}
						else
						{
							sequence << "s" << sequences.size()+1 << "/"; //in case the name of the sequence doesn't have the default format we will give it the right name
						}								
					}	
					else //if we are actually reading the lines that represent the sequence content itself
					{
						sequence << trim(str_line);	//append the different lines in just one
					}
			   }
				else
				{					
					sequences.push_back(sequence.str()); //when the sequence has finished put it in the vector of sequences
				}
			}
			sequences.push_back(sequence.str()); //we have finished reading the file so the last sequence also needs to be stored
			input_file.close();			//close the file from which we were reading
			return 0;					//reading completed succesfully
		}
	}
	catch (runtime_error &e)
	{
		cout<<e.what()<< " read file sequences " <<endl;
	}
	catch(...)
	{
		cout<<"Unknown error in readFileInputSequences"<<endl;
	}
}

/* ****************************************************************************************************************** */

//Printing to the screen the data that has been read for the files
void AlignmentComponents::printDataFromFiles()
{
	try
	{		
		string filename;
		stringstream line;
		cout << endl << endl;
		if (option_alignment == 0)
		{
			line.str("");
			line << "::: Global Alignment with Affine Gap Penalty :::" << endl << string(50,'-') << endl << endl;			
			cout << line.str();
			filename = "globalAlignment.txt";
			ofstream outfile(filename);  
			outfile << line.str();
			outfile.close();
		}
		else
		{	
			line.str("");
			line << "::: Local Alignment with Affine Gap Penalty :::" << endl << string(50,'-') << endl << endl;	
			cout << line.str();
			filename = "localAlignment.txt";
			ofstream outfile(filename);  
			outfile << line.str();
			outfile.close();	
		}
				
		line.str("");
		line << "Scores:    match =" << match << ", mismatch =" << mismatch << ", h =" << h << ", g =" << g << endl << endl;
		generateOutputFile(line.str());	
		for (int i=0; i<sequences.size();i++)
		{
			int pos = sequences[i].find_first_of("/");
			int sequence_length = sequences[i].length()-(pos+1);
			line.str("");
			line << "Sequence " << i+1 << ": \"" << sequences[i].substr(0,pos) << "\", length = " << sequence_length << endl;			
			generateOutputFile(line.str());				
		}	
		line.str("");
		line << endl;
		generateOutputFile(line.str());				
	}
	catch (runtime_error &e)
	{
		cout<<e.what()<< " printFile " <<endl;
	}
	catch(...)
	{
		cout<<"Unknown error in printFile"<<endl;
	}
}


/* ****************************************************************************************************************** */

//Procedure that receives a string and append the string to the output file that has been created for the particular alignment.
//Also it prints the string to the screen
void AlignmentComponents::generateOutputFile(string text)
{
	string filename;

	if(option_alignment == 0)
		filename="globalAlignment.txt";
	else
		filename="localAlignment.txt";
	
	ofstream outfile;  
	outfile.open(filename, ios::app );
	outfile << text;
	outfile.close();
	cout << text;	
}



