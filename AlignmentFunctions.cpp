#include "AlignmentFunctions.h"


//Class constructor that receives as input the parameters for alignment computation
AlignmentFunctions::AlignmentFunctions(int ma, int mi, int h, int g, int option)
{	
	match = ma;
	mismatch = mi;
	h_gap = h;
	g_gap = g;	
	option_alignment = option;
	maximum_score = 0;
	name_sequence1="";
	name_sequence2="";
}


/* ****************************************************************************************************************** */

//Class destructor
AlignmentFunctions::~AlignmentFunctions()
{
}



/* ****************************************************************************************************************** */

//Function that divides the sequence in two parts: name and content
//This is because when I read a sequence from the file I store it like s1/AAAAGTGGGT. So to perform the alignment we jsut need the content
string* AlignmentFunctions::correctSequences(string sequence1, string sequence2)
{
		string *result = new string[2];
		int pos = sequence1.find_first_of("/");
		name_sequence1 = sequence1.substr(0,pos);
		result[0] = sequence1.substr(pos+1,sequence1.length());
		pos = sequence2.find_first_of("/");
		name_sequence2 = sequence2.substr(0,pos);
		result[1] = sequence2.substr(pos+1,sequence2.length());		
		return result;
}


/* ****************************************************************************************************************** */

//Function that initialize each matrix S, D, or I for the global alignment according to option_matrix and according the concepts learned in class
 DP_cell** AlignmentFunctions::initializeMatricesGlobalAlignment(int m, int n,int option_matrix)
{
	//initialize the matrix
	DP_cell** matrix = 0;					//A pointer to pointers to a DP cell.	
	matrix = new DP_cell*[m];				//Matrix is now a pointer to an array of 'rows' pointers.
	
	//define the matrix
	for(int i=0;i<m;i++)
	{
		matrix[i] = new DP_cell[n];		    //the ith array is initialized
		for(int j=0;j<n;j++)				//the i,jth element is defined
		{			
			matrix[i][j].score = 0;					
		}
	}
		

	if((option_matrix==2)||(option_matrix==3))//First cell should be inf in D and I
	{
		matrix[0][0].score = -1000000000;			
	}	

	for(int i=1;i<m;i++)
	{
		if((option_matrix==1)||(option_matrix==3))//Cells in first column, all rows should be inf in S and I
		{
			matrix[i][0].score = -1000000000;			
		}	
		else
		{
			matrix[i][0].score = h_gap+i*g_gap;	//Gap Penalty in D
		}				
		
	}

	for(int j=1;j<n;j++)
	{
		if((option_matrix==1)||(option_matrix==2))//Cells in first row, all columns should be inf in S and D
		{
			matrix[0][j].score = -1000000000;			
		}	
		else
		{
			matrix[0][j].score = h_gap+j*g_gap;	//Gap Penalty in I		
		}			
	}	
	return matrix;
}

/* ****************************************************************************************************************** */

 //Function that initialize each matrix S, D, or I for the local alignment according to option_matrix and according the concepts learned in class
 DP_cell** AlignmentFunctions::initializeMatricesLocalAlignment(int m, int n,int option_matrix)
{

	//initialize the matrix
	DP_cell** matrix = 0;					//A pointer to pointers to a DP cell.
	matrix = new DP_cell*[m];				//Matrix is now a pointer to an array of 'rows' pointers.
	
	//define the matrix
	for(int i=0;i<m;i++)
	{
		matrix[i] = new DP_cell[n];		    //the ith array is initialized
		for(int j=0;j<n;j++)				//the i,jth element is defined
		{			
			matrix[i][j].score = 0;			
		}
	}	

	for(int i=0;i<m;i++)
	{
		if((option_matrix==2)||(option_matrix==3))//Cells in first column, all rows should be inf in D and I
		{
			matrix[i][0].score = -100000000;			
		}			
	}

	for(int j=0;j<n;j++)
	{
		if((option_matrix==2)||(option_matrix==3))//Cells in first row, all columns should be inf in D and I
		{
			matrix[0][j].score = -100000000;			
		}				
	}	
	return matrix;
}


/* ****************************************************************************************************************** */

//Function that returns the cost of a substitution operation that can be a match or a mismatch
int AlignmentFunctions::substitutionEdit(char a, char b)
{
	if(a == b)
		return match;
	else
		return mismatch;
}

/* ****************************************************************************************************************** */

//Function that returns true if a string is contained in another string
bool AlignmentFunctions::findSubstring(string s1, string s2)
{
	bool result = false;
	if (s1.find(s2) != std::string::npos) 
	{
		result = true;
	}
	return result;
}


/* ****************************************************************************************************************** */

//Function that returns the maximum value from the entries of a given vector. Number of item tells how many numbers inside the vector I want to compare
int AlignmentFunctions::findMaximumVector(int* previous_scores, int number_items)
{
	int result;
	result = previous_scores[0];

	for(int i=1;i<number_items;i++) //compare data with max
	{
		result = max(previous_scores[i], result);
	}	
	return result;
}


/* ****************************************************************************************************************** */

//Function that returns the maximum score stored in a matrix. Also it returns the position(i,j) where the maximum value was found
int* AlignmentFunctions::findMaximumMatrix(DP_cell** matrix, int m, int n)
{	
	int* result = new int[3];
	result[0] = 0;
	result[1] = 0;
	result[2] = 0;
	for(int i=0; i<=m; i++)
	{
		for(int j=0; j<=n; j++)
		{
			if(matrix[i][j].score > result[0])
			{
				result[0] = matrix[i][j].score;
				result[1] = i;
				result[2] = j;
			}
		}
	}
	return result;
}


/* ****************************************************************************************************************** */

//Procedure that prints the content of the specified matrix to the screen and to a file
void AlignmentFunctions::printMatrix(DP_cell** matrix, int m, int n)
{
	string filename = "matrix.txt";
	ofstream outfile(filename);  
	
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			cout << matrix[i][j].score << '\t';
			outfile << matrix[i][j].score << '\t';
			
		}
		cout << "Row" << i << endl;
		outfile << "Row" << i << endl;
	}
	outfile.close();	
}

/* ****************************************************************************************************************** */

//Function that copy the content of a specified matrix to another
DP_cell** AlignmentFunctions::copyMatrix(DP_cell** matrix, int m, int n)
{
	DP_cell** result = 0;
	result = new DP_cell*[m];
	for(int i=0;i<m;i++)
	{
		result[i] = new DP_cell[n];	
		for(int j=0;j<n;j++)
		{
			result[i][j].score = matrix[i][j].score;			
		}		
	}	
	return result;
}

/* ****************************************************************************************************************** */

//Function that returns Dynamic Programming Tables S,D,and I after the global dynamic programing recurrences are executed
vector<DP_cell**> AlignmentFunctions::globalAlignAffineGapPenalty(string sequence1, string sequence2)
{	
	int previous_scores[3];			//to store the scores of a particular cell in S, D, and I
	vector<DP_cell**>result_matrices;
	int result;
	
	int m = sequence1.length()+1;
	int n = sequence2.length()+1;

	//initialization of the matrices
	DP_cell** S = initializeMatricesGlobalAlignment(m,n,1);			
	DP_cell** D = initializeMatricesGlobalAlignment(m,n,2);	
	DP_cell** I = initializeMatricesGlobalAlignment(m,n,3);		
	
	//execution of the recurrences to store in each cell of each matrix the maximum score
	for(int i=1;i<m;i++)
	{
		for(int j=1;j<n;j++)
		{
			previous_scores[0] = S[i-1][j-1].score;
			previous_scores[1] = D[i-1][j-1].score;
			previous_scores[2] = I[i-1][j-1].score;						
			result = findMaximumVector(previous_scores,3);	//selecting the maximum score
			S[i][j].score = substitutionEdit(sequence1[i-1],sequence2[j-1])+result;		//Filling matrix S	
			
			previous_scores[0] = S[i-1][j].score+h_gap+g_gap;				
			previous_scores[1] = D[i-1][j].score+g_gap;												
			result = findMaximumVector(previous_scores,2); //selecting the maximum score
			D[i][j].score = result; //Filling matrix D

			previous_scores[0] = S[i][j-1].score+h_gap+g_gap;			
			previous_scores[1] = I[i][j-1].score+g_gap;			
			result = findMaximumVector(previous_scores,2); //selecting the maximum score
			I[i][j].score = result;			//Filling matrix I
		}
	}
			
	/*printMatrix(S,m,n);
	printMatrix(D,m,n);
	printMatrix(I,m,n);*/
	result_matrices.push_back(S);
	result_matrices.push_back(D);
	result_matrices.push_back(I);	
	return result_matrices;	
}

/* ****************************************************************************************************************** */

//Function that returns Dynamic Programming Tables S,D,and I after the local dynamic programing recurrences are executed
vector<DP_cell**> AlignmentFunctions::localAlignAffineGapPenalty(string sequence1, string sequence2)
{
	int previous_scores[3];		//to store the scores of a particular cell in S, D, and I
	int result;					
	vector<DP_cell**>result_matrices;
	
	int m = sequence1.length()+1;
	int n = sequence2.length()+1;

	//initialization of the matrices
	DP_cell** S = initializeMatricesLocalAlignment(m,n,1);			
	DP_cell** D = initializeMatricesLocalAlignment(m,n,2);	
	DP_cell** I = initializeMatricesLocalAlignment(m,n,3);		

	//execution of the recurrences to store in each cell of each matrix the maximum score
	for(int i=1;i<m;i++)
	{
		for(int j=1;j<n;j++)
		{
			previous_scores[0] = S[i-1][j-1].score;
			previous_scores[1] = D[i-1][j-1].score;
			previous_scores[2] = I[i-1][j-1].score;						
			result = findMaximumVector(previous_scores,3);		//selecting the maximum score
			S[i][j].score = substitutionEdit(sequence1[i-1],sequence2[j-1])+result;		//Filling matrix S	
			if(S[i][j].score < 0) 
			{
					S[i][j].score = 0;					//If the score is negative it should be stored in S as 0
			}
			
			previous_scores[0] = S[i-1][j].score+h_gap+g_gap;				
			previous_scores[1] = D[i-1][j].score+g_gap;				
			result = findMaximumVector(previous_scores,2);//selecting the maximum score
			D[i][j].score = result; //Filling matrix D

			previous_scores[0] = S[i][j-1].score+h_gap+g_gap;			
			previous_scores[1] = I[i][j-1].score+g_gap;	
			result = findMaximumVector(previous_scores,2);//selecting the maximum score
			I[i][j].score = result;		//Filling matrix I	
		}
	}
	/*printMatrix(S,m,n);
	printMatrix(D,m,n);
	printMatrix(I,m,n);*/
	result_matrices.push_back(S);
	result_matrices.push_back(D);
	result_matrices.push_back(I);	
	return result_matrices;
}


/* ****************************************************************************************************************** */

//Procedure that performs the traceback for the global alignment
void AlignmentFunctions::traceBackGlobalAlignment(string sequence1, string sequence2, DP_cell** S, DP_cell** D, DP_cell** I)
{		
	int matrix_best_score;	//best score of the global alignment
	string current_matrix;	//matrix that has the best score	
	stringstream s1,s2,sm;	//to store the alignment bases according they have been consumed in the dynamic programing tables	
	int matches=0, mismatches=0,gaps=0,open_gaps=0; //to report at the end the number of matches, mismatches, openin gaps and gaps in general
	
	int m = sequence1.length();
	int n = sequence2.length();

	//Select the maximum calculated score among the ones obtained in S,D and I and set current_matrix as the matrix that has the maximum score
	maximum_score = S[m][n].score;
	current_matrix = "S";
	if(maximum_score < D[m][n].score)
	{
		maximum_score = D[m][n].score;
		current_matrix = "D";
	}
	else if(maximum_score < I[m][n].score)
	{
		maximum_score = I[m][n].score;
		current_matrix = "I";
	}	

	bool band =true;	
	s1.str("");
	s2.str("");
	sm.str("");
	while(band == true)
	{				
		if(current_matrix.compare("S")==0)//substitution operation was performed
		{
			int penalty =substitutionEdit(sequence1[m-1],sequence2[n-1]);
			s1 << sequence1[m-1] << "";
			if(penalty>0) //it means a match is a result of the substitution operation
			{	sm << "|" << "";	
				matches++;
			}
			else //it means a mismatch is a result of the substitution operation
			{	sm << " " << "";	
				mismatches++;
			}
			s2 << sequence2[n-1] << "";
			m--;
			n--;
			if((S[m][n].score+penalty)==S[m+1][n+1].score)
				current_matrix = "S";
			else if((D[m][n].score+penalty)==S[m+1][n+1].score)
				current_matrix = "D";
			else if((I[m][n].score+penalty)==S[m+1][n+1].score)
				current_matrix = "I";
		}
		else if(current_matrix.compare("D")==0)//deletion operation was performed
		{
			s1 << sequence1[m-1] << "";
			sm << " " << "";
			s2 << "_" << "";
			m--;
			if((S[m][n].score+h_gap+g_gap)==D[m+1][n].score)
			{	current_matrix = "S";
				open_gaps++;
			}
			else if((D[m][n].score+g_gap)==D[m+1][n].score)
			{	current_matrix = "D";
				gaps++;
			}		
			else if((I[m][n].score+h_gap+g_gap)==D[m+1][n].score)
			{	current_matrix = "I";
				open_gaps++;
			}
		}
		else if(current_matrix.compare("I")==0)//insertion operation was performed
		{
			s1 << "_" << "";
			sm << " " << "";
			s2 << sequence2[n-1] << "";					
			n--;
			if((S[m][n].score+h_gap+g_gap)==I[m][n+1].score)
			{	current_matrix = "S";
				open_gaps++;
			}
			else if((D[m][n].score+h_gap+g_gap)==I[m][n+1].score)
			{	current_matrix = "D";
				open_gaps++;
			}
			else if((I[m][n].score+g_gap)==I[m][n+1].score)
			{	current_matrix = "I";
				gaps++;
			}
		}

		if(m==0 && n==0)//if we have reached the origin, then stop the retrace process
		{
			band = false;
		}		
	}
	printOutputGlobalAlignment(s1.str(),sm.str(),s2.str());	//print the sequence alignment output for the global alignment
	stringstream line;
	line << "::: Report :::" << endl << string(20,'-') << endl;	//print the report
	generateOutputFile(line.str());
	line.str("");
	line << "Maximum Score: " << maximum_score << endl;
	generateOutputFile(line.str());
	line.str("");
	line << "Matches: " << matches << endl;
	generateOutputFile(line.str());
	line.str("");
	line << "Mismatches: " << mismatches << endl;
	generateOutputFile(line.str());
	line.str("");
	line << "Openin Gaps: " << open_gaps << endl;
	generateOutputFile(line.str());
	line.str("");
	line << "Gaps: " << gaps+open_gaps << endl;
	generateOutputFile(line.str());
	float percent = (float(matches) / sequence1.length())*100;
	line.str("");
	line << "Identities= " << percent << "%" << endl;
	generateOutputFile(line.str());
	percent = (double(gaps+open_gaps) / sequence1.length())*100;
	line.str("");
	line << "Gaps= " << percent  << "%" << endl  << endl;
	generateOutputFile(line.str());

	for(int i=0;i<sequence1.length();i++)//deallocate the space in memory that the dynamic programming tables are consuming
	{
		delete[] S[i];	
		delete[] D[i];	
		delete[] I[i];	
	}
	delete [] S;
	delete [] D;
	delete [] I;
}


/* ****************************************************************************************************************** */


//Procedure that performs the traceback for the local alignment
void AlignmentFunctions::traceBackLocalAlignment(string sequence1, string sequence2, DP_cell** S, DP_cell** D, DP_cell** I, int position_option)
{			
	int matrix_best_score;	//best score of the global alignment
	string current_matrix;	//matrix that has the best score	
	stringstream s1,s2,sm;	//to store the alignment bases according they have been consumed in the dynamic programing tables	
	int matches=0, mismatches=0,gaps=0,open_gaps=0; //to report at the end the number of matches, mismatches, openin gaps and gaps in general	
	int pos_i=0, pos_j=0;	//position(i,j) where the maximum vale was found inside matrix S

	int m = sequence1.length();
	int n = sequence2.length();

	int* values = findMaximumMatrix(S,m,n);
	maximum_score = values[0];
	if(maximum_score != 0)
	{//cout << "MS: " << maximum_score << endl;		
		pos_i = values[1];
		pos_j = values[2];
		current_matrix = "S";
	
		bool band =true;
		s1.str("");
		s2.str("");
		sm.str("");
		while(band == true)
		{			
			if(current_matrix.compare("S")==0)//substitution operation was performed
			{
				if(S[pos_i][pos_j].score!=0)
				{
					int penalty =substitutionEdit(sequence1[pos_i-1],sequence2[pos_j-1]);
					s1 << sequence1[pos_i-1] << "";
					if(penalty>0)//it means a match is a result of the substitution operation
					{	sm << "|" << "";	
						matches++;
					}
					else//it means a mismatch is a result of the substitution operation
					{	sm << " " << "";	
						mismatches++;
					}
					s2 << sequence2[pos_j-1] << "";										
					pos_i--;
					pos_j--;
									
					if((S[pos_i][pos_j].score+penalty)==S[pos_i+1][pos_j+1].score)
						current_matrix = "S";
					else if((D[pos_i][pos_j].score+penalty)==S[pos_i+1][pos_j+1].score)
						current_matrix = "D";
					else if((I[pos_i][pos_j].score+penalty)==S[pos_i+1][pos_j+1].score)
						current_matrix = "I";				
				}				
			}
			else if(current_matrix.compare("D")==0)//deletion operation was performed
			{
				s1 << sequence1[pos_i-1] << "";
				sm << " " << "";
				s2 << "_" << "";								
				pos_i--;
				if((S[pos_i][pos_j].score+h_gap+g_gap)==D[pos_i+1][pos_j].score)
				{	current_matrix = "S";
					open_gaps++;
				}
				else if((D[pos_i][pos_j].score+g_gap)==D[pos_i+1][pos_j].score)
				{	current_matrix = "D";
					gaps++;
				}		
				else if((I[pos_i][pos_j].score+h_gap+g_gap)==D[pos_i+1][pos_j].score)
				{	current_matrix = "I";
					open_gaps++;
				}
			}
			else if(current_matrix.compare("I")==0)//insertion operation was performed
			{
				s1 << "_" << "";
				sm << " " << "";
				s2 << sequence2[pos_j-1] << "";								
				pos_j--;
				if((S[pos_i][pos_j].score+h_gap+g_gap)==I[pos_i][pos_j+1].score)
				{	current_matrix = "S";
					open_gaps++;
				}
				else if((D[pos_i][pos_j].score+h_gap+g_gap)==I[pos_i][pos_j+1].score)
				{	current_matrix = "D";
					open_gaps++;
				}
				else if((I[pos_i][pos_j].score+g_gap)==I[pos_i][pos_j+1].score)
				{	current_matrix = "I";
					gaps++;
				}
			}	
		
			if(S[pos_i][pos_j].score==0)//if we have reached a score of 0 in the S matrix we should stop the retrace
			{	
				band = false;	
			}				
		}	

		printOutputLocalAlignment(s1.str(),sm.str(),s2.str(),pos_i+1,pos_j+1, position_option);	//print the sequence alignment output for the local alignment
		stringstream line;
		line << "::: Report :::" << endl << string(20,'-') << endl;	//print the report
		generateOutputFile(line.str());
		line.str("");
		line << "Maximum Score: " << maximum_score << endl;
		generateOutputFile(line.str());
		line.str("");
		line << "Matches: " << matches << endl;
		generateOutputFile(line.str());
		line.str("");
		line << "Mismatches: " << mismatches << endl;
		generateOutputFile(line.str());
		line.str("");
		line << "Openin Gaps: " << open_gaps << endl;
		generateOutputFile(line.str());
		line.str("");
		line << "Gaps: " << gaps+open_gaps << endl;
		generateOutputFile(line.str());
		float percent = (float(matches) / sequence1.length())*100;
		line.str("");
		line << "Identities= " << percent << "%" << endl;
		generateOutputFile(line.str());
		percent = (double(gaps+open_gaps) / sequence1.length())*100;
		line.str("");
		line << "Gaps= " << percent  << "%" << endl  << endl << endl << endl;
		generateOutputFile(line.str());		

		for(int i=0;i<m;i++)//deallocate the space in memory that the dynamic programming tables are consuming
		{
			delete[] S[i];	
			delete[] D[i];	
			delete[] I[i];	
		}
		delete [] S;
		delete [] D;
		delete [] I;
	}
	else
	{
		cout << "There is no score greater than 0." << endl;
	}
}


/* *************************************************************************************************************** */

//Once we have the sequence obtained in the retrace this procedure prints out the sequence according to a numbedr of positions to show to the screen
//and in the accurate format
void AlignmentFunctions::printOutputGlobalAlignment(string s1,string sm,string s2)
{
	int cont = 0;					//to count how many bases have been printed already
	int consum_s1=0,consum_s2=0;	//to keep track of characters consumed so we print out this value
	bool band = false;
	stringstream s1_start,s2_start,sm_start; //to define the start format of the sequence
	
	s1_start.str("");
	s2_start.str("");
	sm_start.str("");	
	s1_start << name_sequence1 << ":" << 1 << "\t";
	sm_start << "\t";
	s2_start << name_sequence2 << ":" << 1 << "\t";
	reverse(s1.begin(), s1.end()); // we need to reverse since bactrace gives us the alignment sequence in reverse order
	reverse(sm.begin(), sm.end());
	reverse(s2.begin(), s2.end());

	cout << endl;

	for(int i=0;i<s1.length();i++)
	{
		s1_start << s1[i];
		if(s1[i]!='_') consum_s1 = consum_s1+1; //count the number of charracters that have been consumed from Sequence1
		sm_start << sm[i];
		s2_start << s2[i] ;
		if(s2[i]!='_') consum_s2 = consum_s2+1; //count the number of charracters that have been consumed from Sequence2
		cont++;
	
		if(((cont % ALIGNING_POSITIONS) == 0) || (i==s1.length()-1))
		{					
			s1_start  << " " << consum_s1 << endl;
			sm_start  << " " << endl;
			s2_start  << " " << consum_s2 << endl << endl;
			
			generateOutputFile(s1_start.str());
			generateOutputFile(sm_start.str());
			generateOutputFile(s2_start.str());
			s1_start.str("");
			s2_start.str("");
			sm_start.str("");
			s1_start << name_sequence1 << ":" << consum_s1+1 << "\t";
			sm_start << "\t";
			s2_start << name_sequence2 << ":" << consum_s2+1 << "\t";
			band = true;

			if ((consum_s1 >= 10000 )&&(consum_s2 < 10000 )) //to keep format accurately
			{				
				sm_start << '\t';
				s2_start << '\t';
			}
			else if  ((consum_s1 >= 10000 )&&(consum_s2 >= 10000 ))
			{
				sm_start << '\t';
			}
		}		
	}
	sm_start.str("");
	sm_start << endl;
	generateOutputFile(sm_start.str());	//here we add the desired portion of the global oputput sequence to the oputput file and to the screen			
}



/* *************************************************************************************************************** */

//Once we have the sequence obtained in the retrace this procedure prints out the sequence according to a numbedr of positions to show to the screen
//and in the accurate format
void AlignmentFunctions::printOutputLocalAlignment(string s1,string sm,string s2, int initial_pos_s1, int initial_pos_s2, int position_option)
{
	int cont = 0;							//to count how many bases have been printed already
	int consum_s1=0,consum_s2=0;			//to keep track of characters consumed so we print out this value
	bool band = false;
	stringstream s1_start,s2_start,sm_start; //to define the start format of the sequence
	int initial1=0, initial2=0;	//since the alignment is local we need to keep trrack of the correct initial positions in sequence 1 and 2
	s1_start.str("");
	s2_start.str("");
	sm_start.str("");	
	if(position_option==1) //it means we are printing the local alignment for the first best score
	{
		initial1 = initial_pos_s1;
		initial2 = initial_pos_s2;
		s1_start << name_sequence1 << ":" << initial1 << "\t";		
		s2_start << name_sequence2 << ":" << initial2 << "\t";
	}
	if(position_option==2)//it means we are printing the local alignment for the second best score obtained from the subsequence previous to the first best
	{
		initial1 = position_best_local_alignments[0]-initial_pos_s1;
		initial2 = position_best_local_alignments[2]-initial_pos_s2;
		s1_start << name_sequence1 << ":" << initial1 << "\t";		
		s2_start << name_sequence2 << ":" << initial2 << "\t";
	}
	if(position_option==3)//it means we are printing the local alignment for the second best score obtained from the subsequence posterior to the first best
	{
		initial1 = initial_pos_s1+position_best_local_alignments[1];
		initial2 = initial_pos_s2+position_best_local_alignments[3];
		s1_start << name_sequence1 << ":" << initial1 << "\t";		
		s2_start << name_sequence2 << ":" << initial2 << "\t";
	}
	sm_start << "\t";
	if((initial1 >= 10000 )&&(initial2 < 10000 ))
	{				
		sm_start << '\t';
		s2_start << '\t';
	}
	else if((initial1 >= 10000 )&&(initial2 >= 10000 ))
	{
		sm_start << '\t';
	}
	reverse(s1.begin(), s1.end());// we need to reverse since bactrace gives us the alignment sequence in reverse order
	reverse(sm.begin(), sm.end());
	reverse(s2.begin(), s2.end());

	cout << endl;

	for(int i=0;i<s1.length();i++)
	{
		s1_start << s1[i];
		if(s1[i]!='_') consum_s1 = consum_s1+1;
		sm_start << sm[i];
		s2_start << s2[i] ;
		if(s2[i]!='_') consum_s2 = consum_s2+1;
		cont++;
	
		if(((cont % ALIGNING_POSITIONS) == 0) || (i==s1.length()-1))
		{					
			if(band == true)
			{	s1_start  << " " << consum_s1 << endl;
				sm_start  << " " << endl;
				s2_start  << " " << consum_s2 << endl << endl;
			}
			else 
			{//the final position of the local alignment needs to be calculated accurately
				if(position_option == 1) //it means we are printing the local alignment for the first best score
				{
					consum_s1 = (consum_s1+initial1)-1;
					consum_s2 = (consum_s2+initial2)-1;
				}
				if(position_option == 2)//it means we are printing the local alignment for the second best score obtained from the subsequence previous to the first best
				{
					consum_s1 = (consum_s1+initial1)-1;
					consum_s2 = (consum_s2+initial2)-1;
				}
				if(position_option == 3)//it means we are printing the local alignment for the second best score obtained from the subsequence posterior to the first best
				{
					consum_s1 = (consum_s1+initial1)-1;
					consum_s2 = (consum_s2+initial2)-1;
				}
				s1_start  << " " << consum_s1 << endl;
				sm_start  << " " << endl;
				s2_start  << " " << consum_s2 << endl << endl;				
			}
			generateOutputFile(s1_start.str());
			generateOutputFile(sm_start.str());
			generateOutputFile(s2_start.str());
			s1_start.str("");
			s2_start.str("");
			sm_start.str("");
			s1_start << name_sequence1 << ":" << consum_s1+1 << "\t";
			sm_start << "\t";
			s2_start << name_sequence2 << ":" << consum_s2+1 << "\t";
			band = true;

			if ((consum_s1 >= 10000 )&&(consum_s2 < 10000 ))
			{				
				sm_start << '\t';
				s2_start << '\t';
			}
			else if  ((consum_s1 >= 10000 )&&(consum_s2 >= 10000 ))
			{
				sm_start << '\t';
			}
		}		
	}
	sm_start.str("");
	sm_start << endl;
	generateOutputFile(sm_start.str());	//here we add the desired portion of the global oputput sequence to the oputput file and to the screen			

	//At the end of the local alignment we stores the initial and final positions of sequence1 and sequence2 so next time we can remove the local sequence alignment and start a new alignment of the previous of posterior parts of it
	position_best_local_alignments[0] = initial_pos_s1;
	position_best_local_alignments[1] = consum_s1;
	position_best_local_alignments[2] = initial_pos_s2;
	position_best_local_alignments[3] = consum_s2;

}

/* *************************************************************************************************************** */

//Procedure that receives a string and append the string to the output file that has been created for the particular alignment.
//Also it prints the string to the screen
void AlignmentFunctions::generateOutputFile(string text)
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






