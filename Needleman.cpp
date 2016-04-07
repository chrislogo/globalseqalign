#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>

using namespace std;

// consts
const int ROW_LEN = 24;			// size of blosum62 matrix rows
const int GAP_PENALTY = -4;

// structures

// holds the BLOSUM62 matrix
struct Scores_Matrix
{
	char x_aminos[ROW_LEN];				// x-axis aminos
	char y_aminos[ROW_LEN];				// y-axis aminos
	int scores[ROW_LEN][ROW_LEN];		// scores for aminos
};

// the computed matrix
struct Possibilities
{
	int **results;				// for finding highest value
	char **movements;			// for backtracing final result
	int cols;
	int rows;
	string original1;
	string original2;
	string optimal1;
	string optimal2;

	bool file1_first;	// for correct output
};


// function prototypes
void Load_BLOSUM(Scores_Matrix&);
string Extract_Seq(string);
Possibilities Compute_Poss(Scores_Matrix, string, string);
int Find_Pos(Scores_Matrix, char);
int Find_Greatest_Val(int, int, int, int, char&);
void Traceback(Possibilities&);
string Reverse_Str(string);


int main(int argc, char **argv)
{
	Scores_Matrix blosum62;
	Possibilities P;
	string file1;
	string file2;
	string sequence1;
	string sequence2;

	if(argc < 3)
	{
		cout << "Error! Usage: Needleman <file_name> <file_name>" << endl;
		exit(1);
	}
	else
	{
		file1 = argv[1];
		file2 = argv[2];
	}

	// load matrix
	Load_BLOSUM(blosum62);

	

	// extract the 2 sequences
	sequence1 = Extract_Seq(file1);
	sequence2 = Extract_Seq(file2);

	// Compute traceback table
	P = Compute_Poss(blosum62, sequence1, sequence2);

	// follow path to find best alignment
	Traceback(P);

	// x-axis is always longest so if sequence2 was longer, numbers were flipped
	// print in reverse to preserve output of file1 alignment then file2 alignment
	if(P.file1_first)
	{
		cout << P.optimal1 << endl
			 << P.optimal2 << endl; 
	}
	else
	{
		cout << P.optimal2 << endl
			 << P.optimal1 << endl;
	}

	// clean up
	for(int i = 0; i < P.rows; i++)
	{
    	delete[] P.movements[i];
    	delete[] P.results[i];
	}

	delete[] P.movements;
	delete[] P.results; 

	return 0;
}


// load matrix into 2D array
// I don't include the * -4 row and column since the gap penalty is hardcoded
void Load_BLOSUM(Scores_Matrix &blosum62)
{
	int j = 0;
	int i = 3;
	int row_num = 0;
	int column_num = 0;
	string buffer;
	ifstream matrix;

	// open file
	matrix.open("matrix.txt");

	// read first line
	getline(matrix, buffer);

	// load first char
	blosum62.x_aminos[j] = buffer[i];
	j++;

	// load rest of first line
	for(; j < ROW_LEN; j++)
	{
		i += 3;
		blosum62.x_aminos[j] = buffer[i];
	}

	i = 0;

	// based on input file pattern
	for(j = 0; j < ROW_LEN; j++)
	{
		i = 0;
		column_num = 0;
		getline(matrix, buffer);
		blosum62.y_aminos[j] = buffer[i];
		for(i = 2; i < buffer.size();)
		{
			if(buffer[i] == '-')
			{
				i++;
				blosum62.scores[row_num][column_num] = 0 - atoi(&buffer[i]);
			}
			// handles case w matched with w since that is the only double digit score
			else if(isdigit(buffer[i]))
			{
				string temp;
				temp += buffer[i];
				i++;
				temp += buffer[i];
				blosum62.scores[row_num][column_num] = atoi(temp.c_str());
			}
			else
			{
				i++;
				blosum62.scores[row_num][column_num] = atoi(&buffer[i]);
			}
			column_num++;
			i += 2;
		}
		row_num++;
	}

	matrix.close();
}


// get data sequence
string Extract_Seq(string file)
{
	ifstream f;

	f.open(file.c_str());

	string s;

	getline(f, s, '\n');

    f.close();

   	return s;
}


// computes best possible alignment for the 2 data sequences
// Returns Possibilities struct containing 2D traceback table
Possibilities Compute_Poss(Scores_Matrix blosum62, string sequence1, string sequence2)
{
	int max;
	int score;
	int val_d;
	int val_up;
	int val_left;
	char dir;
	int i = 0;
	int j = 0;
	int posx = 0;
	int posy = 0;
	Possibilities P;

	// allows for x-axis to always be longest
	if(sequence1.length() >= sequence2.length())
	{
		P.cols = sequence1.length() + 1;		//x-axis
		P.rows = sequence2.length() + 1;		//y-axis
		P.original1 = sequence1;
		P.original2 = sequence2;
		P.file1_first = true;
	}
	else
	{
		P.cols = sequence2.length() + 1;		//x-axis
		P.rows = sequence1.length() + 1;		//y-axis
		P.original1 = sequence2;
		P.original2 = sequence1;
		P.file1_first = false;
	}

	// init results matrix and movements
	P.results = new int*[P.rows];
	P.movements = new char*[P.rows];

	for(int k = 0; k < P.rows; k++)
	{
		P.results[k] = new int[P.cols];
		P.movements[k] = new char[P.cols];
	}
	
	// always 0
	P.results[0][0] = 0;

	// x marks end for backtrace
	P.movements[0][0] = 'x';

	
	// init first line row of x-axis, only go left
	for(int k = 1; k < P.cols; k++)
	{
		P.results[0][k] = GAP_PENALTY * k;
		P.movements[0][k] = 'l';
	}

	// init first column of y-axis, only up
	for(int k = 1; k < P.rows; k++)
	{
		P.results[k][0] = GAP_PENALTY * k;
		P.movements[k][0] = 'u';
	}

 		
	// compare one element in y(shortest word) to all in x(longest word) first
	for(i = 1; i < P.rows; i++)
	{
		for(j = 1; j < P.cols; j++)
		{
			val_d = P.results[i-1][j-1];
			val_up = P.results[i-1][j];
			val_left = P.results[i][j-1];

			posx = Find_Pos(blosum62, P.original1[j-1]);
			posy = Find_Pos(blosum62, P.original2[i-1]);

			score = blosum62.scores[posx][posy];

			max = Find_Greatest_Val(val_d, val_up, val_left, score, dir);

			P.results[i][j] = max;
			P.movements[i][j] = dir;
 		}
	}
	return P;
}

// since both lists are the same only have to look at
// x_axis for both x and y 
int Find_Pos(Scores_Matrix blosum62, char pos)
{
	for(int i = 0; i < ROW_LEN; i++)
	{
		if(blosum62.x_aminos[i] == pos)
		{
			return i;
		}
	}

}


// computes the correct direction to take for the backtrace of score
int Find_Greatest_Val(int val_d, int val_up, int val_left, int score, char &dir)
{
	int d;
	int up;
	int left;

	d = val_d + score;
	up = val_up + GAP_PENALTY;
	left = val_left + GAP_PENALTY;

	// it is >= because if d == up we want to take diagonal
	if(d >= up)
	{
		if(d >= left)
		{
			dir = 'd';
			return d;
		}
	}
	
	// give up higer precedence to left take == as up
	if(up >= left)
	{
		dir = 'u';
		return up;
	}
	else
	{
		dir = 'l';
		return left;
	}
}


// determines best alignment
void Traceback(Possibilities &P)
{
	int i = P.rows-1;
	int j = P.cols-1;
	int pos1 = P.original1.length()-1;		// for seq1
	int pos2 = P.original2.length()-1;		// for seq2

	while(i >= 0 && j >= 0)
	{
		// found end of strings
		if(P.movements[i][j] == 'x')
		{
			P.optimal1 = Reverse_Str(P.optimal1);
			P.optimal2 = Reverse_Str(P.optimal2);
			return;
		}
		//diagonal
		else if(P.movements[i][j] == 'd') 
		{ 
			P.optimal1 += P.original1[pos1];
			P.optimal2 += P.original2[pos2];
			pos1--;
			pos2--;
			i--;
			j--;
		}
		//up
		else if(P.movements[i][j] == 'u')
		{
			P.optimal1 += "-";
			P.optimal2 += P.original2[pos2];
			pos2--;
			i--;
		}
		//left
		else
		{
			P.optimal1 += P.original1[pos1];
			P.optimal2 += "-";
			pos1--;
			j--;
		}
	}
}


// reverse optimal string to correct the order
string Reverse_Str(string optimal)
{
	int i = optimal.length()-1;
	string temp;

	for(; i >= 0; i--)
	{
		temp += optimal[i];
	}

	return temp;
}