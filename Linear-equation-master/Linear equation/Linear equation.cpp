#pragma once
#include <stdio.h>
#include <iostream>
#include <vector>

using namespace std;

typedef struct Matrix
{
public:
	int maxRow;
	int maxCol;
	vector<vector<float>> matrixdata;
	Matrix(int row, int col)
	{
		//matrixdata의 내부를 num의 사이즈로 재설정하며, 0으로 초기화함
		matrixdata.resize(col, vector<float>(row, 0));
		maxRow = row;
		maxCol = col;
	}
}MATRIX;

void getCofactor(Matrix& A, Matrix& temp, int p, int q, int n)
{
	int i = 0, j = 0;

	// Looping for each element of the matrix 
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			//  Copying into temporary matrix only those element 
			//  which are not in given row and column 
			if (row != p && col != q)
			{
				temp.matrixdata[i][j++] = A.matrixdata[row][col];

				// Row is filled, so increase row index and 
				// reset col index 
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
int determinant(Matrix& A, int n)
{
	int D = 0; // Initialize result

	//  Base case : if matrix contains single element
	if (n == 2)
		return ((A.matrixdata[0][0] * A.matrixdata[1][1]) - (A.matrixdata[1][0] * A.matrixdata[0][1]));

	Matrix temp(A.maxRow, A.maxCol); // To store cofactors

	int sign = 1;  // To store sign multiplier

	 // Iterate for each element of first row
	for (int f = 0; f < n; f++)
	{
		// Getting Cofactor of A[0][f] 
		getCofactor(A, temp, 0, f, n);
		D += sign * A.matrixdata[0][f] * determinant(temp, n - 1);

		// terms are to be added with alternate sign 
		sign = -sign;
	}

	return D;
}

// Function to get adjoint of A[N][N] in adj[N][N]. 
void adjoint(Matrix& A, Matrix& adj)
{
	if (A.matrixdata.size() == 1)
	{
		adj.matrixdata[0][0] = 1;
		return;
	}

	// temp is used to store cofactors of A[][] 
	int sign = 1;

	Matrix temp(A.maxRow, A.maxCol);

	for (int i = 0; i < A.matrixdata.size(); i++)
	{
		for (int j = 0; j < A.matrixdata.size(); j++)
		{
			// Get cofactor of A[i][j] 
			getCofactor(A, temp, i, j, A.matrixdata.size());

			// sign of adj[j][i] positive if sum of row 
			// and column indexes is even. 
			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the 
			// transpose of the cofactor matrix 
			adj.matrixdata[j][i] = (sign) * (determinant(temp, A.matrixdata.size() - 1));
		}
	}
}

// Function to calculate and store inverse, returns false if 
// matrix is singular 
bool inverse(Matrix& A, Matrix& inverse)
{
	// Find determinant of A[][] 
	int det = determinant(A, A.matrixdata.size());
	if (det == 0)
	{
		cout << "Singular matrix, can't find its inverse";
		return false;
	}

	// Find adjoint 
	Matrix adj(A.maxRow, A.maxCol);

	adjoint(A, adj);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
	for (int i = 0; i < A.matrixdata.size(); i++)
		for (int j = 0; j < A.matrixdata.size(); j++)
			inverse.matrixdata[i][j] = adj.matrixdata[i][j] / float(det);

	return true;
}

// Generic function to display the matrix.  We use it to display 
// both adjoin and inverse. adjoin is integer matrix and inverse 
// is a float. 
void display(Matrix A)
{
	for (int i = 0; i < A.matrixdata.size(); i++)
	{
		for (int j = 0; j < A.matrixdata.size(); j++)
			cout << A.matrixdata[i][j] << " ";
		cout << endl;
	}
}

void main()
{
	Matrix A(3, 3);
	A.matrixdata = { {1, 0, 5},
					 {2, 1, 6},
					 {3, 4, 0} };

	Matrix adj(3, 3);  // To store adjoint of A[][] 
			  
	Matrix inv(3, 3); // To store inverse of A[][] 

	//Matrix A(4, 4);
	//A.matrixdata = { {5, -2, 2, 7},
	//				{1, 0, 0, 3},
	//				{-3, 1, 5, 0},
	//				{3, -1, -9, 4} };

	//Matrix adj(4, 4);  // To store adjoint of A[][] 
	//		  
	//Matrix inv(4, 4); // To store inverse of A[][] 

	cout << "Input matrix is :\n";
	display(A);

	//cout << determinant(A, A.matrixdata.size());

	cout << "\nThe Adjoint is :\n";
	adjoint(A, adj);
	display(adj);

	//cout << "\nThe Inverse is :\n";
	//if (inverse(A, inv))
	//	display(inv);
}