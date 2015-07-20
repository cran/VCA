/************************************************************************************************
 *																								*
 * Compute A- and B-matrix in the Square Root and Abbreviated Doolittle methods for computing   *
 * the coefficient matrix C in ANOVA Type 1 estimation of variance components.          		*
 * Implementation of the SWEEP-Operator for obtaining a generalized inverse of a matrix and for *
 * computing ANOVA sum of squares.																*
  *				  																			    *			   																			        *
 * Author:			Dr. Andre Schuetzenmeister											        *
 *																								*
 *					Roche Diagnostics GmbH														*
 *					DXREBA       																*
 *					Nonnenwald 2																*
 *					82377 Penzberg / Germany													*
 *																								*
 *					Phone: +49 8856 60 6621														*
 *					mailto:andre.schuetzenmeister@roche.com										*
 *																								*
 *				  																			    *
 * Last modified:	2015-May-19																	*
 *																								*
 ************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

void getAmatBmatSparse(	double* xEl, int* iEl, int* pEl, int* NobsCol, int* ncol, int* nrow, double* tol, 
						double* Amat, double* Bmat, double* amat, double* Cmat, int* NVC, int* asgn, 
						int* DF);
						
void getAmatBmat(	double* mm, int* ncol, int* nrow, double* tol, double* Amat, double* Bmat, double* amat,
					double* Cmat, int* NVC, int* asgn, int* DF);
					
void Tsweep(double* M, int* k, double* thresh, int* NumK, int* nr, int* LC, double* tol, double *SSQ);

					
void TsweepFull(double* M, int* nr, double* tol);

double Round(double val, double Pow);

double Round(double val, double Pow)
{
	val = val * Pow;
	val = round(val);
	val = val/Pow;
	return val;
}


/*
	Implements the "Square Root and Abbreviated Doolittle Method", specifically, computes
	matrices 'Amat' and 'Bmat' from which the coefficient matrix C is derived.
	
	This is the dense-matrix implementation of the algorithm.
	
	mm			(double*) model matrix of the full model
	ncol		(int*) pointer to integer specifying the number of columns of 'mm'
	nrow		(int*) pointer to integer specifying the number of rows of 'mm'
	tol			(double*) pointer to double holding the  numerical tolerance values used to distinguish a value from zero
	Amat		(double*) pointer to double array where 'Amat' is stored
	Bmat		(double*) pointer to double array where 'Bmat' is stored
	amat		(double*) pointer to double array where 'amat' is stored
	Cmat		(double*) pointer to double array where matrix C is stored
	NVC			(int*) pointer to integer storing the (n)umber of (v)ariance (c)omponents
	asgn		(int*) pointer to integer array specifying the columns (pEl) where i-th (v)ariance (c)omponent is stored as integer 'i'
	DF			(int*) pointer to integer array specifying the ANOVA (d)egrees of (f)reedom
*/

void getAmatBmat(	double* mm, int* ncol, int* nrow, double* tol, double* Amat, double* Bmat, double* amat,
					double* Cmat, int* NVC, int* asgn, int* DF)
{
	double tmp;
	double Bmu[(*ncol)+1];
	double Amu[(*ncol)+1];	
	double rsums[(*ncol)+1];								// vector of row sums of amat
	int i,j,k,l, sgn;
	*ncol = *ncol - 1;
	double Pow = 1/(*tol);
	

	for(i=0;i<((*ncol)+1);i++)						// over columns of 'mm' 
	{
		tmp = 0;
		for(j=0;j<(*nrow);j++)
		  	tmp = tmp + mm[i*(*nrow)+j];
		Amu[i] = tmp;								// Amu contains columns sums
		if(i == 0)
			Bmu[0] = 1;								// Bmu contais standardized columns sums
		else
		    Bmu[i] = tmp/(*nrow);
	}

	for(i=0; i<(*ncol); i++)						// each column ...	
	{
		for(j=i; j<(*ncol); j++)					// ... with each other column	
		{
			tmp = 0;
			
			for(k=0; k<(*nrow); k++)									// multiply k-th element of col (i+1) with
			{
				tmp = tmp + mm[(i+1)*(*nrow)+k] * mm[(j+1)*(*nrow)+k];		// k-th element of col (j+1) and sum up results
			}
			tmp = tmp - (Amu[i+1] * Bmu[j+1]);
			Amat[j*(*ncol)+i] = tmp;
			
			if(i > 0)
			{
				tmp = 0;
				
				for(k=0; k<i; k++)
				{
					tmp = tmp + (Amat[i*(*ncol)+k] * Bmat[j*(*ncol)+k]);
				}
					
				Amat[j*(*ncol)+i] = Amat[j*(*ncol)+i] - tmp;
			}
			
			if((i == j) && (Amat[i*(*ncol)+i] < *tol))		// set all elements to zero if diagonal element is zero
			{
				for(k=0; k<(*ncol); k++)
				{
					Amat[k*(*ncol)+i] = 0;
					Bmat[k*(*ncol)+i] = 0;
					amat[k*(*ncol)+i] = 0;
				}
				break;								// break j-loop 
			}
			
			Bmat[j*(*ncol)+i] = Amat[j*(*ncol)+i] / Amat[i*(*ncol)+i];			
			amat[j*(*ncol)+i] = Round(Amat[j*(*ncol)+i], Pow) * Round(Bmat[j*(*ncol)+i], Pow);			
			sgn = Amat[j*(*ncol)+i] < 0 ? -1 : 1;			
			amat[j*(*ncol)+i] = sgn * sqrt(amat[j*(*ncol)+i]);			
		}
	}

	// now build elements of the C-matrix

	for(i=0;i<(*NVC);i++)							// over rows of C
	{
		for(j=i;j<(*NVC);j++)						// over columns of C
		{
			for(k=0;k<(*ncol);k++)					// over rows of amat
			{
				if( (rsums[k] < *tol) & (rsums[k] > Pow) )		// row sum is zero
					break;
					
				for(l=0;l<(*ncol);l++)				// over columns of amat
				{
					if((asgn[k]-1)==i && (asgn[l]-1==j))
					{					
						Cmat[j*(*NVC)+i] = Cmat[j*(*NVC)+i] + amat[l*(*ncol)+k]*amat[l*(*ncol)+k];		// square of that element
					}
				}
			}
			Cmat[j*(*NVC)+i] = Cmat[j*(*NVC)+i] / DF[i];
		}
	}
			
	return;
}

/*
	Implements the "Square Root and Abbreviated Doolittle Method", specifically, computes
	matrices 'Amat' and 'Bmat' from which the coefficient matrix C is derived.
	
	This is the sparse-matrix implementation of the algorithm.
	
	xEl			(double*) vector storing elements of the sparse matrix
	iEl			(int*) pointer to integer vector storing non-zero elements of the sparse-matrix (see R-class CsparseMatrix)
	pEl			(int*) pointer to integer vector storing the index in 'x' of the first non-zero element per column
	NobsCol		(int*) pointer to integer vector storing the number of non-zero elements per column
	ncol		(int*) pointer to integer specifying the number of columns of the sparse matrix
	nrow		(int*) pointer to integer specifying the number of rows of the sparse matrix
	tol			(double*) pointer to double holding the  numerical tolerance values used to distinguish a value from zero
	Amat		(double*) pointer to double array where 'Amat' is stored
	Bmat		(double*) pointer to double array where 'Bmat' is stored
	amat		(double*) pointer to double array where 'amat' is stored
	Cmat		(double*) pointer to double array where matrix C is stored
	NVC			(int*) pointer to integer storing the (n)umber of (v)ariance (c)omponents
	asgn		(int*) pointer to integer array specifying the columns (pEl) where i-th (v)ariance (c)omponent is stored as integer 'i'
	DF			(int*) pointer to integer array specifying the ANOVA (d)egrees of (f)reedom
	 
*/

void getAmatBmatSparse(	double* xEl, int* iEl, int* pEl, int* NobsCol, int* ncol, int* nrow, 
						double* tol, double* Amat, double* Bmat, double* amat,
						double* Cmat, int* NVC, int* asgn, int* DF)
{
	double tmp, root=0.;
	double Bmu[(*ncol)+1];
	double Amu[(*ncol)+1];	
	int i,j,k,l,sgn;

	*ncol = *ncol - 1;
	double Pow = 1/(*tol);


	for(i=0; i<((*ncol)+1); i++)							// over columns of 'mm' 
	{
		tmp = 0;
		for(j=0;j<NobsCol[i];j++)
		  	tmp = tmp + xEl[pEl[i]+j];
		Amu[i] = tmp;										// Amu contains columns sums
		if(i == 0)
			Bmu[0] = 1;										// Bmu contais standardized columns sums
		else
		    Bmu[i] = tmp/(*nrow);
	}
	
	for(i=0; i<(*ncol); i++)								// each column ...	
	{
		for(j=i; j<(*ncol); j++)							// ... with each other column	
		{
			tmp = 0.;
			
			for(k=0; k<NobsCol[i+1]; k++)											// multiply k-th element of col (i+1) with
			{				
				for(l=0; l<NobsCol[j+1]; l++)
				{
					if(iEl[pEl[i+1]+k] < iEl[pEl[j+1]])								// if row index in col i < 1st row index of col j leave loop 
						break;
								
					if(iEl[pEl[i+1]+k] == iEl[pEl[j+1]+l])							// both columns have non-zero element in same row
					{
						tmp = tmp + xEl[pEl[i+1]+k] * xEl[pEl[j+1]+l];				// multiply them and add result
					}					
				}
			}
			
			tmp = tmp - (Amu[i+1] * Bmu[j+1]);
			Amat[j*(*ncol)+i] = tmp;

			if(i > 0)
			{
				tmp = 0;
				
				for(k=0; k<i; k++)
					tmp = tmp + (Amat[i*(*ncol)+k] * Bmat[j*(*ncol)+k]);
					
				Amat[j*(*ncol)+i] = Amat[j*(*ncol)+i] - tmp;
			}
			
			if((i == j) && (Amat[i*(*ncol)+i] < *tol))		// set all elements to zero if diagonal element is zero
			{
				for(k=0; k<(*ncol); k++)
				{
					Amat[k*(*ncol)+i] = 0;
					Bmat[k*(*ncol)+i] = 0;
					amat[k*(*ncol)+i] = 0;
				}
				break;								// break j-loop 
			}
			
			Bmat[j*(*ncol)+i] = Amat[j*(*ncol)+i] / Amat[i*(*ncol)+i];			
			amat[j*(*ncol)+i] = Round(Amat[j*(*ncol)+i], Pow) * Round(Bmat[j*(*ncol)+i], Pow);			
			sgn = Amat[j*(*ncol)+i] < 0 ? -1 : 1;		
			root = sqrt(amat[j*(*ncol)+i]);	
			amat[j*(*ncol)+i] = sgn * root;		
			
			//rsums[i] = rsums[i] + root;
		}		
	}

	// now build elements of the C-matrix
	
	for(i=0;i<(*NVC);i++)							// over rows of C
	{
		for(j=i;j<(*NVC);j++)						// over columns of C
		{
			for(k=0;k<(*ncol);k++)					// over rows of amat
			{									
				for(l=0;l<(*ncol);l++)				// over columns of amat
				{
					if((asgn[k]-1)==i && (asgn[l]-1==j))
					{					
						Cmat[j*(*NVC)+i] = Cmat[j*(*NVC)+i] + amat[l*(*ncol)+k]*amat[l*(*ncol)+k];		// square of that element
					}
				}
			}
			Cmat[j*(*NVC)+i] = Cmat[j*(*NVC)+i] / DF[i];
		}
	}
		
	return;
}






/* 
   Implementation of the sweep operator operating on the transpose of matrix 'M'. This ensures that rows of M
   are stored in consecutive parts of the memory speeding up sweeping a bit. 
   Function sweeps complete matrix and updates the vector of ANOVA sum of squares SSQ and the vector indicating
   linearily depending columns LC.
   The sweep-operator is modified leaving out operations on those elements of 'M' not required for 
   computing the ANOVA sum of squares. Thus, the result is not as described in Goodnight (1979), generating
   the inverted matrix M and the vector of fixed effects b_hat.
   
   M			(double*) pointer to first element to double array storing the matrix elements
   k			(int*) pointer to first element of integer array, value 'i' corresponds to i-th variance component, 
   				indicating which columns belong to which variance component (variable)
   thresh		(double*) pointer to double storing the threshold value used to check for linear dependencies among
   				columns of M
   NumK			(int*) pointer to integer storing number of columns to consider
   nr			(int*) pointer to integer storing the number of rows (columns) of the matrix
   LC			(int*) pointer to first element of an integer vector indicating linear dependence of the i-th column
                to preceeding columns as 1
   tol			(double*) pointer to double value specifying the numerical equivalence to zero to be used
   SSQ			(double*) pointer to double array where ANOVA Type-I sum of squares will be stored, i-th element
   				corresponds to i-th variable in the model (with one or multiple elements)
*/

void Tsweep(double* M, int* k, double* thresh, int* NumK, int* nr, int* LC, double* tol, double* SSQ)
{
	double ESS, CSS=0., Rsq, D, B;
	int offset, i, j, n, l=0, NN;					// l runs over variables
	
	NN=(*nr) * (*nr) - 1;

	for(n=0; n<*NumK; n++)							// over columns k to be swept
	{
		offset = n*(*nr);

		if( fabs(M[offset+n]) < *tol )
		{				
			LC[n] = 1;								// indicate linear dependency of current column
			
			if(n == (*NumK-1) || k[n] != k[n+1])	// next variable
			{
				SSQ[l] = CSS;
				l += 1;
			}
			continue;	
		}
		
		ESS = M[offset+n];							// M[k,k]
		D = ESS;
	
		for(i=n; i<*nr; i++)						// adapt all rows
		{					
		  	if(i == n)
		  	{
		  		for(j=n; j<*nr; j++)
		  		{
					M[offset+j] = M[offset+j]/D;
				}
		  	}
		  	else
		  	{
		  		B = M[i*(*nr)+n];
		  		if(fabs(B) < *tol)					// won't do anything
		  			continue;
		  				  		
		  		for(j=n; j<*nr; j++)
		  		{
						M[i*(*nr)+j] = M[i*(*nr)+j] - B * M[n*(*nr)+j];
				}
			}
		}

		CSS = M[NN];								// M[nr,nr]
		Rsq = (CSS - ESS) / CSS;
		if(Rsq > (1-*thresh))
		{
			LC[n] = 1;
		}
		
		if(n == (*NumK-1) || k[n] != k[n+1])		// either last column overall or corresponding to specific variable
		{
			SSQ[l] = CSS;
			l += 1;
		}
	}
	return;
}


/* 
   Implementation of the sweep operator for matrix inversion operating on the transpose of matrix 'M'. 
   This ensures that rows of M are stored in consecutive parts of the memory speeding up sweeping a bit. 
      
   M			(double*) pointer to first element to double array storing the matrix elements
   nr			(int*) pointer to integer storing the number of 
   double		(double*) pointer to double value specifying the numerical equivalence to zero
*/

void TsweepFull(double* M, int* nr, double* tol)
{
	double D, B;
	int offset, i, j, n;					// l runs over variables
	
	for(n=0; n<*nr; n++)							// over columns k to be swept
	{
		offset = n*(*nr);

		if( fabs(M[offset+n]) < *tol )
		{	
			for(j=0; j<*nr; j++)
			{
			  	M[j*(*nr)+n] = 0;					// set k-th t(row) zero
			  	M[offset+j] = 0;					// set k-th t(column) zero
			}
			
			continue;	
		}
				
		D = M[offset+n];
		
		for(j=0; j<*nr; j++)						// adjust n-th row
  		{
			M[offset+j] = M[offset+j]/D;
		}

		for(i=0; i<*nr; i++)						// adapt all rows
		{					
		  	if(i != n)
		  	{		  		
		  		B = M[i*(*nr)+n];
		  		if(fabs(B) < *tol)
		  			continue;
		  		
		  		for(j=0; j<*nr; j++)
		  		{
						M[i*(*nr)+j] = M[i*(*nr)+j] - B * M[n*(*nr)+j];
				}
				M[i*(*nr)+n] = (-1)*B/D;
			}
		}
		M[offset+n] = 1/D;
	}
	return;
}

