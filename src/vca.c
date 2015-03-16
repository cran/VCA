/************************************************************************************************
 *																								*
 * Compute A- and B-matrix in the Square Root and Abbreviated Doolittle methods for computing   *
 * the coefficient matrix C in ANOVA Type 1 estimation of variance components.          		*
 *			   																			        *
 * Author:			Dr. André Schützenmeister											        *
 *					SW Engineering & Data Processing											*
 *																								*
 *					Roche Diagnostics GmbH														*
 *					DXREBC..6164																*
 *					Nonnenwald 2																*
 *					82377 Penzberg / Germany													*
 *																								*
 *					Phone: +49 8856 60 6621														*
 *					mailto:andre.schuetzenmeister@roche.com										*
 *																								*
 *				  																			    *
 * Last modified:	2014-July-18																*
 *																								*
 ************************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>
#ifndef PI
#define PI 3.141592653589793115998
#endif
typedef enum {false, true} bool;
#ifndef EPS
#define EPS 1.e-12		/* precision to be used, not used right now */
#endif

void getAmatBmat(double* mm, int* ncol, int* nrow, double* tol, double* Amat, double* Bmat);
void multVectors(double* mm1, double* mm2, int* row1, int* row2, int* ind1, int* ind2, int* len, double* res);

/*
	mm1			(double*) pointer to double array representing matrix 1
	mm2			(double*) pointer to double array representing matrix 2
	row1		(int*) pointer to integer indicating whether i-th row (1) or j-th column (0) of 'mm1' should be used
	row2		(int*) pointer to integer indicating whether i-th row (1) or j-th column (0) of 'mm2' should be used
	ind1		(int*) pointer to integer specifying the i-th row or j-th column (depends) on 'row1') of 'mm1'
	ind2		(int*) pointer to integer specifying the i-th row or j-th column (depends) on 'row2') of 'mm2'
	len			(int*) pointer to integer specifying the length of both vectors to be matrix multiplied
	res			(double*) pointer to double where the result should be stored
*/

void multVectors(double* mm1, double* mm2, int* row1, int* row2, int* ind1, int* ind2, int* len, double* res)
{
	double tmp = 0;
	
	for(int i=0; i<(*len); i++)
	{
		if(*row1 == 1 && *row2 == 1)				/* both rows */
		{
			tmp = tmp + mm1[i*(*len)+(*ind1)] * mm2[i*(*len)+(*ind2)];
		}
		else if(*row1 != 1 && *row2 != 1) 			/* both columns */
		{
			tmp = tmp + mm1[(*ind1)*(*len)+i] * mm2[(*ind2)*(*len)+i];
		}
		else if(*row1 == 1 && *row2 != 1)			/* row of mm1 and column of mm2 */
		{
			tmp = tmp + mm1[i*(*len)+(*ind1)] * mm2[(*ind2)*(*len)+i];
		}
		else										/* column of mm1 and row of mm2 */
		{
			tmp = tmp + mm1[(*ind1)*(*len)+i] * mm2[i*(*len)+(*ind2)];
		}
	}
	*res = tmp;
}

/*
	Implements the "Square Root and Abbreviated Doolittle Method", specifically, computes
	matrices 'Amat' and 'Bmat' from which the coefficient matrix C is derived.
	
	mm			(double*) model matrix of the full model
	ncol		(int*) pointer to integer specifying the number of columns of 'mm'
	nrow		(int*) pointer to integer specifying the number of rows of 'mm'
	tol			(double*) pointer to double holding the  numerical tolerance values used to distinguish a value from zero
	Amat		(double*) pointer to double array where 'Amat' is stored
	Bmat		(double*) pointer to double array where 'Bmat' is stored
*/

void getAmatBmat(double* mm, int* ncol, int* nrow, double* tol, double* Amat, double* Bmat)
{
	double tmp;
	double Bmu[(*ncol)+1];
	double Amu[(*ncol)+1];	
	int i,j,k;
	
	*ncol = *ncol - 1;
	
	for(i=0;i<((*ncol)+1);i++)						/* over columns of 'mm' */
	{
		tmp = 0;
		for(j=0;j<(*nrow);j++)
		  	tmp = tmp + mm[i*(*nrow)+j];
		Amu[i] = tmp;
		if(i == 0)
			Bmu[0] = 1;
		else
		    Bmu[i] = tmp/(*nrow);
	}
	
	for(i=0; i<(*ncol); i++)
	{
		for(j=i; j<(*ncol); j++)
		{
			tmp = 0;
			for(int k=0; k<(*nrow); k++)
				tmp = tmp + mm[(i+1)*(*nrow)+k] * mm[(j+1)*(*nrow)+k];
			tmp = tmp - (Amu[i+1] * Bmu[j+1]);
			Amat[j*(*ncol)+i] = tmp;
			
			if(i > 0)
			{
				tmp = 0;
				
				for(k=0; k<i; k++)
					tmp = tmp + (Amat[i*(*ncol)+k] * Bmat[j*(*ncol)+k]);
					
				Amat[j*(*ncol)+i] = Amat[j*(*ncol)+i] - tmp;
			}
			
			if((i == j) && (Amat[i*(*ncol)+i] < *tol))
			{
				for(k=0; k<(*ncol); k++)
				{
					Amat[k*(*ncol)+i] = 0;
					Bmat[k*(*ncol)+i] = 0;
				}
				break;								/* break j-loop */
			}
			
			Bmat[j*(*ncol)+i] = Amat[j*(*ncol)+i] / Amat[i*(*ncol)+i];
		}
	}
	return;
}
