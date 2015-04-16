/************************************************************************************************
 *																								*
 * Compute A- and B-matrix in the Square Root and Abbreviated Doolittle methods for computing   *
 * the coefficient matrix C in ANOVA Type 1 estimation of variance components.          		*
 *			   																			        *
 * Author:			Dr. André Schützenmeister											        *
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
 * Last modified:	2015-04-07																    *
 *																								*
 ************************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


void getAmatBmat(double* mm, int* ncol, int* nrow, double* tol, double* Amat, double* Bmat);

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
