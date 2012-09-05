#ifndef MY_MEM_H
#define MY_MEM_H

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string>
using std::string;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//							Templates get mem      
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


template <class VarType>
void getMem1D(VarType **ptr, int i1, int i2, string name, bool setZero = false)
{
	VarType *ar = new VarType[i2-i1+1];

	if(ar == NULL)
		throw(-1);

	if(setZero)
	  for (int i = 0; i < i2-i1+1; ++i) 
	    ar[i] = (VarType)0.0;
		//memset(ar, 0, sizeof(VarType)*(i2-i1+1));

	ar -= i1;
	*ptr = ar;
}

template <class VarType>
void getMem2D(VarType ***ptr, int i1, int i2, int j1, int j2, string name, bool setZero = false)
{
	VarType **ar = new VarType*[i2-i1+1];

	if(ar == NULL)
    throw(-1);

	ar -= i1;

	for( int i=i1; i<=i2; i++ )
		getMem1D(&ar[i],j1,j2,name,setZero);

	*ptr = ar;
}

template <class VarType>
void getMem3D(VarType ****ptr, int i1, int i2, int j1, int j2, int k1, int k2, string name, bool setZero = false)
{
	VarType ***ar = new VarType**[i2-i1+1];

	if(ar == NULL)
    throw(-1);

	ar -= i1;

	for( int i=i1; i<=i2; i++ )
		getMem2D(&ar[i],j1,j2,k1,k2,name,setZero);

	*ptr = ar;
}

template <class VarType>
void getMem4D(VarType *****ptr, int i1, int i2, int j1, int j2, int k1, int k2, int l1, int l2, string name, bool setZero = false)
{
	VarType ****ar = new VarType***[i2-i1+1];

	if(ar == NULL)
    throw(-1);

	ar -= i1;

	for( int i=i1; i<=i2; i++ )
		getMem3D(&ar[i],j1,j2,k1,k2,l1,l2,name,setZero);

	*ptr = ar;
}

template <class VarType>
void getMem5D(VarType ******ptr, int i1, int i2, int j1, int j2, int k1, int k2, int l1, int l2, int m1, int m2, string name, bool setZero = false)
{
	VarType *****ar = new VarType****[i2-i1+1];

	if(ar == NULL)
    throw(-1);

	ar -= i1;

	for( int i=i1; i<=i2; i++ )
		getMem4D(&ar[i],j1,j2,k1,k2,l1,l2,m1,m2,name,setZero);

	*ptr = ar;
}

template <class VarType>
void getMem6D(VarType *******ptr, int i1, int i2, int j1, int j2, int k1, int k2, int l1, int l2, int m1, int m2, int n1, int n2, string name, bool setZero = false)
{
	VarType ******ar = new VarType*****[i2-i1+1];

	if(ar == NULL)
    throw(-1);

	ar -= i1;

	for( int i=i1; i<=i2; i++ )
		getMem5D(&ar[i],j1,j2,k1,k2,l1,l2,m1,m2,n1,n2,name,setZero);

	*ptr = ar;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//							Templates Free      
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template<class VarType> 
void freeMem1D(VarType *ptr, int i1, int i2)
{
	delete [] (ptr+i1);
	
}

template<class VarType> 
void freeMem2D(VarType **ptr, int i1, int i2, int j1, int j2)
{
	for( int i=i1; i<=i2; i++)
		delete [] (ptr[i] + j1);

	delete [] (ptr+i1);
}

template<class VarType> 
void freeMem3D(VarType ***ptr, int i1, int i2, int j1, int j2, int k1, int k2)
{
	for( int i=i1; i<=i2; i++)
	{
		for( int j=j1; j<=j2; j++)
			delete [] (ptr[i][j] + k1);

		delete [] (ptr[i] + j1);
	}
	delete [] (ptr+i1);
}

template<class VarType> 
void freeMem4D(VarType ****ptr, int i1, int i2, int j1, int j2, int k1, int k2, int l1, int l2)
{
	for( int i=i1; i<=i2; i++)
	{
		for( int j=j1; j<=j2; j++)
		{
			for( int k=k1; k<=k2; k++)
				delete [] (ptr[i][j][k] + l1);

			delete [] (ptr[i][j] + k1);
		}
		delete [] (ptr[i] + j1);
	}
	delete [] (ptr+i1);
}

template<class VarType> 
void freeMem5D(VarType *****ptr, int i1, int i2, int j1, int j2, int k1, int k2, int l1, int l2, int m1, int m2)
{
	for( int i=i1; i<=i2; i++)
	{
		for( int j=j1; j<=j2; j++)
		{
			for( int k=k1; k<=k2; k++)
			{
				for( int l=l1; l<=l2; l++)
					delete [] (ptr[i][j][k][l] + m1);

				delete [] (ptr[i][j][k] + l1);
			}
			delete [] (ptr[i][j] + k1);
		}
		delete [] (ptr[i] + j1);
	}
	delete [] (ptr+i1);
}

template<class VarType> 
void freeMem6D(VarType ******ptr, int i1, int i2, int j1, int j2, int k1, int k2, int l1, int l2, int m1, int m2, int n1, int n2)
{
	for( int i=i1; i<=i2; i++)
	{
		for( int j=j1; j<=j2; j++)
		{
			for( int k=k1; k<=k2; k++)
			{
				for( int l=l1; l<=l2; l++)
				{
					for( int m=m1; m<=m2; m++)
						delete [] (ptr[i][j][k][l][m] + n1);

					delete [] (ptr[i][j][k][l] + m1);
				}
				delete [] (ptr[i][j][k] + l1);
			}
			delete [] (ptr[i][j] + k1);
		}
		delete [] (ptr[i] + j1);
	}
	delete [] (ptr+i1);
}

#endif




