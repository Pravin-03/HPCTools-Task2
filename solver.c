#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double *generate_matrix(int size)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);  //Pas un pointeur vers une matrice mais un pointeur vers vecteur colonne
    srand(1); // repete les valeurs aléatoires

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
}

void print_matrix(const char *name, double *matrix, int size)
{
    int i, j;
    printf("matrix: %s \n", name);

    for (i = 0; i < size; i++)
    {
            for (j = 0; j < size; j++)
            {
                printf("%f ", matrix[i * size + j]);
            }
            printf("\n");
    }
}

int check_result(double *bref, double *b, int size) {
    int i;
    for(i=0;i<size*size;i++) {
        if (bref[i]!=b[i]) return 0;
    }
    return 1;
}

void transp(double *matrix,int size)         //transposition of a matrix
{
	double temp;
	int i,j;
	for(i=0;i<size;i++)
	{
		for(j=0;j<i;j++)
		{
			temp=matrix[i*size+j]; //Because as it is defined, A is not a matrix n*n but a vector with n² rows
			matrix[i*size+j]=matrix[j*size+i];
			matrix[j*size+i]=temp;
		}
	}
}

double *mult(double *A,double *B,int size)  //multiplication of 2 matrix
{
	int i,j,k;
	double *C = (double *)malloc(sizeof(double) * size * size); //initialisaiton of the matrix multiplication
	for(i=0;i<size*size;i++)
	{
		C[i]=0.0;
	}
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			for(k=0;k<size;k++)
			{
				C[i*size+j]=C[i*size+j]+A[i*size+k]*B[k*size+j];
			}
		}
	}
	return C;
}

double *mat_rot(int size, int i, int j, double c, double s) // for QR algorithm
{
	int k;
	double *M = (double *)malloc(sizeof(double) * size * size);
	for(k=0;k<size*size;k++)
	{
		M[k]=1;
	}
	M[i*size+i]=c;
	M[j*size+j]=c;
	M[i*size+j]=-s;
	M[j*size+i]=s;
	return M;
}

void QR_givens(double *Q, double *R, int size) // QR method with Givens rotation
{
	int j,i;
	for(j=0;j<size-1;j++)
	{
		for(i=size-1;i>j;i--)
		{
			int position1=i;
			int position2=i-1;
			double r;
			r=sqrt(R[position1*size+j]*R[position1*size+j]+R[position2*size+j]*R[position2*size+j]);
			double *G;
			G=generate_matrix(size);
			G=mat_rot(size,position1,position2,R[position2*size+j]/r,R[position1*size+j]/r);
			//print_matrix("G",G,size);
			R=mult(G,R,size);
			transp(G,size);                                     // in order to use it in the following multiplication
			Q=mult(Q,G,size);
		}
	}
	if(R[(size-1)*size+(size-1)]<0)
	{
		R[(size-1)*size+(size-1)]=R[(size-1)*size+(size-1)]*(-1);
		for(i=0;i<size;i++)
		{
			Q[i*size+2]=Q[i*size+2]*(-1);
		}
	}
	//print_matrix("Rint",R,size);
	//print_matrix("Qint",Q,size);
}

double *resolve(double *R,double *B, int size) //solve AX=B
{
	double sum;
	int i,j,k;
	double *X = (double *)malloc(sizeof(double) * size * size);
	for(i=0;i<size*size;i++)
	{
		X[i]=0.0;
	}
	for(i=size-1;i>=0;i--)
	{
		for(j=0;j<size;j++)
		{
			if (i==size-1)
			{
				if (R[i*size+i]!=0)
				{
					X[i*size+j]=B[i*size+j]/R[i*size+i];
				}
				else
				{
					X[i*size+j]=0;
				}
			}
			else
			{
				sum=0.0;
				for(k=size-1;k>i;k--)
				{
					sum=sum+R[i*size+k]*X[k*size+i];
				}
				if (R[i*size+i]!=0)
				{
					X[i*size+j]=(B[i*size+j]-sum)/R[i*size+j];
				}
				else
				{
					X[i*size+j]=0.0;
				}
			}
		}
	}
	return X;
}

void main(int argc, char *argv[])
{
   int size = atoi(argv[1]);

   double *A, *Aref;
   double *B, *Bref;

   A = generate_matrix(size);
   Aref = generate_matrix(size);        
   B = generate_matrix(size);
   Bref = generate_matrix(size);
        
   //print_matrix("A", A, size);
   //print_matrix("B", B, size);
	
	int i;

	double *R,*Rref;
	double *Q,*Qref;
	
	R=generate_matrix(size);
	Q=generate_matrix(size);

	
	for(i=0;i<size*size;i++) // Initialization of Q and R
	{
		R[i]=A[i];
		Q[i]=0.0;
	}
	for(i=0;i<size;i++)
	{
		Q[i*size+i]=1.0;
	}
	
	//print_matrix("Q before",Q,size);
	//print_matrix("R before",R,size);
	QR_givens(Q,R,size);                       // Don't know why but Q and R stay the same after QR_givens
	//print_matrix("Q after",Q,size);
	//print_matrix("R after",R,size);
	transp(Q,size);
	double *C;
	C=mult(Q,B,size);
	double *X;
	X=resolve(R,C,size);
	//print_matrix("X",X,size);
}