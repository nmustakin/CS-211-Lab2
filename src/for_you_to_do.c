#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 30;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */

	int i, j, k, t;
	double temp;  
	for(i = 0; i < n-1; i++){
		/*** pivoting ***/
		int maxInd = i; 
		double max = fabs(A[i*n+i]);
		for(t = i+1; t< n; t++){
			if(fabs(A[t*n + i]) > max){
				maxInd = t;
				max = fabs(A[t*n + i]); 
			}
		} 	
		
		if(max ==0){
			printf("LUfactorization failed: coefficient matrix is singular"); 
			return -1; 
		}
		else{
			if(maxInd != i){
				/***save pivoting information***/
				int temp1 = ipiv[i];
				ipiv[i] = ipiv[maxInd];
				ipiv[maxInd] = temp1;

				/***Swap rows***/
				for(k=0; k<n; k++){
					temp = A[i*n + k];
					A[i*n +k] = A[maxInd*n + k];
					A[maxInd*n + k] = temp; 
				} 
			}
		}
		/***factorization***/
		for(j = i+1; j<n; j++){
			A[j*n+i] /= A[i*n+i];
			for(k = i+1 ; k<n; k++){
				A[j*n+k] -= A[j*n + i]*A[i*n +k];
			}
		}
	}
		

    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
	int i, j; 
	if(UPLO == 'L'){
		// A = Lower triangular matrix (matlab A) 
		// B = values (matlab b) -- save values of y here 
		// ipiv = order of x (matlab pvt)`	

		double * y = (double*)malloc(n*sizeof(double)); 
		y[0] = B[ipiv[0]]; 
		for(i = 1; i<n; i++){
			y[i] = B[ipiv[i]]; 
			for(j= 0; j<i; j++){
				y[i] -= y[j]*A[i*n + j]; 
			}
		}	
		memcpy(B, y, n*sizeof(double));
		free(y); 
	}
	else if(UPLO == 'U'){
		// A = Uppoer triangular matrix (matlab A)
		// B = values of y from L 
		// ipiv = not needed?
		B[n-1] = B[n-1]/A[(n-1)*n + n-1];
		for(i = n-2; i > -1; i--){
			for(j=i+1; j<n; j++){
				B[i] -= B[j]*A[i*n+j]; 
			}
			B[i]/= A[i*n+i];
		} 	

	}
	else{
		printf("Unrecognized option for substitution");
		return; 
	}
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */

	int i1, j1, k1;
    /* B x B mini matrix multiplications */
    for (i1 = i; i1 < i+b && i1<n; i1+=3){
	  	for (j1 = j; j1 < j+b && j1<n; j1+=3) {
      	   	register int t = i1*n+j1; 
      	    register int tt = t+n; 
      	    register int ttt = tt+n; 
      	    register double c00 = C[t]; 
     	    register double c01 = C[t+1];
      	    register double c02 = C[t+2]; 
      	    register double c10 = C[tt]; 
      	    register double c11 = C[tt+1]; 
      	    register double c12 = C[tt+2];
      	    register double c20 = C[ttt];
      	    register double c21 = C[ttt+1];
      	    register double c22 = C[ttt+2];

			for(k1=k; k1<k+b && k1<n; k1+=3){
				register int ta = i1*n+k1; 
				register int tta = ta+n;
				register int ttta = tta+n;  
				register int tb = k1*n+j1;
				register int ttb = tb+n; 
				register int tttb = ttb+n; 

				register double a00 = A[ta];
				register double a10 = A[tta];
				register double a20 = A[ttta]; 
				register double b00 = B[tb];
				register double b01 = B[tb+1];
				register double b02 = B[tb+2]; 

				c00 -= a00*b00; c01 -= a00*b01; c02 -= a00*b02;
				c10 -= a10*b00; c11 -= a10*b01; c12 -= a10*b02; 
				c20 -= a20*b00; c21 -= a20*b01; c22 -= a20*b02; 

				a00 = A[ta+1]; a10 = A[tta+1]; a20 = A[ttta+1];
				b00 = B[ttb]; b01 = B[ttb+1]; b02 = B[ttb+2];

				c00 -= a00*b00; c01 -= a00*b01; c02 -= a00*b02;
				c10 -= a10*b00; c11 -= a10*b01; c12 -= a10*b02; 
				c20 -= a20*b00; c21 -= a20*b01; c22 -= a20*b02;

				a00 = A[ta+2]; a10 = A[tta+2]; a20 = A[ttta+2];
				b00 = B[tttb]; b01 = B[tttb+1]; b02 = B[tttb+2];

				c00 -= a00*b00; c01 -= a00*b01; c02 -= a00*b02;
				c10 -= a10*b00; c11 -= a10*b01; c12 -= a10*b02; 
				c20 -= a20*b00; c21 -= a20*b01; c22 -= a20*b02;
			}

			C[t] = c00;
			C[t+1] = c01;
			C[t+2] = c02; 
			C[tt] = c10; 
			C[tt+1] = c11;
			C[tt+2] = c12; 
			C[ttt] = c20; 
			C[ttt+1] = c21;
			C[ttt+2] = c22; 
		}
 	}

    return;
}

void bijk(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    int i1, j1, k1;
    for (i1 = i; i1 < i + b && i1 < n; i1++)
    {
        for (j1 = j; j1 < j + b && j1 < n; j1++)
        {
            register double r = C[i1 * n + j1];
            for (k1 = k; k1 < k + b && k1 < n; k1++)
            {
                r -= A[i1 * n + k1] * B[k1 * n + j1];
            }
            C[i1 * n + j1] = r;
        }
    }
	return; 
}


/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
	 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
	int ib, i, j, k, maxInd;
	double max; 

	/// loop for ib = 1 to n-1 with step b

	for(ib = 0; ib < n; ib+=b){ 
		for(i = ib; i<ib+b && i < n; i++){
			/// Pivot
			/// 1) Find pivot row k, cloumn broadcast
			/// 2) Swap rows k and i in block column, broadcast row k 
			maxInd = i; 
			max = fabs(A[i*n + i]); 
			for(j = i+1; j<n; j++){
				if(fabs(A[j*n + i]) > max){
					maxInd = j; 
					max = fabs(A[j*n+i]);
				}
			}
			if(max ==0){
				return -1;
			}
			else{
				if(maxInd != i){
					/// save pivoting info 
					int temp = ipiv[i];
					ipiv[i] = ipiv[maxInd];
					ipiv[maxInd] = temp;
					/// swap rows

					/// 5) Broadcase swap info left and right
					/// 6) Apply all row swaps to other columns
					for(k=0; k<n; k++){
						double temp = A[i*n + k];
						A[i*n +k] = A[maxInd*n + k];
						A[maxInd*n + k] = temp; 
					} 

				}

			}
			/// Factorization	
			/// 3) A(i+1:n, i) = A(i+1:n, i)/A(i,i)
			/// 4) A(i+1:n, i+1:end) -= A(i+1:n, i)*A(i, i+1:end)  /// Use dgemm here? 
			for(j = i+1; j<n; j++){
				A[j*n + i] /= A[i*n + i]; 
				for(k = i+1; k<ib+b && k<n; k++){
					A[j*n + k] -= A[j*n + i] * A[i*n + k];
				}

			}

		}
		
		/// Update A(ib:end, end+1:n) 
 		/// 8) A(ib:end, end+1:n) = LL_inv * A(ib:end, end+1:n) 
		for(i = ib; i< ib+b && i < n; i++){
  			for(k = ib; k<i; k++){
				register double r = A[i*n + k];
				for(j = ib+b; j<n; j++){
					A[i*n + j] -= A[k*n + j] * r; 
				}				
			}
			/*
			for(j = ib + b; j<n ; j++){
				double register sum = 0; 
				for(k = ib; k<i; k++){
					sum += A[i*n + k] * A[k*n +j];
				}
				A[i*n + j] -= sum;
			}
			*/
		}

		/// Update A(end+1:n, end+1:n) 
		/// MM of green = green - blue * pink 
		
		for( i = ib+b; i<n; i+= b){
			for(j = ib+b; j<n; j+=b) mydgemm(A, A, A, n, i, j, ib, b);  
		}

	}	
    return 0;
}

