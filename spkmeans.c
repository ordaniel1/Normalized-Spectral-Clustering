#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*an assertion - in a case of error (in a memory allocation, for example),
 * print an error message and exit*/
#define  MY_ASSERT(x) if (!x) { printf("An Error Has Occured"); exit(0); }


enum goal{spk, wam, ddg, lnorm, jacobi}; /* the options for the goal argument*/

struct map {
    const char *str;
    enum goal value;
} mapEnum[] = { /*match input goal (string) and it's value (enum)*/
        {"spk", spk },{ "wam", wam },
        {"ddg", ddg },{"lnorm", lnorm },
        {"jacobi", jacobi}
};


struct Indices{ /*a struct that will be used in different ways in the program*/
    int i, j;
    double off;
};
typedef struct Indices Indices;


struct Eigen{ /*a struct that will contain an eigenValue and it's position in
                a diagonal matrix*/
    int index;
    double eigenValue;
};
typedef struct Eigen Eigen;



void transpose(double **P, double **Pt, int r, int c);
void printMatrix(double **a, int n, int k);
void printArray(double *a, int n);
double** create2DArray(int numOfRows, int numOfCols);
Indices MaxOffDiagonal(double **a, int n);
double computeTheata(double **a, int i, int j);
double computeT(double theata);
void updateA(double **a, double *irow, double *jrow ,int n,
             int i, int j, double c, double s);

int Jacobi(double **a, double **Vt ,int n);
void initialV(double **v, int n);
void updateV(double **v, double c, double s, int i, int j ,int n);
Indices rotateA(double **a, int n, Indices maxIndices, double c, double s);
int cmpfunc(const void *a, const void *b);
int findK(double **a, double**Vt, int n, int k);
void calcW(double **a, double **w ,int n, int d);
double GaussianRBF(double *v, double *u, int d);
double sumVector(double *v, int n);
void calcD(double **w, double **d, int n, int inputGoal);
void calcLnorm(double *d, double**w, int n);
Indices readData(const char *filename, double **pointsInFile);
double** runSPK(int n, int d, int k, int inputGoal,
                double** points, int python);
double** runJacobiAndKmeans(double** lnorm, int n, int k, int python);
void printD(double *d, int n);
void createT(double **Vt, double **T, int n, int k);
double calacEuclideanDistance(double *b, double *c, int d);
double** kmeans(double **T,int n,int k, int d);
int kmeanspp(double *points, double *cluster, int n, int d, int k);
void kmeansIterations(double** dataPoints, double** centroids,
                      double **oldCentroids, int n, int k, int d);

void copyData(double **target, double **source, int numOfRows, int numOfcols);
void addEigenValues(double **a, double **Vt, int n);
void printJacobi(double **Vt, int n);



int main(int argc, char *argv[]) {
    int n,d, k, i; /* n- number of points in the Data points,
    d - dimension of points in the Data points, k- the input k,
    i- loop's counter variable*/

    char *arg; /*will contain the string of the input goal and
                                   then the name of input file*/

    int inputGoal; /* will be the value of the input goal*/
    double **pointsInFile,**points,**output;/*contain data points and results*/
    Indices nd; /*a struct s.t nd.i=n and nd.j=d after reading the input file*/

    if(argc!=4){ /* the user must pass 3 arguments*/
        printf("Invalid Input!");
        exit(0);
    }

    arg=argv[1];
    k=atoi(arg); /*convert argv[1] to an integer.
                  we can assume that the argument is valid*/

    arg=argv[2];
    inputGoal=0; /*initialize the value of inputGoal*/

    /*find the real value of the input goal by string comparison*/
    for(i=0; i<5; i++){ /* there are 5 options*/
        if (strcmp(arg, mapEnum[i].str)==0){
            inputGoal=i;
            break;
        }
    }

    pointsInFile=create2DArray(50,50);
    nd = readData(argv[3], pointsInFile);
    n = nd.i;
    d = nd.j;

    points= create2DArray(n,d);
    copyData(points, pointsInFile, n, d);/*hold the data in a smaller array*/
    free(*pointsInFile); /*free memory*/
    free(pointsInFile);

    /*calculate the output according to the inputGoal*/
    output=runSPK(n,d,k,inputGoal,points, 0);

    /*print the output according to the inputGoal*/
    switch(inputGoal){
        case 0:break; /*spk*/
        case 1: printMatrix(output,n,n); break; /*wam*/
        case 2: printD(output[0],n); break;  /*ddg*/
        case 3: printMatrix(output,n,n); break; /*lnorm*/
        case 4: printJacobi(output,n); break;  /*jacobi*/
    }

    /*free the memory and return*/
    free(*output);
    free(output);
    return 0;
}


/* double** runSPK(int n, int d, int k, int inputGoal, double** points,
                int python)
 *
 * input:
 * int n-number of points in the data points,
 * int d- dimension of points in the data points
 * int k- the argument 'k' that passed by the user
 * int inputGoal- the value of the goal that passed by the user
 * points- 2D-array (nxd) with the data points
 *
 * output:
 * 2D-array (pointer) - the the desired result of the program,
 * according to the goal.
 */
double** runSPK(int n, int d, int k, int inputGoal, double** points,
                int python) {
    double **ddm;  /*variables for memories that will be allocated*/
    double **w, **Vt;
    /*int i; loop's counter variable*/
    if (inputGoal == 4) { /* if gaol is jacobi*/
        Vt= create2DArray(n+1,n);/*will contain the results*/
        Jacobi(points, Vt, n);/*preform jacobi's algorithm*/
        addEigenValues(points, Vt, n);/*add the eigen values in a's diagonal
                                        to the n'th row of Vt*/
        free(*points); /*free memory*/
        free(points);
        return Vt;
    }
    /* if goal!=4 we have to calculate the Weighted Adjacency Matrix*/
    w= create2DArray(n,n);
    calcW(points, w, n, d);
    free(*points); /* we don't need the data points anymore - free!*/
    free(points);
    if (inputGoal != 1) { /*if goal!=wam*/
        ddm = create2DArray(2,n);
        calcD(w, ddm, n, inputGoal); /*calculate the Diagonal Degree Matrix*/
        if (inputGoal == 2) { /*if goal=ddg*/
            free(*w);
            free(w);
            return ddm;
        }
        else { /*goal is 'lnorm' or 'spk', so calc Normalized Graph Laplacian*/
            calcLnorm(ddm[1], w, n); /*now w is Normalized Graph Laplacian*/
            free(*ddm);
            free(ddm);
        }
    }
    if (inputGoal == 1 || inputGoal == 3) { /*goal is 'wam' or 'lnorm'*/
        return w;
    }

    /* if we get this line - goal=spk*/
    /*return spk Kmeans result (or matrix T of eigenVectors if python==1)*/
    return runJacobiAndKmeans(w,n,k,python);/*w is Normalized Graph Laplacian*/

}




/*double** runJacobiAndKmeans(double** lnorm, int n, int k, int python)
 *
 * input:
 * double** lnorm- Normalized Graph Laplacian (2D-array with size (nxn))
 * int n- size of lnorm
 * int k- the input k that was passed by the user
 *
 * this method computes jacobi algorithm on Normalized Graph Laplacian
 * and kmeans algorithm if python=0, and returns the k centroids.
 * (otherwise - the method computes and returns T matrix of first k
 * eigenVectors that calculated by Jacobi's algorithm).
 *
 * (if k=0 the method computes a new k according to the Eigengap Heuristic)
 *
 */
double** runJacobiAndKmeans(double** lnorm, int n, int k, int python){
    int i; /*loop's counter variable*/

    /* variables for memories that will be allocated*/
    double **originVt, **Vt, **T, *Tnums, **kmeansResult;

    Vt= create2DArray(n,n);/*will hold the eigenvectors*/
    Jacobi(lnorm, Vt, n); /*preform jacobi's algorithm */


    k = findK(lnorm, Vt, n, k);  /*sort eigenvectors and perform
    The Eigengap Heuristic if k=0. now lnorm[0] is pointing to the contiguous
    block of Vt's values*/

    if(k==0||k==n){
        printf("Invalid Input!");
        exit(0);
    }

    originVt= calloc(1,sizeof (double*));
    MY_ASSERT(originVt);/*if calloc failed- print an error message and exit*/

    /*save a pointer to the contiguous block of Vt's values*/
    originVt[0]=lnorm[0];

    lnorm[0]=NULL;
    free(*lnorm);
    free(lnorm);
    if (python==1) { /*if the function was called from Python-C API*/
        Tnums= calloc((n*k)+1, sizeof (double ));
        MY_ASSERT(Tnums);/*if calloc failed- print an error message and exit*/
        T= calloc(n, sizeof (double *)); /*creat 2D-array T*/
        MY_ASSERT(T); /*if calloc failed- print an error message and exit*/
        for(i=0;i<n+1;i++) {
            T[i]=Tnums+k*i;
        }
        T[n][0]=k; /*save the value of k*/
        createT(Vt, T, n, k);
        free(*originVt);
        free(originVt);
        free(Vt);
        return T; /*T contains the first k eigenVectors int its column*/
    }
    T= create2DArray(n,k);
    createT(Vt, T, n, k); /*T contains the first k eigenVectors int its column*/
    free(*originVt);
    free(originVt);
    free(Vt);
    /*compute k-means, dimension of data points is k*/
    kmeansResult=kmeans(T, n, k, k);
    free(*T);
    free(T);
    return kmeansResult;
}




/*Indices readData(const char *filename, double **pointsInFile)
 *
 * input: char* file name, 2D-Array double** pointsInFile (50x50)
 *
 * this method reads data from the file and save it in pointsInFile.
 *
 * output: Indices struct s.t Indices.i=n (number of points),
 * Indices.j=d (dimension of points)
 *
 */
Indices readData(const char *filename, double **pointsInFile){
    FILE* file;
    int d, n, j;
    double value /* will contain a real value that was read from file*/;
    char ch; /* will contain a char that was read from file like: ',' or '\n'*/
    Indices result; /*a struct that will be returned and will contain the
                      values of n-number of points in file, and d - the
                      dimension of points in file.*/

    d=0; /*initialize d - the dimension of points in file*/
    n=0; /*initialize n - number of points in file*/
    j=0; /*a loop's counter variable to determine d*/

    file = fopen(filename, "r");
    MY_ASSERT(file); /*if opening failed - print an error message and exit*/

    while (fscanf(file,"%lf%c", &value, &ch)!=EOF){
        pointsInFile[n][j]=value; /*save value in pointsInFile 2D-array*/
        j++; /*increment the dimension*/
        if (ch=='\n'){ /*if we get end of line in file*/
            d=j; /*dimension is j*/
            j=0; /*j=0- now we will read a new line, or finish the reading*/
            n++; /*increment the number of points in file by 1*/

        }
    }
    fclose(file); /*close the file*/

    /*save n+1 and d in result's fields*/
    result.i=n+1;
    result.j=d;
    return  result;
}




/*void addEigenValues(double **a, double **Vt, int n)
 *
 * input:
 * double **a - a diagonal matrix (2D-array (nxn))
 * double **Vt - a matrix that contains eigenVectors (2D-array ((nx1)xn))
 * int n- size of 'a' and 'Vt'.
 *
 * this method save in the n'th row of Vt the diagonal elements of a
 * (the eigenvalues)
 *
 */
void addEigenValues(double **a, double **Vt, int n){
    int i; /*a loop's counter varaibale*/
    for(i=0; i<n; i++){
       Vt[n][i]=a[i][i]; /*save a[i][i] in Vt[n][i]*/
    }
}



/*void printJacobi(double **Vt, int n)
 *
 * input:
 * double **Vt - a matrix that contains eigenvectors in the it's rows
 * and eigenvalues in it's n-th row (Vt is a 2D-array ((nx1)xn)).
 * int n - number of columns in Vt (n+1 is the number of rows).
 *
 * this method prints the eigenvalues in one row, and then prints
 * the eigenvectors - an eigenvector in a row.
 *
 */
void printJacobi(double **Vt, int n){
    printArray(Vt[n], n); /*print the last row of Vt - the eigenvalues*/
    printf("\n");/*skip a row*/
    printMatrix(Vt,n,n); /*print the eigenvectors*/

}



/* double** kmeans(double **T,int n,int k, int d)
 *
 * input:
 * double **T - data points (2D-array (nxd))
 * int n - number of points in the data points - number of "rows" in T
 * int k - number of desired clusters
 *
 * this method perform the k-means algorithm on T - the data points,
 * and returns a 2D-array that contains the K centroids of the k clusters.
 */
double** kmeans(double **T,int n,int k, int d) {

    /*initialize*/
    double **centroids, **oldCentroids; /*will contain the current and previous
                            centroids during the iterations of the algorithm*/
    int i, j, counter; /*loop's counter variables */

    centroids= create2DArray(k,d+1); /* allocate 2D-array (k, d+1):
    k "rows" for k centroids, d+1 "columns" s.t the first d columns will
    contain the coordinates of the centroid, and the last column will
    contain the number of points that belong to the cluster */

    /*initialize the centroids as the first k points in T (the data points).
    this loop saves the coordinates of each centroid in the first d 'columns'
    of 'centroids' (each centroid in a different row)*/
    counter = 0;
    for (i = 0; i < k * (d + 1); i++) {
        for (j = 0; j < d; j++) {
            (*centroids)[i] = (*T)[i - counter];
            i++;
        }
        counter++;
    }

    oldCentroids= create2DArray(k, d+1);/*allocate 2D-array (k, d+1):
    k "rows" for k centroids, d+1 "columns" s.t the first d columns will
    contain the coordinates of the previous centroids, and the last column will
    contain the number of points that belong to the cluster*/


    /*perform the k-means algorithm*/
    kmeansIterations(T,centroids,oldCentroids,n,k,d);


    free(*oldCentroids); /*free the previous centroids*/
    free(oldCentroids);
    printMatrix(centroids, k, k); /*print the centroids*/
    return centroids; /*return the centroids that will
                      be freed immediately in the main function*/
}

/*int kmeanspp(double *points, double *cluster, int n, int d, int k)
 *
 *input:
 * double *points - 1-D contiguous array (size: n*d) that contains the
 *                  data points
 * double *cluster - 1-D contiguous array (size: k*(d+1)) that contains the
 *                  K initial centroids.
 * int n - number of points in the data points.
 * int d - the dimension of the points in the data points
 * int k - the number of desired clusters
 *
 * this method perform the kmeans algorithm s.t after the calculations,
 * 'clusters' contains (contiguously) the k centroids that the algorithm found.
 *
 * output: 0
 */
int kmeanspp(double *points, double *cluster, int n, int d, int k){

    /*initialize*/
    int i; /*loop's counter variables */
    double **matrixPoints; /*array of pointers that will point to the
                            contiguous array 'points' in order to work
                             with the data points as a 2D-Array*/

    double **centroids; /*array of pointers that will point to the
                        contiguous array 'cluster' in order to work
                        with the data points as a 2D-Array*/

    double **oldCentroids; /*will contain the previous centroids
                             during the iterations of the algorithm*/


    matrixPoints=calloc(n, sizeof(double*));
    MY_ASSERT(matrixPoints); /*if calloc failed- print an error message
                                and exit*/

    /*initial the pointers in matrixPoints s.t they will point to the right
     positions in 'points' in order to work with the data points as
     a 2D-Array*/
    for(i=0; i<n; i++){
        matrixPoints[i]=points+i*d;
    }

    centroids=calloc(k, sizeof(double*));
    MY_ASSERT(centroids); /*if calloc failed- print an error message
                                and exit*/

    /*initialize the pointers in centroids s.t they will point to the right
     positions in 'cluster' in order to work with the centroids as
     a 2D-Array*/
    for(i=0; i<k; i++){
        centroids[i]=cluster+(i*(d+1));
    }

    /* allocate 2D-array (k, d+1):
    k "rows" for k centroids, d+1 "columns" s.t the first d columns will
    contain the coordinates of the previous centroids, and the last column will
    contain the number of points that belong to the cluster*/
    oldCentroids=create2DArray(k, d+1);

    /*perform the k-means algorithm*/
    kmeansIterations(matrixPoints, centroids, oldCentroids, n, k, d);

    /*free all the memory that we allocated*/
    free(matrixPoints);
    free(*oldCentroids);
    free(oldCentroids);
    free(centroids);
    return 0;
}




/*void kmeansIterations(double** dataPoints, double** centroids,
                      double **oldCentroids, int n, int k, int d)
 *
 * input:
 * double** dataPoints - 2D-array (nxd) of data points
 * double** centroids - 2D-array (nx(d+1)) of k centroids
 * double **oldCentroids - 2D-array (nx(d+1))
 * int n - number of points in the data points
 * int k - number of desired clusters
 * int d - the dimension of points in the data points
 *
 * this method perform the K-means algorithm.
 * after the calculations, 'centroids' will contain the K centroids
 * that the algorithm found.
 *
 */
void kmeansIterations(double** dataPoints, double** centroids,
                      double **oldCentroids, int n, int k, int d){
    int i, j, m; /*loop's variables*/
    int clustersIndex, equalCentroids; /*clusteIndex- will contain
    index of a cluster during the calculation, and equalCentroids - will
    count the number of centroids that equal to check if the algorithm
    get the stop condition*/

    double distance, minDistance, epsilon;
    /*distance - to check the distance between point and centroids
                 and between centroids
    minDistance - save the minimal distance between point and a centroid,
                    in order to determine the cluster for the point*/

    epsilon=0.0000000001; /*to check equality (epsilon approximation)*/
    for(i=0; i<300;i++){/*this is the main loop of the algorithm MAX_ITER=300*/
        for(m=0;m<k*(d+1); m++){/*save the current centroids in 'oldCentroids*/
            (*oldCentroids)[m]=(*centroids)[m];
            (*centroids)[m]=0;/*now each centroid in 'centroids' is zero*/
        }
        for(j=0; j<n; j++){ /*for each point i dataPoint*/
            /*initialize minDistance as the distance
                between point and the first centroid*/
            minDistance=calacEuclideanDistance(dataPoints[j],
                                                    oldCentroids[0],d);
            clustersIndex=0; /*clustersIndex is the index of first centroid*/
            for(m=1; m<k; m++){ /*for each centroid (1 to k-1)*/
                /*calculate the distance between point and centroid*/
                distance=calacEuclideanDistance(dataPoints[j],
                                                    oldCentroids[m],d);
                if(distance<minDistance){ /*if th distance is smaller*/
                    minDistance=distance; /*update minDistance*/
                    clustersIndex=m; /*update clustersIndex*/
                }
            }
            /*we found the cluster for the point*/
            for(m=0; m<d; m++){ /*sum the points' coordinates in the cluster*/
                centroids[clustersIndex][m]+=dataPoints[j][m];
            }
            /*increment the number of points in the cluster by 1*/
            centroids[clustersIndex][d]+=1;
        }

        equalCentroids=0; /*initialize to be 0'*/
        /*compute the new centroids by divide the sum of each coordinate by
        the number of points in the cluster*/
        for(m=0; m<k; m++){
            for(j=0; j<d; j++){
                centroids[m][j]=(centroids[m][j]/centroids[m][d]);
                }


            /*the last column of each centroid (number of points) is zero*/
            centroids[m][d]=0;
        }

        /*check if all of the new centroids are equal to the older centroids*/
        for(m=0; m<k; m++){
            distance=calacEuclideanDistance(oldCentroids[m], centroids[m],d);
            if(distance<epsilon){ /*if the centroids are equal*/
                equalCentroids++; /*increment the numer of equalCentroids*/
            }
        }
        if(equalCentroids==k){ /*the break condition*/
            break;/*all of the new centroids are equal to the older centroids*/
        }
    }
}




/*double calacDistance(double *b, double *c, int d)
 * input: 2 d-size arrays of doubles (vectors) and int d (size of arrays)
 * output: euclidean squared distance between the two input vectors
 */
double calacEuclideanDistance(double *b, double *c, int d){
    int i; /*loop's variable*/
    double res=0; /*initialize the result value to be zero*/
    for(i=0; i<d; i++){ /* for each coordinate*/
        /*calcaulate the (b[i]-c[i])^2 and add to res*/
        res+=(*(b+i)-*(c+i))*(*(b+i)-*(c+i));
    }
    return res;
}




/*
 *input: 2 d-size arrays of doubles (vectors) and int d (size of arrays)
 *output: distance between the 2 vectors according to GaussianRBF
 */
double GaussianRBF(double *v, double *u, int d){
    double euclideanDistance;
    euclideanDistance=calacEuclideanDistance(v,u,d);
    return exp((-1*sqrt(euclideanDistance))/2);
}




/*void calcW(double **a, double **w ,int n, int d)
 *
 *input:
 * double **a- 2D (nxn) array - the data points,
 * double **w- 2D (nxn) array - w,
 * int n- number of points (number 'rows' in a),
 * int d - dimension of points (number of 'columns' in a)
 *
 * this method computes W - the Weighted Adjacency Matrix:
 * w will contain the values of the Weighted Adjacency Matrix.
 */
void calcW(double **a, double **w ,int n, int d){
    int i,j; /*loop's variables*/
    for (i=0; i<n; i++){ /*for every point p in a*/
        for (j=0; j<n; j++){ /*for every point q!=p in a*/
            if (i!=j){ /*calculate the distance between p and q
                                            using GaussianRBF*/
                w[i][j]=GaussianRBF(a[i], a[j], d);
            }
            else{ /*if points are equal (i=j)*/
                w[i][i]=0;
            }
        }
    }
}




/*void calcD(double **w, double **d, int n, int inputGoal)
 *
 * input:
 * double **w - Weighted Adjacency Matrix,
 * double **d - 2D array d (2xn),
 * int n - the size of the Weighted Adjacency Matrix,
 * int inputGoal - the value of the goal that passed by the user
 *
 * This method computes  D- the Diagonal Degree Matrix,
 * using the Weighted Adjacency Matrix.
 * after the calculation, d contains in it's first 'row'
 * the diagonal elements of D, and if goal!=ddg,
 * it's second 'row' contains the diagonal element of D^(-0.5).
 */
void calcD(double **w, double **d, int n, int inputGoal){
    int i;
    double x;

    /*there are two options here, in order to avoid calculation when goal=ddg*/
    if (inputGoal==2){
        for (i=0; i<n; i++) {
            x = sumVector(w[i], n);/*sum of elements in i's row in w*/
            d[0][i] = x; /*the value of D[i][i] is saved in d[0][i] */
        }
    }
    else{
        for (i=0; i<n; i++){
            x=sumVector(w[i],n);/*sum of elements in i's row in w*/
            d[0][i]=x; /*the value of D[i][i] is saved in d[0][i] */
            d[1][i]= 1/sqrt(x);/*the value of D^(-0.5)[i][i] is saved
                                                           in d[1][i]*/
        }
    }

}




/*void printD(double *d, int n)
 *
 * input:
 * double *d - 1-D array
 * int n - size of the array
 *
 * this method prints the Diagonal Degree Matrix.
 *
 */
void printD(double *d, int n) {
    int i, j; /*loop's variables*/
    for (i = 0; i < n - 1; i++) { /*for each 'row'*/
        for (j = 0; j < i; j++) { /*print the off diagonal elements - zeros*/
            printf("0.0000,");
        }
        /*print D[i][i]*/
        /*if D[i][i] is a very small negative value - print as zero*/
        if (d[i] < 0 && d[i] > -0.00005) {
            printf("0.0000,");
        } else {
            printf("%.4f,", d[i]);
        }

        for (j=i+1; j<n-1; j++) { /*print the off diagonal elements - zeros*/
            printf("0.0000,");
        }
        printf("0.0000\n");
    }

    /*print last row - same as before, without skip a row in the end*/
    for (j = 0; j < n - 1; j++) {
        printf("0.0000,");
    }
    if (d[n - 1] < 0 && d[n - 1] > -0.00005) {
        printf("0.0000,");
    } else {
        printf("%.4f,", d[n - 1]);
    }
}




/* double sumVector(double *v, int n)
 * input:
 * double *v - 1-D array v and int n - number of elements in v
 * output: sum of elements in v
 */
double sumVector(double *v, int n){
    int i; /*loop's variable*/
    double sum; /*the result's variable*/
    sum=0;
    for (i=0; i<n; i++){ /*add v[i] to sum*/
        sum+=v[i];
    }
    return sum;
}




/*void calcLnorm(double *d, double**w, int n)
 *
 * input:
 * double *d - an array (size:n) that contains the diagonal elements
 *              of D^(-0.5), while D is Diagonal Degree Matrix.
 * double **w - a 2D-array (size: nxn) - Weighted Adjacency Matrix
 * int n - size of d and w.
 *
 * after the calculations, w contains the values of Normalized Graph Laplacian
 */
void calcLnorm(double *d, double**w, int n){
    /* calac C=W(D^-0.5) */
    int i, j;
    for (i=0; i<n; i++){
        for(j=0; j<n; j++){
            w[i][j]=w[i][j]*d[j];
        }
    }
    /*calc I-D^(-0.5)C */
    for (i=0; i<n; i++){
        for(j=0; j<n; j++){
            w[i][j]=(-1)*w[i][j]*d[i];
        }
        w[i][i]=1;
    }



}

/*Indices MaxOffDiagonal(double **a, int n)
 *
 * input:
 * double **a - a symmetric matrix (2D array (nxn))
 * int n - size of a.
 *
 * output: an Indices struct s.t:
 * a[struct.i][struct.j]=off-diagonal element of 'a' with
 *                       the largest absolute value.
 * struct.off= off(a)^2= the sum of squares of all off-diagonal elements of 'a'
 *
 */
Indices MaxOffDiagonal(double **a, int n){
    int i,j; /*loop's variables*/
    double sum; /*will contain off(a)^2*/
    Indices indices; /*the struct that will be returned*/
    indices.i=0; /*start from the a[0][1]*/
    indices.j=1;
    sum=0; /*initialize sum to be 0*/
    for(i=0; i<n; i++){ /*for each row*/
        /*sum+=a[i][i]*a[i][i];*/
        for(j=i+1; j<n; j++){ /*for each position above the diagonal
                                (a is symmetric matrix so we don't need to
                                look at elements under the diagonal)*/

            sum+=2*a[i][j]*a[i][j]; /*add the square of the element to sum
 *                                   twice - because it appears twice in a*/

            /*if the absolute value of a[i][j] is bigger - update indices */
            if (fabs(a[i][j])>fabs(a[indices.i][indices.j])){
                indices.i=i;
                indices.j=j;
            }
        }
    }
    indices.off= sum; /*ndiced.off=off(a)^2*/

    return indices;

}




/*double computeTheata(double **a, int i, int j)
 *
 *input:
 *double **a - a symmetric matrix (2D array)
 * int i, int j s.t a[i][j]= off-diagonal element of 'a' with
 * the largest absolute value,  and a[i][j]!=0.
 *
 * output: value of theata parameter for jacobi's algorithm.
 *
 */
double computeTheata(double **a, int i, int j){
    return ((a[j][j]-a[i][i])/(2*a[i][j]));
}




/*double computeT(double theata)
 *
 * input: double theata - a parameter from Jacobi's algorithm.
 * output: value of t parameter in for Jacobi's algorithm.
 *
 */
double computeT(double theata){
    double x;
    if (theata<0){
        x=-1;
    }
    else{
        x=1;
    }
    /*calculations according to the formula*/
    return x/(fabs(theata)+(sqrt(pow(theata,2)+1)));
}




/* double** create2DArray(int numOfRows, int numOfCols)
 *
 * this method creates contiguous 2D-Array
 *
 * input:
 * int numOfRows - number of 'rows' that are desired in the array
 * int numOfCols - number of 'columns' that are desired in the array
 *
 * output: a pointer to array with numOfRows pointers that point to
 * positions in contiguous array of doubles (size: numOfRows*numOfCols)
 *
 */
double** create2DArray(int numOfRows, int numOfCols){
    double *elements; /* contiguous array of (numOfRows*numOfCols) doubles*/
    double **matrix; /* contiguous array of (numOfRows) pointers to double*/
    int i; /*loop's variable*/
    elements=calloc(numOfRows*numOfCols, sizeof (double ));
    MY_ASSERT(elements); /*if calloc failed print an error message and exit*/
    matrix=calloc(numOfRows, sizeof (double*));
    MY_ASSERT(matrix); /*if calloc failed print an error message and exit*/
    for(i=0; i<numOfRows; i++){ /*create the pointers according to the number
                                                        of rows and columns*/
        matrix[i]=elements+i*numOfCols;
    }
    return matrix; /*return a pointer*/
}




/* void transpose(double **P, double **Pt, int r, int c)
 *
 * input:
 * double **P - 2D array (size : rxc)
 * double **Pt - 2D array (size : cxr)
 * int r - number of 'rows' in P
 * int c - number of 'columns' in P
 *
 * after the calculation Pt is the transpose of P
 */
void transpose(double **P, double **Pt, int r, int c) {
    int i, j; /*loop's variables*/
    for (i = 0; i < r; i++) { /*for ech row in P*/
        for (j = 0; j < c; j++) { /*for each column in c*/
            Pt[j][i] = P[i][j]; /*transpose !*/
        }
    }
}




/*void printArray(double *a, int n){
 *
 * input:
 * double *a - 1D-Array a, int n - number of elements in a.
 * This method prints a's elements in one row, seperated by ','
 *
 */
void printArray(double *a, int n){
    int i; /*loop's variable*/
    for(i=0; i<n-1; i++){ /*for each of the first n-1 element*/
        if(a[i]<0&&a[i]>-0.00005){ /*if a[i] is a very small negative number*/
            printf("0.0000,"); /*print it as zero*/
        }
        else{
            printf("%.4f,", a[i]);
        }

    }

    /*print the last element - same as before, but without ',' in the end*/
    if(a[n-1]<0&&a[n-1]>-0.00005){
        printf("0.0000");
    }
    else{
        printf("%.4f", a[n-1]);
    }

}




/* void printMatrix(double **a, int n, int k)
 *
 * input:
 * double **a- 2D-Array
 * int n - number of 'rows' in a
 * int k - number of 'columns' int a
 * This method prints 'a' as a 2D-Matrix, row by row, seperated by '\n'.
 *
 */
void printMatrix(double **a, int n, int k) {
    int i; /*loop's variable*/

    for(i=0; i<n-1; i++){ /*for each of the first n-1 rows*/
        printArray(a[i],k); /*print the i'th row*/
        printf("\n"); /*skip a roe*/
    }
    printArray(a[n-1],k); /*print the last row without skip a row in the end*/

}









/*void updateV(double **v, double c, double s, int i, int j ,int n)
 *
 * input:
 * double **v - 2D-Array (nxn)
 * double c, double s - elements in matrix P from Jacobi's algorithm.
 * int i, int j - position in matrix P from Jacobi's algorithm, s.t:
 * P[i][i]=c, P[i][j]=s, P[j][i]=-s, P[j][j]=c
 * (we don't need to hold p explicitly, because all of the other diagonal
 * elements are 1, and all of the other off-diagonal elements of P equal to 0).
 *
 * after this calculation, v is equal to v*p
 */
void updateV(double **v, double c, double s, int i, int j ,int n){
    int k;
    double temp; /*will contain a value of element in v that will be update
                during the calculation but we'll still need it for other
                calculations*/

    /*update only the c elements of v*/
    for(k=0; k<n; k++){ /*for each row in v */
        temp=v[k][i]; /*save v[k][i] in temp*/
        /*update relevant elements in v*/
        v[k][i]=c*temp-s*v[k][j];
        v[k][j]=s*temp+c*v[k][j];
    }
}




/*void updateA(double **a, double *irow, double *jrow,
             int n, int i, int j, double c, double s)
 *
 * input:
 * double **a - symmetric matrix (2D-array (nxn))
 * double *irow - the i-th 'row' of a (1D array)
 * double *jrow - the j-th 'row' of a (1D array)
 * int n - size of a's 'rows' and a's 'columns'.
 * int i, int j - positions in 'a' s.t a[i][j]=
 * the off-diagonal element with the largest absolute value.
 *
 * double c, double s - element in P.
 *
 * after the calculations,  A'=(Pt)AP (a will contains the values of A').
 * (we don't need to hold p explicitly, because all of the other diagonal
 * elements are 1, and all of the other off-diagonal elements of P equal to 0).
 *
 */
void updateA(double **a, double *irow, double *jrow,
             int n, int i, int j, double c, double s){
    int r; /*loop's variable*/
    /*update only the relevant elements and avoid matrix multiplication*/
    for(r=0; r<n; r++){
        a[r][i]=c*irow[r]-s*jrow[r];
        a[i][r]=c*irow[r]-s*jrow[r];
        a[r][j]=c*jrow[r]+s*irow[r];
        a[j][r]=c*jrow[r]+s*irow[r];
    }
    a[i][i]=(c*c)*irow[i]+(s*s)*jrow[j]-2*s*c*irow[j];
    a[j][j]=(s*s)*irow[i]+(c*c)*jrow[j]+2*s*c*irow[j];
    a[i][j]=0;
    a[j][i]=0;
}

/*void initialV(double **v, int n)
 *
 * input:
 * double **v - 2D-array (size: nxn)
 * int n - size of v's 'rows' and 'columns'
 *
 * After the calculation - v=I
 *
 */
void initialV(double **v, int n){
    int i; /*loop's variable*/
    for(i=0; i<n; i++) { /* for each row*/
        v[i][i] = 1;
    }

}




/*Indices rotateA(double **a, int n, Indices maxIndices, double c, double s)
 *
 * input:
 * double **a - double **a - symmetric matrix (2D-array (nxn)),
 * int n - size of a's 'rows' and 'columns',
 * Indices maxIndices s.t.
 * a[maxIndices.i][maxIndices.j]=the off-diagonal element with the
 *                                 largest absolute value,
 * double c and double s - elements in matrix P from Jacobi's algorithm
 * after the calculations,  A'=(Pt)AP (a will contains the values of A').
 *
 * output: an Indices struct s.t:
 * A'[struct.i][struct.j]=off-diagonal element of 'A'' with
 *                       the largest absolute value.
 * struct.off= off(A')^2= the sum of squares of all off-diagonal
 * elements of 'A''
 */
Indices rotateA(double **a, int n, Indices maxIndices, double c, double s){
    int i; /*loop's variable */
    double *irow, *jrow; /*will contain the i an j rows of a */
    Indices maxIndices2; /*the output struct*/
    irow = calloc(n, sizeof(double));
    MY_ASSERT(irow);/*if calloc failed, print an error message and exit*/
    jrow = calloc(n, sizeof(double));
    MY_ASSERT(jrow);/*if calloc failed, print an error message and exit*/
    for (i = 0; i < n; i++) {
        irow[i] = a[maxIndices.i][i]; /*save the i-th row of a in irow*/
        jrow[i] = a[maxIndices.j][i]; /*save the j-th row of a in jrow*/
    }
    /*update a to be a' - A'=(Pt)Ap*/
    updateA(a, irow, jrow, n, maxIndices.i, maxIndices.j, c, s);

    /*calc fields of maxIndices2 that will contain the position of the
      off-diagonal element of A' with the largest absolute value,
      and off(A')^2*/
    maxIndices2= MaxOffDiagonal(a,n);

    /*free memory that we allocated*/
    free(irow);
    free(jrow);

    return maxIndices2;

}


/*int Jacobi(double **a, double **Vt ,int n)
 * Method that computes Jacobi's Algorithm
 * input: a symetric matrix 'a' (2D-array) with size of (nxn),
 * a 2D-Array 'Vt' with size of (nxn) and
 * an int 'n' - the size of 'a'.
 * after the calculations, a is a diagonal matrix, it's diagonal contains it's
 * eigenvalues and Vt's rows contain the eigenvectors.
 *
 * output: 0;
 */

int Jacobi(double **a, double **Vt ,int n) {
    int iter;/* counter of jacobi's iterations*/
    double **v; /*matrix v that it's columns will be the eigenvectors*/
    double epsilon, theata, t, c, s; /*parameters for the computations */

    Indices maxIndices, maxIndices2; /* structs that will contain the position
                                    (i,j fields) of the off-diagonal element
                                    with the largest absolute value and the
                                    square of off diagonal elements (off field)
                                    of a and a'*/

    epsilon = pow(10, -15);
    v= create2DArray(n,n); /*creates a 2D-array*/
    initialV(v, n); /*v[i][i]=1 for 0<=i<=n-1 */
    maxIndices = MaxOffDiagonal(a, n); /*find the position of the off-diagonal
                                    element with the largest absolute value and
                                    the square of off diagonal elements of a*/

    if(a[maxIndices.i][maxIndices.j]==0){ /*if 'a' is already diagonal- stop*/
        return 0;
    }
    /* compute s and c according to the algorithm*/
    theata = computeTheata(a, maxIndices.i, maxIndices.j);
    t = computeT(theata);
    c = 1 / (sqrt((t * t + 1)));
    s = t * c;
    /*update v (v=p1)*/
    v[maxIndices.i][maxIndices.i] = c;
    v[maxIndices.i][maxIndices.j] = s;
    v[maxIndices.j][maxIndices.i] = -s;
    v[maxIndices.j][maxIndices.j] = c;
    for (iter = 0; iter < 100; iter++) {

        /*Update a to be a': a'=(Pt)aP and find the the position of the
         off-diagonal element with the largest absolute value and the
         square of off diagonal elements of a' */
        maxIndices2=rotateA(a, n, maxIndices, c, s);

        if (((maxIndices.off-maxIndices2.off)<=epsilon)||
            (a[maxIndices2.i][maxIndices2.j]==0)){
            break; /*if we get the stop condition
                        or a is already diagonal, then break*/
        }
        /*here we didn't get the stop condition, therefore we update V */
        maxIndices=maxIndices2;
        theata = computeTheata(a, maxIndices.i, maxIndices.j);
        t = computeT(theata);
        c = 1 / (sqrt((t * t + 1)));
        s = t * c;
        updateV(v, c, s, maxIndices.i, maxIndices.j, n);
    }
    transpose(v,Vt,n, n); /*Vt is the transpose of v*/
    free(*v);
    free(v);
    return 0;

}


/*int cmpfunc(const void *a, const void *b)
 *
 * a comparison function for a stable sorting in qsort
 * input: 2 'Eigen' structs a and b.
 * output: 1 if a's eigenValue field is bigger than b's eigenValue field,
 *        -1 otherwise.
 */
int cmpfunc(const void *a, const void *b){
    Eigen eigenA, eigenB;
    double x;
    /*int temp;*/
    eigenA=*(Eigen*)a;
    eigenB=*(Eigen*)b;
    x=eigenA.eigenValue-eigenB.eigenValue;
    if (x>0){
        return 1;

    }
    if(x<0){
        return -1;
    }
    else{ /*if there is equality of eigen values,
           then, sort by index*/
        if(eigenA.index<eigenB.index){
            return -1;
        }
        if(eigenA.index>eigenB.index){
            return 1;
        }
        else{
            return 0;
        }

    }
}




/*int findK(double **a, double**Vt, int n, int k)
 *
 * input:
 * a diagonal (nxn) matrix 'a' (2D-Array),
 * a (nxn) matrix 'Vt' (2D-Array), that it's rows
 * are eigenVectors, (s.t the i's row belongs to the eigenValue a[i][i]
 * for 0<=i=<n-1)
 * int k.
 *
 * This method sort the rows of Vt (eigenVectors) according to the
 * increasing of the diagonal elements of a.
 *
 * output: if k!=0 - k is returned, and if k==0, this method computes a new k
 * according to The Eigengap Heuristic, and returns it.
 */
int findK(double **a, double**Vt, int n, int k){
    int i; /* a loop's counter*/
    int newk,m; /*if k==0, newk will be computed, m will be floor(n/2)*/
    double delta ,maxdelta; /*deltas of The Eigengap Heuristic calculations*/

    Eigen *eigenValues; /*an array that will contain n 'Eigen' structs -
    each struct has fields of eigenValue and index of the eigenVlaue
    according to the position i on a's diagonal*/

    eigenValues= calloc(n, sizeof(Eigen));
    MY_ASSERT(eigenValues);
    for (i=0;i<n; i++){ /*build the array*/
        Eigen eigen;
        eigen.index=i;
        eigen.eigenValue=a[i][i];
        eigenValues[i]=eigen;
    }


    qsort(eigenValues, n, sizeof(Eigen), cmpfunc); /*sort the array*/

    free(*a); /* free the contiguous memory that contains the values
                of matrix a */
    for(i=0; i<n; i++){
        a[i]=Vt[i]; /* a's pointers are pointing now to the contiguous memory
                     that contains the values of matrix Vt */
    }

    /*update Vt's rows according to the order of the eigenValues*/
    for(i=0; i<n; i++){
        Vt[i]=a[eigenValues[i].index];
    }

    if(k!=0){
        free(eigenValues);
        return k;
    }
    /*compute newK according to the Eigengap Heuristic*/
    m= (int)floor(n/2);
    if(n==1||n==2){
        free(eigenValues);
        return 1; /*only one gap if n==2*/
    }
    /*maxdelta is the first gap untill we will find a bigger gap*/
    maxdelta=fabs(eigenValues[0].eigenValue-eigenValues[1].eigenValue);
    newk=1; /* the index of the first delta*/
    for(i=1; i<=m; i++){
        delta=fabs(eigenValues[i].eigenValue-eigenValues[i+1].eigenValue);
        if (maxdelta<delta){
            maxdelta=delta;
            newk=i+1; /* we count the gaps from 1, not from 0
                            so we increment i by 1*/
        }
    }



    free(eigenValues);
    return newk;
}




/* void normalizeVector(double *v, int n)
 * input:
 * double *v - a non-zero Vector (1D array) v,
 * int n - size of v
 * after the calculations, v is a normalized vector.
 */
void normalizeVector(double *v, int n){
    double sum, squareRootSum;
    int i; /*loop's variables */
    sum=0;
    /*calc the sum of squares of v's elements*/
    for (i=0; i<n; i++){ /*for each element in v*/
        sum+=(v[i])*(v[i]); /*add it's square to sum*/
    }
    squareRootSum= sqrt(sum); /*sqrt of sum of squares of v's elements*/
    for(i=0; i<n; i++){ /* for each element in v*/
        v[i]=v[i]/squareRootSum; /*divide it by squareRootSum*/
    }
}




/* void normalizeMatrix(double **a, int n, int k)
 * input:
 *
 * double **a - Matrix a (2D-array),
 * int n=number of 'rows' in a,
 * int k=number of 'columns' in a.
 * after the calculations - a's rows are normalized
 */
void normalizeMatrix(double **a, int n, int k){
    int i; /*loop's variable*/
    for(i=0; i<n; i++){ /*for ecah row*/
        normalizeVector(a[i], k); /*normalize the row*/
    }
}




/* void createT(double **Vt, double **T, int n, int k)
 *
 * input:
 * double **Vt - 2D (nxn) Array,
 * double **T - 2D (nxk) array T,
 * int n - number of 'rows' in T
 * int k - number of 'columns' in T
 *
 * after the calculations T's rows will be normalization of first K columns of V
 * (that are the first K rows of Vt)
 */
void createT(double **Vt, double **T, int n, int k){
    transpose(Vt, T, k, n); /*we need only the first k rows of Vt - T size will be (nxk)*/
    normalizeMatrix(T, n, k); /*normalize T's rows*/

}




/*void copyData(double **target, double **source, int numOfRows, int numOfcols)
 *
 * input:
 * double **target - 2D array,
 * double **source - 2D array,
 * int numOfRows - number of 'rows' in target
 * which <=number of 'rows' in source,
 *
 * int numOfcols - number of 'columns' in target
 * which <=number of 'columns' in source
 *
 * after the calculation - target contains the values of
 * first numOfRows rows and first numOfcols columns of source.
 *
 */
void copyData(double **target, double **source, int numOfRows, int numOfcols){
    int i,j;
    for(i=0; i<numOfRows; i++){ /*for each row*/
        for(j=0; j<numOfcols; j++){ /*for each column*/
            target[i][j]=source[i][j]; /*save data*/
        }
    }

}













