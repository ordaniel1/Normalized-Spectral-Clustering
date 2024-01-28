#ifndef SPKMEANS_H
#define SPKMEANS_H

/*prototypes of functions in spkmeans.c that are called from spkmeansmodule.c*/
double** runSPK(int n, int d, int k, int inputGoal, double** points, int python);
int kmeanspp(double *points, double *cluster, int n, int d, int k);
double** create2DArray(int numOfRows, int numOfCols);

#endif /*SPKMEANS_H*/

/*an assertion - in a case of error (in a memory allocation, for example),
 * print an error message and exit*/
#define MY_ASSERT(x) if (!x) { printf("An Error Has Occured"); exit(0); }