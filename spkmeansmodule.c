#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <assert.h>
#include "spkmeans.h"
#include <stdlib.h>
#include <stdio.h>


/* static PyObject* spk(PyObject *self, PyObject *args)
 *
 * input (from python):
 * points_list - python list (1-dimnension) of data points
 * int n - number of points in points_list
 * int d- dimension of points in points_list
 * int k- the argument 'k' that passed by the user
 * int inputGoal- the value of the goal that passed by the user
 *
 * output:
 * python list of floats that contains the desired result according
 * to the goal.
 *
 */
static PyObject* spk(PyObject *self, PyObject *args){
    PyObject *points_list, *python_res; /*save python's data-input and output*/
    int n, d, k, inputGoal; /*save python's argument */
    double **points, **res; /*will save c's data*/
    int i,j; /*loop's variables*/

    if (!PyArg_ParseTuple(args, "Oiiii", &points_list,  &n, &d, &k
                                                    ,&inputGoal)){ /*parsing*/
        printf("An Error Has Occured"); /*wrong arguments! */
        exit(0); /*exit the program*/
    }
    /*allocate memory to save data points from 'points_list' as doubles*/
    points=create2DArray(n,d); /* in order to wrok with 2D-array*/
    for (i=0; i<n; i++){
        for(j=0; j<d; j++){
            PyObject *num; /*num will be points_list[i*d+j]:*/
            num=PyList_GetItem(points_list, ((i*d)+j));
            points[i][j]=PyFloat_AsDouble(num); /*convert to double and save*/
        }
    }
    /*calc the desire value according to goal*/
    res=runSPK(n,d,k,inputGoal,points,1);/*points was freed in runSPK*/
    switch (inputGoal) { /*return result to python, according to the goal*/
        case 0: /*goal=spk, then res=T*/
            k=res[n][0]; /*k that was calculated is the last element in T*/
            python_res = PyList_New(n*k); /*python list to save T's values*/
            for(i=0; i<n*k; i++) {
                /*convert element to python float*/
                PyObject *python_float = Py_BuildValue("d", (*res)[i]);
                /*save element in python_res*/
                PyList_SetItem(python_res, i, python_float);
            }
            break;
        case 2: /*goal=ddg, then res=D (2xn)*/
            python_res = PyList_New(n); /*python list to save D[0]'s values*/
            for (i=0; i<n; i++) {/*convert elements to python float and save*/
                PyObject *python_float = Py_BuildValue("d", (*res)[i]);
                PyList_SetItem(python_res, i, python_float);
            }
            break;
        case 4:/*goal=ddg, then res=Vt ((n+1)xn)-eigenvalues and eigenvectors*/
            python_res=PyList_New(n*(n+1)); /*python list to save Vt's values*/
            for (i=0;i<n*(n+1);i++){
                /*convert element to python float and save in python_res*/
                PyObject *python_float = Py_BuildValue("d", (*res)[i]);
                PyList_SetItem(python_res, i, python_float);
            }
            break;
        default: /*goal=wam or goal=lnorm, the res in (nxn) matrix*/
            python_res=PyList_New(n*n); /*python list to save values of res*/
            for (i = 0; i < n*n; i++) {
                /*convert element to python float and save in python_res*/
                PyObject *python_float = Py_BuildValue("d", (*res)[i]);
                PyList_SetItem(python_res, i, python_float);
            }
            break;
    }
    free(*res); /*free memory that we've allocated in runSPK*/
    free(res);
    return python_res; /*return the result to python*/
    }



/* static PyObject* fit(PyObject *self, PyObject *args)
 *
 * this function called from python program, when goal=spk
 * and perform spkmeans++ algorithm.
 *
 * input (from python):
 * points_list - python list (1-dimnension) of data points
 * centroids_list - python list (1-dimnension) of first k
 *      cenrtoids, that were chosen according to kmeans++ algorithm.
 * int n - number of points in points_list
 * int d- dimension of points in points_list
 * int k- the desired number of clusters.
 *
 *
 * output:
 * python list of floats that contains the k centroids that
 * were calculated by kmeans algorithm.
 *
 */
static PyObject* fit(PyObject *self, PyObject *args){
    PyObject *points_list; /*save python's data points (input)*/
    PyObject *centroids_list; /*save python's centroids (input)*/
    int n,d,k, i,j,counter; /*arguments and loop's variables*/
    double *points; /*will save data points as c's data*/
    double *clusters; /*will save centroids as c's data*/
    if (!PyArg_ParseTuple(args, "OOiii", &points_list, &centroids_list,
                                                      &n, &d, &k)){ /*parsing*/
        printf("An Error Has Occured"); /*wrong arguments! */
        exit(0); /*exit the program*/
    }
    /*save the data points (points_list values) in 1D array of doubles*/
    points=calloc(n*d, sizeof(double));
    MY_ASSERT(points); /*if calloc failed- print an error message and exit*/
    for(i=0; i<n*d; i++){
        PyObject *num;
        num=PyList_GetItem(points_list, i); /*num will be points_list[i]:*/
        points[i]=PyFloat_AsDouble(num); /*convert to double and save*/

    }

    /*save the centroids (from centroids_list) in 1D array of doubles*/
    clusters=calloc(k*(d+1), sizeof(double));
    MY_ASSERT(clusters); /*if calloc failed- print an error message and exit*/
    counter=0;
    /*save centroids' coordinated in clusters (with skipping one 'cell' after
    centroid-the last cell will be used for caclculations of kmeans algorithm*/
    for(i=0; i<k*(d+1); i++){  /*for each centroids*/
        for(j=0; j<d; j++){ /*for each coordinate*/
            PyObject *num; /*num will be the j's coordinate of a centroid*/
            num=PyList_GetItem(centroids_list, i-counter);
            clusters[i]=PyFloat_AsDouble(num); /*convert to double and save*/
            i++;
        }
        counter++;
    }
    /*perform kmeans algorithm on the data. now clusters contains the centroids
    that were found by kmeans algorith*/
    kmeanspp(points, clusters, n, d, k);

    free(points);/*free memory that we've allocated*/

    /*python list to save the final centroids*/
    PyObject* python_res = PyList_New(k*(d+1));

    /*convert every value of clusters to python float and save it in python_res*/
    for(i=0; i<k*(d+1); i++){
        PyObject* python_float=Py_BuildValue("d", clusters[i]);
        PyList_SetItem(python_res, i, python_float);
    }
    free(clusters); /*free memory we've allocated*/
    return python_res; /*return the result to python*/

}




/*
 * This array tells Python what methods this module has.
 * We will use it in the next structure
 */
static PyMethodDef capiMethods[]={
        {"spk",  /* a Python method name that will be used */

                /* the C-function that implements the Python
                 function and returns static PyObject* */
                (PyCFunction)spk,

                /*flags indicating parameters accepted for this function*/
                METH_VARARGS,

                /*The docstring for the function */
                PyDoc_STR("result of spk module according to the input goal")},


        {"fit",  /* a Python method name that will be used */

                 /* the C-function that implements the Python
                 function and returns static PyObject* */
                (PyCFunction)fit,
                /*flags indicating parameters accepted for this function*/
                METH_VARARGS,

                /*The docstring for the function */
                PyDoc_STR("Centroids of Kmeans++ Algorithm")},

        {NULL, NULL, 0, NULL} /* The last entry must be all NULL*/
};




/*This initiates the module using the above definitions*/
static struct PyModuleDef moduledef={
        PyModuleDef_HEAD_INIT,
        "myspkmeans", /* name of module */
        NULL,
        -1,
        capiMethods /*the PyMethodDef array from before containing
                     the methods of the extension */
};




/*
 * The PyModuleDef structure, in turn, must be passed
 * to the interpreter in the module’s initialization function.
 * The initialization function must be named PyInit_name(),
 * where name is the name of the module and should match
 * what we wrote in struct PyModuleDef.
 * This should be the only non-static item defined in the module file
 */
PyMODINIT_FUNC
PyInit_myspkmeans(void)
{
    PyObject *m;
    m=PyModule_Create(&moduledef);
    if(!m){
        return NULL;
    }
    return m;
}
