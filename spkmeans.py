import myspkmeans as mk
import numpy as np
import pandas as pd
import sys

#define a dictionary to determine goal's value
goals={"spk":0, "wam":1, "ddg": 2, "lnorm":3, "jacobi": 4};

"""
input:
int k - number of clusters
points - 2-D numpy array of points

choose the inices of first k centroids according to Kmeans++ algorithm
output: tuple - list of indices of the centroids that were chosen, 
        and list of the centroids.
"""
def kmeansStep1(k, points):
    z=0
    n=points.shape[0] #number of points
    d=points.shape[1] #dimension of points
    np.random.seed(0)
    indices=[] #this list will contain the indices of the first k centroids
    random_index=np.random.choice(n) #choose a random index from 0 to n-1
    indices.append(random_index)
    u1=points[random_index] #u1 is the point that it's index is random_index we
                            # chose.
    centroids=np.zeros(shape=(k,d)) #create a 2D zero numpy array (size: kxd)
    centroids[0]=u1 #first row in centroids is u1
    z=1 #increment z
    D = np.zeros(shape=(1, n)) #create a 1D zero numpy array (size: 1xn)
    while(z<k): #until we will choose more k-1 centroids
        for i in range(0,n): #for each point
            di=np.sum((centroids[z-1]-points[i])**2)
            if (di<D[0][i] or z<=1):
                D[0][i]=di #di=min(points[i]-centroids[j]) for 0<=j<=z
        sum=np.sum(D, axis=1) #sum of elements in D
        propbability=np.divide(D,sum) #divide each D[i] by sum

        # choose indexramdomly according to the probiblity we calaculated
        random_index = np.random.choice(n, p=propbability[0])

        #add a new point to centoroids, according to the index that was chosen
        centroids[z]=points[random_index]

        #add the index that was chosen to indices
        indices.append(random_index)
        z+=1 #increment z by 1

    #convert centoroids to 1D python list
    oneDimCent=np.reshape(centroids, (np.product(centroids.shape),))
    res=oneDimCent.tolist()

    return (indices,res)


"""
input:
vector - row of 2-D numpy array
int n - size of vector
This method prints the row's elements, seperated by ","
"""
def printRow(vector, n):
    for i in range(0,n-1): #print first n-1 elements
        if (vector[i] < 0 and vector[i] > -0.00005): #if element has a small
                                                     #negative value
            vector[i] = 0.0000 #print it as zero
        print("%.4f,"%vector[i],end="");

    # print last element - same as before, without "," in the end
    if (vector[n-1] < 0 and vector[n-1] > -0.00005):
        vector[n-1] = 0.0000
    print("%.4f"%vector[n-1],end="")



"""
input:
res - 2-D numpy array (nxn)
int n - size of res

This method prints the rows of res, row by row.
"""
def printResult(res, n):
    for i in range(0,n-1):
        printRow(res[i],n);
        print("\n", end="");
    printRow(res[n-1],n); #print last row




#the program flow

k=int(sys.argv[1]) #covert k input that that passed by the user to int
goal=goals[str(sys.argv[2])] #choose goal value (int) according to the key (string)

df1=pd.read_csv(sys.argv[3], delimiter=",", header=None) #read the file

points=df1.to_numpy() #convert the data-frame to numpy array
n=points.shape[0] #n=number of points
d=points.shape[1] #d=dimension of points

#convert the points to 1-D python list
oneDimPoints=np.reshape(points, (np.product(points.shape),));
listOfpoints=oneDimPoints.tolist();

#calc the desrired result acoording to the goal through C.
newData=mk.spk(listOfpoints, n, d,k,goal);

if goal==0: #goal=spk
    k = int((len(newData)) / n) #new data is T (size: n*k)
    T = np.reshape(newData, (-1, k)) #convert T to n*k matrix

    # find first k centroids and indices
    indices, listOfcentroids=kmeansStep1(k,T)
    #print(*indices, sep=',')

    #perform kmean algorithm with newData, with the first k centroids we chose
    centroids = mk.fit(newData, listOfcentroids, n, k, k)

    #save the centroids in 2D - numpy array.
    res = np.reshape(centroids, (-1, k + 1))[0:, 0:k]
    print(*indices, sep=',') #print the indices of the first centroids
    printResult(res, k) #print the final centroids

elif goal==4: #goal=jacobi
    # save the eigenvalues and eigenvectors in 2D - numpy array.
    res = np.reshape(newData, (-1, n))
    printRow(res[n],n) #print eigenvalues
    print("\n", end="")
    printResult(res,n) #print eigenvectors

elif goal==2: #goal=ddg
    a = np.zeros((n, n), float);
    np.fill_diagonal(a, newData);
    printResult(a,n) #print D

else: #goal=wam or goal=lnorm
    res = np.reshape(newData, (-1, n))
    printResult(res, n) #print Weighted Adjacency Matrix or
                        # Normalized Graph Laplacian















