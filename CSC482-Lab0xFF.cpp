#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>

enum Algs {
    F_BRUTEFORCE,
    F_GREEDY,
    F_ANTCOLONY
};

typedef struct Node {
    double x;
    double y;
} Node;

#define VERBOSE false

const int test = F_BRUTEFORCE;
long long unsigned int busyCount;
const long long int N_max = 10000; // adjust as needed, keep small for debugging
const int N = 10;
const int C_MAXVAL = 10;
double costMatrix[N_max][N_max] = { 0 };    /* Matrix associated with the weight of a path */
int indexes[N_max] = { 0 };                 /* array associated with tracking indexes for creating unique brute force paths */
Node coords[N_max] = { 0 };                 /* array for tracking X, Y values for Euclidean costMatrix */
double shortestPathCost = -1;               /* tracks the shortest path in a given matrix */
int shortestPath[N_max] = { 0 };

using namespace std;

void doBusyWork(void)
{
    for (int k = 0; k < N; k++)
        busyCount++;
}

void printMatrix(double costMatrix[][N_max], long long int numV)
{
    for (long long i = 0; i < numV; i++) {
        for (long long k = 0; k < numV; k++) {
            printf("%5.2f", costMatrix[i][k]);
        }
        printf("\n\n");
    }
    printf("\n\n");
}


/* Function to swap values at two pointers */
void swap(int* x, int* y)
{
    int temp;
    temp = *x;
    *x = *y;
    *y = temp;

  //  printf("x: %f, t: %f\n", *x, *y);
}

/* Function to print permutations of string
This function takes three parameters:
1. String
2. Starting index of the string
3. Ending index of the string. */
void TspBruteForce(int indexes[N_max], double costMatrix[][N_max], int l, int r)
{

    int i;
    if (l == r) {
        double cost = 0;
        int iteration = 0;
        int index = 0;

        while (iteration <= r) {
            if (VERBOSE)
                printf("%3d ", indexes[iteration]);
            cost += costMatrix[index][indexes[iteration]];
            index = indexes[iteration];
            iteration++;
        }
        if (VERBOSE)
            printf("Total cost: %f\n\n", cost);
        if (shortestPathCost < 0) {
            shortestPathCost = cost;
            memcpy(shortestPath, indexes, sizeof(indexes)*r);
        } else if (shortestPathCost > cost) {
            shortestPathCost = cost;
            memcpy(shortestPath, indexes, sizeof(indexes)*r);
        }
    }
    else
    {
        for (i = l; i < r; i++)
        {
            swap((indexes[l]), (indexes[i]));
            TspBruteForce(indexes, costMatrix, l + 1, r);
            swap((indexes[l]), (indexes[i]));
        }
    }
}



void TspGreedy(double costMatrix[][N_max], int numV)
{
    /* find shortest, check if index is already used, if not use as next index
       
    */
    double smallest = C_MAXVAL + 1;
    int index = 0;
    int tempIndex = 0;
    int visited[N_max] = { 0 };
    int path[N_max] = { 0 };
    visited[0]++;


    for (long long int j = 0; j < numV-1; j++) {
        for (long long int k = 0; k < numV; k++) {
            if (k != index) {
                if (costMatrix[index][k] < smallest && visited[k] == 0) {
                    smallest = costMatrix[index][k];
                    tempIndex = k;
                }
            }
        }
        index = tempIndex;
        visited[index]++;
        path[j+1] = index;
        smallest = C_MAXVAL + 1;
        printf("Selected: %d\n", index);
        printf("Visited: ");
        for (int l = 0; l < numV; l++) {
            if (visited[l] > 0)
                printf("%d, ", l);
        }
    }
    printf("Generation done. Final path: ");
    for (int k = 0; k < numV; k++)
        printf("%d -> ", path[k]);
    printf("0");
    printf("\n\n");

}

void TspAntColony(void)
{

}

void GenerateRandomCostMatrix(double costMatrix[][N_max], int numV, int maxValue)
{
    int num;

    for (long long int i = 0; i < numV; i++) {
        for (long long int j = 0; j < numV; j++) {
            if (i == j)
                costMatrix[i][j] = 0;
            else {
                num = rand() % maxValue;
                costMatrix[i][j] = num;
                costMatrix[j][i] = num;
            }
        }
    }
}


void GenerateRandomEuclideanCostMatrix(double costMatrix[][N_max], int numV, int maxValue)
{
   /* Given a number of vertices, and a max x& y coordinate, generate randome x, y 
   coordinates for each vertext, then populate the cost matrix by calculating the 
   distance between each pair of points. */

    for (long long int k = 0; k < numV; k++) {
        coords[k].x = rand() % maxValue;
        coords[k].y = rand() % maxValue;
    }

    for (long long int l = 0; l < numV; l++) {
        for (long long int m = 0; m < numV; m++) {
            if (l != m)
                costMatrix[l][m] = sqrt(pow(coords[m].x - coords[l].x, 2) + pow(coords[m].y - coords[l].y, 2));
        }
    }
}

void GenerateRandomCircularGraphCostMatrix()
{

}




int main(int argc, char** argv) {

    double trialsTime_max = .250; // in seconds
    long long int trialsCount_max = 1000000,
        N_min = 1,
        trial;
    clock_t splitTimeStamp, trialSetStart;
    const long long  Array_max = 100000000;
    double splitTime, trialSetCount, trialSetTime, dummySetCount, dummySetTime, averageTrialTime, averageDummyTrialTime, estimatedTimePerTrial;

    double times[100] = { 0 };
    int index = 1;
    srand(time(NULL));

    // If you are printing a results table directly, print your header line here.
    printf("+-------------------------------------------------------------------------------------------+\n");
    printf("| %20s | %20s | %20s | %20s |\n", "N", "Time", "Exprmntl Dbl Ratio", "Expected Dbl Ratio");
    printf("+-------------------------------------------------------------------------------------------+\n");
    // power of 2 | N | Measured Time | Measured Doubling Ratio | Expected Doubling Ratio |Busy Count | Measured Time / Busy Count
    // For each size N of input we want to test -- typically start N at 1 and double each repetition
    //for (long long int n=1; n<N_max; n=2*n ) {
    for (long long int n = 4; n < N_max; n++) {
        /********* any preparations, test input generation, to be used for entire trial set ********/

            /* test matrix */

 //   GenerateRandomCostMatrix(costMatrix, N, C_MAXVAL);
        GenerateRandomEuclideanCostMatrix(costMatrix, n, C_MAXVAL);
        if (VERBOSE)
            printMatrix(costMatrix, n);


        splitTime = 0.0;
        // get timestamp before set of trials are run:
        trialSetStart = clock();
        // For each trial trialNumber=1 to Number of Trials per input size:
        for (trial = 0; trial < trialsCount_max && splitTime < trialsTime_max; trial++) {
            /******** any preparations, test input generation, for a single trial run *********/
            busyCount = 0;
            /**** >>>> Call the algorithm function we are testing, using the generated test input <<<<< ******/




            switch (test) {
            case F_BRUTEFORCE:
                for (int k = 0; k < n - 1; k++) {
                    indexes[k] = k + 1;
                }

                TspBruteForce(indexes, costMatrix, 0, n - 1);

                if (VERBOSE) {
                    printf("Shortest Path: ");
                    for (int k = 0; k < n; k++) {
                        if (k == n - 1)
                            printf("%d\n\n", shortestPath[k]);
                        else
                            printf("%d -> ", shortestPath[k]);
                    }
                    printf("Shortest Path Cost: %f\n", shortestPathCost);
                }
                shortestPathCost = -1;
                break;
            case F_GREEDY:
                TspGreedy(costMatrix, n);
                break;
            case F_ANTCOLONY:
                TspAntColony();
                break;
            default:
                break;
            }


            /******* do any "clean up" after running your algorithm ********/

            // get split time -- "split time" just means the total time so far for this tiral set
            splitTimeStamp = clock(); // 
            // split time is the difference between current split timestamp and the starting time stamp for trial set
            splitTime = (splitTimeStamp - trialSetStart) / (double)CLOCKS_PER_SEC; // CLOCK_PER_SEC define time.h 
        }
        trialSetCount = trial; // value of trial when loop ends is how many we did
        trialSetTime = splitTime; // total time for trial set is the last split time
        averageTrialTime = trialSetTime / trialSetCount; // this is the average tiem per trial, including any prep/overhead



        /********* NOW DO A "DUMMY" TRIAL SET TO ESTIMATE THE OVERHEAD TIME ********/
        /* We can just repeat this trialSetCount times, which we know should be faster than above */

        splitTime = 0.0;
        // get timestamp before set of dummy trials are run:
        trialSetStart = clock();
        for (trial = 0; trial < trialSetCount && splitTime < trialsTime_max; trial++) {

            /******** any preparations, test input generation, for a single trial run *********/

            /**** DO NOT Call the algorithm function!!! ******/

            /******* do any "clean up" after running your algorithm ********/

            // get split time -- "split time" just means the total time so far for this tiral set
            splitTimeStamp = clock(); // 
            // split time is the difference between current split timestamp and the starting time stamp for trial set
            splitTime = (splitTimeStamp - trialSetStart) / (double)CLOCKS_PER_SEC; // CLOCK_PER_SEC define time.h 
        }
        dummySetCount = trial; // value of trial when loop ends is how many we did, should be same as trialsSetCount
        dummySetTime = splitTime; // total time for dummy trial set is the last split time
        averageDummyTrialTime = dummySetTime / dummySetCount; // this is the average tiem per dummy trial, including any prep/overhead


        estimatedTimePerTrial = averageTrialTime - averageDummyTrialTime; // should be a decent measurement of time taken to run your algorithm

        times[index] = estimatedTimePerTrial;

        if (n == 1)
            printf("| %20llu | %20.2f | %20.2f | %20.2f |\n", n, times[index], times[index] / times[index-1], n);
        else {
            switch (test) {
            case F_BRUTEFORCE:
                printf("| %20llu | %20.2f | %20.2f | %20.2f |\n", n, times[index], times[index] / times[index - 1], pow(n,2)/pow(n-1,2));
                break;
            case F_GREEDY:
                printf("| %20llu | %20.2f | %20.2f | %20.2f |\n", n, times[index], times[index] / times[index - 1], 2);
                break;
            case F_ANTCOLONY:
                printf("| %20llu | %20.2f | %20.2f | %20.2f |\n", n, times[index], times[index] / times[index - 1], 2);
                break;
            default:
                break;
            }
        }
        index++;

        busyCount = 0;

        /************ save and/or print your results here to generate your table of results **************/
        // You should be able to print one row of your results table here.
        // Calculate doubling ratios and any other results for desired columns in table.
        // Consider additional columns you include in a more detailed table that may help with debugging or just making sense of results.
        // May be a good idea to print a row of results to output terminal window or a file, so if your program crashes you still have results to look at
        // use the fflush(stdout) or fflush(filePointer) or equivalent in C++, to force output so nothing is lost if a crash occurs
        // can also save data to array(s) for later printing/processing

    }
}

//////////////////////////////////
/////////////////////////////////