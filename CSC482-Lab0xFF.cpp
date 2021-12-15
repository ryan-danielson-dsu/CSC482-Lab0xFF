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
    F_ANTCOLONY,
    F_ANTCOLONY_VS_EXACT,
    F_GREEDY_VS_EXACT,
    F_ALL
};
    
typedef struct Node {
    double x;
    double y;
} Node;

#define VERBOSE true
#define SHOWTABLE false
#define EUCLIDEAN false

const int TEST = F_GREEDY_VS_EXACT;
const int TIME_STEPS = 1000;
const int NUM_ANTS = 100;
const float DECAY = 0.01;
const float PHEROMONE_FACTOR = 1;
const int RADIUS = 100;
long long unsigned int busyCount;
const long long int N_max = 1000; // adjust as needed, keep small for debugging
const int N = 250;
const int C_MAXVAL = 10000;
double costMatrix[N_max][N_max] = { 0 };    /* Matrix associated with the weight of a path */
int indexes[N_max] = { 0 };                 /* array associated with tracking indexes for creating unique brute force paths */
Node coords[N_max] = { 0 };                 /* array for tracking X, Y values for Euclidean costMatrix */
double shortestPathCost = -1;               /* tracks the shortest path in a given matrix */
int shortestPath[N_max] = { 0 };


using namespace std;

void doBusyWork(void)
{
    for (int k = 0; k < N; k++)
        for (int j = 0; j < N; j++)
            busyCount++;
}

void printMatrix(double costMatrix[][N_max], long long int numV)
{
    for (long long i = 0; i < numV; i++) {
        for (long long k = 0; k < numV; k++)
            printf("%8.2f", costMatrix[i][k]);
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
}

/* 
 * BruteForce based on the string permutation algorithm from Geeks for Geeks:
 * https://www.geeksforgeeks.org/write-a-c-program-to-print-all-permutations-of-a-given-string/
 * Recursively operates on a given set of indexes (representing nodes), summing the cost
 * by the costMatrix. Returns shortest path after all permutations have been calculated.
 * Uses l to indicate left index, where r is the base case for when all permutations have been
 * discovered. Swaps all indexes in indexes, using left index as starting node to swap.
 */
void TspBruteForce(int indexes[N_max], double costMatrix[][N_max], int l, int r)
{
    int i;
    if (l == r) {
        double cost = 0;
        int iteration = 0;
        int index = 0;

        while (iteration <= r) {
            cost += costMatrix[index][indexes[iteration]];
            index = indexes[iteration];
            iteration++;
        }

        if (shortestPathCost < 0) {
            shortestPathCost = cost;
            memcpy(shortestPath, indexes, sizeof(indexes)*r);
        } else if (shortestPathCost > cost) {
            shortestPathCost = cost;
            memcpy(shortestPath, indexes, sizeof(indexes)*r);
        }
    } else {
        for (i = l; i < r; i++) {
            swap((indexes[l]), (indexes[i]));
            TspBruteForce(indexes, costMatrix, l + 1, r);
            swap((indexes[l]), (indexes[i]));
        }
    }
}

/* 
 * The purpose of this function is to find the worst-case path in a given matrix. Same as BruteForce
 * except tests for largest value instead of smallest.
 */
void TspBruteForceWorstCase(int indexes[N_max], double costMatrix[][N_max], int l, int r)
{
    int i;
    if (l == r) {
        double cost = 0;
        int iteration = 0;
        int index = 0;

        while (iteration <= r) {
            cost += costMatrix[index][indexes[iteration]];
            index = indexes[iteration];
            iteration++;
        }
        if (shortestPathCost < 0) {
            shortestPathCost = cost;
            memcpy(shortestPath, indexes, sizeof(indexes) * r);
        }
        else if (shortestPathCost < cost) {
            shortestPathCost = cost;
            memcpy(shortestPath, indexes, sizeof(indexes) * r);
        }
    } else {
        for (i = l; i < r; i++) {
            swap((indexes[l]), (indexes[i]));
            TspBruteForceWorstCase(indexes, costMatrix, l + 1, r);
            swap((indexes[l]), (indexes[i]));
        }
    }
}

/* 
 * Greedy algorithm takes the shortest available path from any given node, regardless of where
 * this may lead us.
 */

void TspGreedy(double costMatrix[][N_max], int numV)
{
    /*
     *   find shortest, check if index is already used, if not use as next index
     */

    double smallest = C_MAXVAL + 1;
    int index = 0;
    int tempIndex = 0;
    int visited[numV + 1] = { 0 };
    int path[numV + 1] = { 0 };
    
    /* visit home node default */
    visited[0]++;

    for (long long int j = 0; j < numV - 1; j++) {
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
        path[j + 1] = index;
        smallest = C_MAXVAL + 1;
    }
    double cost = 0;
    int pathIndex = 0;
    int acc = 1;
    while (acc < numV+1) {
        cost += costMatrix[pathIndex][path[acc]];
        pathIndex = path[acc];
        acc++;   
    }
    cost += costMatrix[pathIndex][0];


    if (VERBOSE) {
        printf("[Greedy]     Best path cost: %f\n", cost);
        printf("[Greedy]     Best path: ");
        for (int l = 0; l < numV; l++)
            if (visited[l] > 0)
                printf("%d -> ", path[l]);
        printf("0\n");
    }
}


void TspAntColony(double costMatrix[][N_max], int timeSteps, int numV, float pheromoneFactor, int m, float decayFactor)
{
    int visited[numV] = { 0 };
    int currentPath[numV] = { 0 };
    int bestPath[numV] = { 0 };
    double pheromoneMatrix[numV][numV] = { 0 };
    double newPheromoneMatrix[numV][numV] = { 0 };

    const int homeNode = 0;
    int currentNode = 0;
    int nextNode = 0;
    int maxUnchangedTimeSteps = 50;
    int unchangedTimeSteps = 0;
    double totalAttraction = 0;
    double cummulativeProbability = 0;
    double minPathCostSoFar = 0;
    double pathCost = 0;
    double edgeSelectionProbability = 0;


    for (long long int step = 0; step < timeSteps; step++) {
        if (unchangedTimeSteps > maxUnchangedTimeSteps) {
            break;

        }
        memset(newPheromoneMatrix, 0, sizeof(newPheromoneMatrix));

        for (long long int j = 0; j < m; j++) {
            /* upkeep for those ants */
            pathCost = 0;
            currentPath[0] = homeNode;
            memset(visited, 0, sizeof(visited));
            visited[homeNode] = 1;
            
            doBusyWork();

            
            if (step > 0) {
                for (int k = 1; k < numV; k++) {
                    currentNode = currentPath[k - 1];
                    totalAttraction = 0;

                    for (nextNode = 0; nextNode < numV - 1; nextNode++) {
                        if (!visited[nextNode] && nextNode != currentNode) {
                            if (costMatrix[currentNode][nextNode] > 0.0) {
                                totalAttraction += (1 + pheromoneMatrix[currentNode][nextNode]) / costMatrix[currentNode][nextNode];
                            }
                        }
                    }

                    /* probability of choosing an edge = attraction / sum of attraction of all remaining available edges */
                    double q = (double)rand() / (double)RAND_MAX;
                    cummulativeProbability = 0;

                    for (nextNode = 0; nextNode < numV - 1; nextNode++) {
                        if (!visited[nextNode] && nextNode != currentNode) {
                            if (costMatrix[currentNode][nextNode] > 0.0 && totalAttraction > 0) {
                                edgeSelectionProbability = ((1 + pheromoneMatrix[currentNode][nextNode]) / costMatrix[currentNode][nextNode]) / totalAttraction;
                                cummulativeProbability += edgeSelectionProbability;
                            }
                            else {
                                edgeSelectionProbability = (1 / costMatrix[currentNode][nextNode]) / totalAttraction;
                                cummulativeProbability += edgeSelectionProbability;
                            }
                        }
                        if (q < cummulativeProbability)
                            break;
                    }
                    currentPath[k] = nextNode;
                    visited[nextNode] = 1;
                }
            }

            /* update path cost for new path */
            for (long long int u = 0; u < numV-1; u++)
                pathCost += costMatrix[currentPath[u]][currentPath[u + 1]];
            pathCost += costMatrix[currentPath[numV-1]][homeNode];

            /* Test if path cost is less than predecessors */
            if (pathCost < minPathCostSoFar || minPathCostSoFar == 0) {
                minPathCostSoFar = pathCost;
                for (int k = 0; k < numV; k++)
                    bestPath[k] = currentPath[k];
            } else {
                unchangedTimeSteps++;
            } 

            /* add new pheromone for path */
            for (long long int n = 0; n < numV-1; n++) {
                currentNode = currentPath[n];
                nextNode = currentPath[(n+1) % numV];
                newPheromoneMatrix[currentNode][nextNode] = (newPheromoneMatrix[currentNode][nextNode] + pheromoneFactor) / pathCost;
            }
         
        }

        /* decay pheromone matrix */
        for (long long int k = 0; k < numV-1; k++) {
            for (long long int h = 0; h < numV-1; h++) {
                pheromoneMatrix[k][h] = pheromoneMatrix[k][h] * decayFactor;
                pheromoneMatrix[k][h] = pheromoneMatrix[k][h] + newPheromoneMatrix[k][h];
            }
        }
    }
    if (VERBOSE) {
        printf("[Ant Colony] Best path cost: %f\n", minPathCostSoFar);
        printf("[Ant Colony] Best path: ");
        for (int p = 0; p < numV; p++)
            printf("%d -> ", bestPath[p]);
        puts("0");
    }
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

void GenerateRandomCircularGraphCostMatrix(double costMatrix[][N_max], int numV, int maxValue, int radius)
{
    double stepAngle = 2 * M_PI / numV;
    int nextNode;
    Node temp;
    int used[numV] = { 0 };
    int bestPath[numV] = { 0 };


    /* Create the coords by step */
    for (long long int s = 0; s < numV; s++) {
        coords[s].x = radius * sin(s * stepAngle);
        coords[s].y = radius * cos(s * stepAngle);
    }

    

    /* best path will be 0,1,2,3...N-1,N before shuffle */
    for (long long int k = 0; k < numV; k++)
        bestPath[k] = k;

    /*shuffle the coordinate indexes around */
    
    for (long long int s = 1; s < numV/2; s+=2) {
        nextNode = 1 + (rand() % numV - 1);
        while (used[nextNode] >= 1 || nextNode == 0 || s == nextNode)
            nextNode = 1 + (rand() % numV - 1);

        temp = coords[nextNode];
        coords[nextNode] = coords[s];
        coords[s] = temp; 

        swap(bestPath[nextNode], bestPath[s]);

        used[nextNode]++;
        used[s]++;
    }  

    double shortestDistance = -1;
    /* calculate distances and populate costMatrix */
    for (long long int l = 0; l < numV; l++) {
        for (long long int m = 0; m < numV; m++) {
            if (l != m) {
                costMatrix[l][m] = sqrt(pow(coords[m].x - coords[l].x, 2) + pow(coords[m].y - coords[l].y, 2));
                if (costMatrix[l][m] < shortestDistance || shortestDistance == -1)
                    shortestDistance = costMatrix[l][m];
            }
        }
    }

    /* Output for testing(comparison purposes */
    if (VERBOSE) {
        printf("\nBest path: ");
        for (long long int k = 0; k < numV; k++) {
            printf("%d -> ", bestPath[k]);
        }
        printf("0\n");

        printf("Best path coords: ");
        for (long long int k = 0; k < numV - 1; k++)
            printf("(%.2f,%.2f) -> ", coords[k].x, coords[k].y);
        printf("(%.2f,%.2f)\n", coords[numV].x, coords[numV].y);

        printf("Expected Cost: (%f * %d) = %f\n\n", shortestDistance, numV, (double)(shortestDistance * numV));
    }
}

unsigned long long int factorial(unsigned int n)
{
    unsigned long long int res = 1;
    for (unsigned long long int i = 2; i <= n; i++)
        res *= i;
    return res;
}



int main(int argc, char** argv) {

    double trialsTime_max = .250; // in seconds
    long long int trialsCount_max = 1,// 1000000,
        N_min = 1,
        trial;
    clock_t splitTimeStamp, trialSetStart;
    const long long  Array_max = 100000000;
    double splitTime, trialSetCount, trialSetTime, dummySetCount, dummySetTime, averageTrialTime, averageDummyTrialTime, estimatedTimePerTrial;

    double times[N_max] = { 0 };
    int index = 1;
    srand(time(NULL));

    // If you are printing a results table directly, print your header line here.
    printf("+-------------------------------------------------------------------------------------------+\n");
    printf("| %20s | %20s | %20s | %20s |\n", "N", "Time", "Exprmntl Dbl Ratio", "Expected Dbl Ratio");
    printf("+-------------------------------------------------------------------------------------------+\n");
    // power of 2 | N | Measured Time | Measured Doubling Ratio | Expected Doubling Ratio |Busy Count | Measured Time / Busy Count
    // For each size N of input we want to test -- typically start N at 1 and double each repetition
    //for (long long int n=1; n<N_max; n=2*n ) {
    for (long long int n = 3; n < N_max; n++) {
        /********* any preparations, test input generation, to be used for entire trial set ********/

        /* test matrix */


        splitTime = 0.0;
        // get timestamp before set of trials are run:
        trialSetStart = clock();
        // For each trial trialNumber=1 to Number of Trials per input size:
        for (trial = 0; trial < trialsCount_max && splitTime < trialsTime_max; trial++) {
            /******** any preparations, test input generation, for a single trial run *********/
            busyCount = 0;
            /**** >>>> Call the algorithm function we are testing, using the generated test input <<<<< ******/




            switch (TEST) {
            case F_GREEDY_VS_EXACT:
                for (int sqrTrials = 0; sqrTrials < 10; sqrTrials++) {
                    if (EUCLIDEAN)
                        GenerateRandomEuclideanCostMatrix(costMatrix, n, RADIUS);
                    else
                        GenerateRandomCircularGraphCostMatrix(costMatrix, n, C_MAXVAL, RADIUS);

                    if (!SHOWTABLE && VERBOSE) {
                        puts("+---------------------------------------------------------------------------------------------+");
                        printf("N: %d\n", n);
                        puts("+---------------------------------------------------------------------------------------------+");
                    }
                    for (int k = 0; k < n - 1; k++) {
                        indexes[k] = k + 1;
                    }

                    TspBruteForce(indexes, costMatrix, 0, n - 1);
                    printf("[BruteForce] Best path cost: %f \n", shortestPathCost);
                    printf("[BruteForce] Best path: 0 -> ");
                    for (int k = 0; k < n - 1; k++)
                        printf("%d -> ", shortestPath[k]);
                    puts("0");

                    memset(shortestPath, 0, sizeof(shortestPath));
                    shortestPathCost = -1;
                    TspGreedy(costMatrix, n);
                }
                break;
            case F_ANTCOLONY_VS_EXACT:
                for (int sqrTrials = 0; sqrTrials < 10; sqrTrials++) {
                    if (EUCLIDEAN)
                        GenerateRandomEuclideanCostMatrix(costMatrix, n, RADIUS);
                    else
                        GenerateRandomCircularGraphCostMatrix(costMatrix, n, C_MAXVAL, RADIUS);

                    if (!SHOWTABLE && VERBOSE) {
                        puts("+---------------------------------------------------------------------------------------------+");
                        printf("N: %d\n", n);
                        puts("+---------------------------------------------------------------------------------------------+");
                    }
                    for (int k = 0; k < n - 1; k++) {
                        indexes[k] = k + 1;
                    }

                    TspBruteForce(indexes, costMatrix, 0, n - 1);
                    printf("[BruteForce] Best path cost: %f \n", shortestPathCost);
                    printf("[BruteForce] Best path: 0 -> ");
                    for (int k = 0; k < n - 1; k++)
                        printf("%d -> ", shortestPath[k]);
                    puts("0");

                    memset(shortestPath, 0, sizeof(shortestPath));
                    shortestPathCost = -1;
                    TspAntColony(costMatrix, TIME_STEPS, n, PHEROMONE_FACTOR, NUM_ANTS, DECAY);

                }
                break;
            case F_BRUTEFORCE:
                if (EUCLIDEAN)
                    GenerateRandomEuclideanCostMatrix(costMatrix, n, RADIUS);
                else
                    GenerateRandomCircularGraphCostMatrix(costMatrix, n, C_MAXVAL, RADIUS);

                for (int k = 0; k < n - 1; k++) {
                    indexes[k] = k + 1;
                }

                TspBruteForce(indexes, costMatrix, 0, n - 1);

                if (VERBOSE) {
                    printf("BrutForce shortest Path: 0 -> ");
                    for (int k = 0; k < n; k++) {
                        if (k == n - 1)
                            printf("%d\n", shortestPath[k]);
                        else
                            printf("%d -> ", shortestPath[k]);
                    }
                    printf("BruteForce shortest Path Cost: %f\n\n", shortestPathCost);
                    printMatrix(costMatrix, n);
                    if (!SHOWTABLE)
                        puts("+---------------------------------------------------------------------------------------------+");
                }

                shortestPathCost = -1;
                memset(shortestPath, 0, sizeof(shortestPath));
                break;
            case F_GREEDY:
                if (EUCLIDEAN)
                    GenerateRandomEuclideanCostMatrix(costMatrix, n, RADIUS);
                else
                    GenerateRandomCircularGraphCostMatrix(costMatrix, n, C_MAXVAL, RADIUS);

                TspGreedy(costMatrix, n);
                if (!SHOWTABLE && VERBOSE)
                    puts("+---------------------------------------------------------------------------------------------+");
                break;
            case F_ANTCOLONY:
                if (EUCLIDEAN)
                    GenerateRandomEuclideanCostMatrix(costMatrix, n, RADIUS);
                else
                    GenerateRandomCircularGraphCostMatrix(costMatrix, n, C_MAXVAL, RADIUS);

                TspAntColony(costMatrix, TIME_STEPS, n, PHEROMONE_FACTOR, NUM_ANTS, DECAY);
                break;
            case F_ALL:

                /* This case only exists for side by side algorithm comparisons */
                if (EUCLIDEAN)
                    GenerateRandomEuclideanCostMatrix(costMatrix, n, RADIUS);
                else
                    GenerateRandomCircularGraphCostMatrix(costMatrix, n, C_MAXVAL, RADIUS);

                printf("+-------------------------\n");
                TspAntColony(costMatrix, TIME_STEPS, n, PHEROMONE_FACTOR, NUM_ANTS, DECAY);
               
                
                for (int k = 0; k < n - 1; k++) {
                    indexes[k] = k + 1;
                }

                TspBruteForce(indexes, costMatrix, 0, n - 1);
                printf("[BruteForce] Best path cost: %f \n", shortestPathCost);
                printf("[BruteForce] Best path: 0 -> ");
                for (int k = 0; k < n - 1; k++)
                    printf("%d -> ", shortestPath[k]);
                puts("0");

                memset(shortestPath, 0, sizeof(shortestPath));
                shortestPathCost = -1;
               
                for (int k = 0; k < n - 1; k++) {
                    indexes[k] = k + 1;
                }
                TspBruteForceWorstCase(indexes, costMatrix, 0, n - 1);
                printf("[WorstCase]  Worst path cost: %f \n", shortestPathCost);
                printf("[WorseCase]  Worst path: 0 -> ");
                for (int k = 0; k < n - 1; k++)
                    printf("%d -> ", shortestPath[k]);
                puts("0");

                memset(shortestPath, 0, sizeof(shortestPath));
                shortestPathCost = -1;

                TspGreedy(costMatrix, n);

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
        
        if (SHOWTABLE) {
            if (n > 3) {
                switch (TEST) {
                case F_BRUTEFORCE:
                    printf("| %20llu | %20.6f | %20.6f | %20.6f |\n", n, times[index], times[index] / times[index - 1], (double)(factorial(n)/factorial(n-1)));
                    break;
                case F_GREEDY:
                    printf("| %20llu | %20.6f | %20.6f | %20.6f |\n", n, times[index], times[index] / times[index - 1], pow(n, 2)/pow(n-1,2));
                    break;
                case F_ANTCOLONY:
                    printf("| %20llu | %20.6f | %20.6f | %20.6f |\n", n, times[index], times[index] / times[index - 1], pow(n-1, 2) / pow(n - 2, 2));
                    break;
                default:
                    break;
                }
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