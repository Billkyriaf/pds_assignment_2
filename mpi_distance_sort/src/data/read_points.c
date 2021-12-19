#include <stdlib.h>
#include <stdio.h>
#include "read_points.h"

/**
 * Reads point from a file and writes them to the pointsArr. Every process reads different sections from the file that is
 * why the rank of the process and the number of points are provided.
 *
 * @param fileName         The path to the data file
 * @param pointsArr        The array that holds all the points and their info for each process
 * @param rank             The rank of the process that calls the function
 * @param worldSize        The total number of processes
 * @param pointsDimension  The dimension of each point (number of coordinates)
 * @param pointsPerProcess The number of points each process has to manage
 *
 */
Point *readFromFile(char *fileName, int rank, int worldSize, int *pointsDimension, int *pointsPerProcess){
    FILE *file = NULL;  // File pointer
    file = fopen(fileName, "rb");  // Open the file for reading

    // Check if the file is open
    if (file == NULL){
        // If not write error to the logger and exit
        // TODO add to log file
        exit(-1);

    } else {
        int64_t totalPoints;  // The total number of points
        int64_t d;  // The dimension of the space

        // Read the dimension and the number of points
        fread(&d, sizeof(int64_t), 1, file);
        *pointsDimension = (int)d;

        // Read the total points from the file
        fread(&totalPoints, sizeof(int64_t), 1, file);

        // Calculate the number of points per process
        *pointsPerProcess = (int)totalPoints / worldSize;

        // printf("Dimension: %ld, Points per process: %ld\n", *pointsDimension, *pointsPerProcess);

        // Allocate memory for the points array
        Point *pointsArr = (Point *) malloc (*pointsPerProcess * (sizeof (Point)));

        // Skip the first rank * numberOfPoints lines. This is how every process ends up with different points.
        // All the process read from the same file, but they all read different sections of the file.
        fseek(file, rank * (*pointsDimension) * (*pointsPerProcess) * (long int)sizeof (double), SEEK_CUR);

        // Start reading the next numberOfPoints lines
        for (int i = 0; i < *pointsPerProcess; ++i) {

            // Allocate memory for the coordinates array
            pointsArr[i].coordinates = (double *) malloc (*pointsDimension * sizeof (double));

            // Read the coordinates fread function can read multiple numbers at a time and stores them to an array
            fread(pointsArr[i].coordinates, sizeof(double), *pointsDimension, file);
        }

        // Close the file
        fclose(file);

        // Return the reference to the points array
        return pointsArr;
    }
}