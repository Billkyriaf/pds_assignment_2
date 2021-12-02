#include <stdlib.h>
#include <stdio.h>
#include "read_points.h"

/**
 * Reads point from a file and writes them to the pointsArr. Every process reads different sections from the file that is
 * why the rank of the process and the number of points are provided.
 * @param fileName
 * @param pointsArr
 * @param rank
 * @param worldSize
 *
 */
void readFromFile(char *fileName, Point *pointsArr, int rank, int worldSize, int *pointsDimension, int *pointsPerProcess){

    FILE *file = NULL;  // File pointer
    file = fopen(fileName, "r");  // Open the file for reading

    // Check if the file is open
    if (file == NULL){
        // If not write error to the logger and exit
        // TODO add to log file
        return;
    } else {
        int totalPoints;
        char *line = NULL;
        size_t len = 0;
        ssize_t read;

        // TODO parse initial info for the file (dimension, number of points)
        // ...

        // TODO calculate the number of points per process
        // ...


        // Allocate memory for the points array
        pointsArr = (Point *) malloc (*pointsPerProcess * (sizeof (Point)));

        // Skip the first rank * numberOfPoints lines
        for (int i = 0; i < rank * (*pointsPerProcess); ++i) {
            getline(&line, &len, file);
        }

        // Start reading the next numberOfPoints lines
        for (int i = 0; i < *pointsPerProcess; ++i) {
            read = getline(&line, &len, file);

            if (read != -1){
                // TODO Parse the line data
                // ...
            } else {
                // If not write error to the logger and exit
                // TODO add to log file

                // Close the file
                fclose(file);

                // Deallocate any memory left from line
                free(line);


                return;
            }

            // Allocate memory for the coordinates array
            pointsArr[i].coordinates = (float *) malloc(*pointsDimension * sizeof (float));

            // TODO insert the coordinates
            for (int j = 0; j < *pointsDimension; ++j) {
                pointsArr[i].coordinates[j] = 0;  // FIXME update the code
            }
        }

        // Close the file
        fclose(file);

        // Deallocate any memory left from line
        free(line);
    }
}