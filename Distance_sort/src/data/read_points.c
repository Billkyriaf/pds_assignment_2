#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
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

        // Parse initial info for the file (dimension, number of points)
        getline(&line, &len, file);
        totalPoints = (int) strtol(line, NULL, 10);

        getline(&line, &len, file);
        *pointsDimension = (int) strtol(line, NULL, 10);


        // Calculate the number of points per process
        *pointsPerProcess = totalPoints / worldSize;


        // Allocate memory for the points array
        pointsArr = (Point *) malloc (*pointsPerProcess * (sizeof (Point)));

        // Skip the first rank * numberOfPoints lines
        for (int i = 0; i < rank * (*pointsPerProcess); ++i) {
            getline(&line, &len, file);
        }

        // Start reading the next numberOfPoints lines
        for (int i = 0; i < *pointsPerProcess; ++i) {

            // Read the next line from the file
            read = getline(&line, &len, file);

            // Allocate memory for the coordinates array
            pointsArr[i].coordinates = (float *) malloc(*pointsDimension * sizeof (float));

            if (read != -1){
                // Points to the first char of the line initially
                char *p = line;

                // For every number in the string
                for(int j = 0; j < *pointsDimension; j++) {

                    // If the first character of the string is digit or sign
                    if ( isdigit(*p) || ( (*p=='-'||*p=='+') && isdigit(*(p+1)) )) {

                        // Insert the coordinates to the array
                        // strtof() reads the hole number and points p to the next non digit char in the string
                        pointsArr[i].coordinates[j] = strtof(p, &p);

                    } else {
                        // If *p is a non digit char...

                        // Repeat the j because no number found yet
                        j--;

                        // Go to the next char
                        p++;
                    }
                }

            } else {
                break;
            }
        }

        // Close the file
        fclose(file);

        // Deallocate any memory left from line
        free(line);
    }
}