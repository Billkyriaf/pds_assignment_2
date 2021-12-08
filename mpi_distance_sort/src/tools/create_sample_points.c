#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define RANGE 200.0

/**
 * Reads a data file and prints the points for debugging
 * @param filename The file to read
 */
void testRead(char *filename){
    // Test read the file and print the numbers
    int64_t nPoints;
    int64_t dim;

    // Open the file for read
    FILE *fh = fopen (filename, "rb");

    // If the file is there
    if (fh != NULL) {

        // Read the dimension and the number of points
        fread(&dim, sizeof(int64_t), 1, fh);
        fread(&nPoints, sizeof(int64_t), 1, fh);

        // Allocate memory for one point
        double *point = (double *)malloc(sizeof(double) * dim);

        for(int i = 0; i < nPoints; i++){
            // Read the next dim x 64 bits from the file (this reads all the coordinates for one point)
            fread (point, dim * sizeof (double), dim, fh);

            // Print the point's coordinates
            for(int j = 0; j < dim; j++){
                printf("%.10f ", point[j]);
            }

            printf("\n");
        }

        // Deallocate memory and close the file
        free(point);
        fclose (fh);
    }
}

/**
 * Creates a binary file with random points.
 *
 * The command line arguments must be:
 *     1. Space dimension
 *     2. Number of points
 *     3. Path to the file to create + file name
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv){

    // Check for the arguments
    if(argc != 4){
        printf("Wrong number of arguments. Required:\n");
        printf("    1. Space dimension\n");
        printf("    2. Number of points\n");
        printf("    3. Path to the file to create + file name\n");
        return -1;
    }

    // Init the random number generator
    srand(time(NULL));


    // Both dimension and numberOfPoints are 64bit integers
    int64_t spaceDimension= atoi(argv[1]);
    int64_t numberOfPoints = atoi(argv[2]);

    // Create the new binary file
    FILE *fh = fopen (argv[3], "wb");

    // If the file is successfully created start writing
    if (fh != NULL) {

        // The first 128 bits are information about the points (2x64bit integers)
        fwrite (&spaceDimension, sizeof (int64_t), 1, fh);
        fwrite (&numberOfPoints, sizeof (int64_t), 1, fh);

        // Set the range for the random numbers
        double div = RAND_MAX / RANGE;

        // Create the actual points
        for (int j = 0; j < numberOfPoints; j++){
            for(int i = 0; i < spaceDimension; i++){

                double coordinate = -100 + (rand() / div);
                // printf("%.10f ", coordinate);  // Uncomment for debug

                // Write the double to the file
                fwrite(&coordinate, sizeof(double), 1, fh);
            }

            // printf("\n");  // Uncomment for debug
        }
        // Close the file
        fclose (fh);
    }

    // testRead(argv[3]);  // Uncomment for debug

    return 0;
}