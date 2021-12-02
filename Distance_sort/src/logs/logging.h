#ifndef DISTANCE_SORT_LOGGING_H
#define DISTANCE_SORT_LOGGING_H
#include <stdio.h>


typedef struct log {
    char logFileName[1024];  /// the/path/to/the/logFile (includes timestamp)
    FILE *file;  /// The pointer to the file
    int rank;  /// The rank of the process
} Log;

void createFile(Log logFile);

void appendToFile(Log logFile);

void sortLogs(Log logFile);

void printLogs(Log logFile);

#endif
