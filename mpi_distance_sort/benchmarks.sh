#!/bin/bash

# Welcome message
function welcome () {
    echo
    echo " __  __  _____  _____     ____                      _                              _"
    echo "|  \/  ||  __ \|_   _|   |  _ \                    | |                            | |"
    echo "| \  / || |__) | | |     | |_) |  ___  _ __    ___ | |__   _ __ ___    __ _  _ __ | | __ ___"
    echo "| |\/| ||  ___/  | |     |  _ <  / _ \| '_ \  / __|| '_ \ | '_ \` _ \  / _\` || '__|| |/ // __|"
    echo "| |  | || |     _| |_    | |_) ||  __/| | | || (__ | | | || | | | | || (_| || |   |   < \__ \\"
    echo "|_|  |_||_|    |_____|   |____/  \___||_| |_| \___||_| |_||_| |_| |_| \__,_||_|   |_|\_\|___/"
    echo
    echo
}


# Help message
function usage() {
    echo "Usage $0 [OPTION]..."
    echo
    echo "Runs a series of tests with the MPI executable and the data files found in the directory specified"
    echo
    echo "OPTIONS:"
    echo "    -d {directory}              The directory containing the datasets."
    echo "                                This is mandatory and the data files must"
    echo "                                have .bin extension"
    echo "    -v true || false || one     Verbose output for the executable. Values:"
    echo "                                    true: All runs are verbose "
    echo "                                    false: All runs are non verbose"
    echo "                                    one: Only the first run is verbose (Default)"
    echo "    -n {processes}              The number of processes. Must be a power of 2!! Default is 8."
    echo "    -h                          Help message for the executable"
    echo "    -u                          Usage for this script"
    echo
    echo "Example usage:"
    echo "    $0 -v -d ./datasets"
    echo
}

# Checks if a number is power of 2
function isPowerOfTwo() {
    declare -i n=$1
    (( n > 0 && (n & (n - 1)) == 0 ))
}

# Variables
DIR=""
VERBOSITY="v"
VERBOSITY_FLAG="true"
HELP=""
PROCESSES=8

# parse the flags
while getopts ":v:n:d:hu" options; do
    case "${options}" in
    d)
        # Set the directory
        DIR=${OPTARG}
        ;;
    v)
        # If -v is true VERBOSITY will be v else it will be nv
        if [ "${OPTARG}" == "true" ]; then
            VERBOSITY="v"
            VERBOSITY_FLAG="false"
        elif [ "${OPTARG}" == "false" ]; then
            VERBOSITY="nv"
            VERBOSITY_FLAG="false"
        elif [ "${OPTARG}" == "one" ]; then
            VERBOSITY="v"
            VERBOSITY_FLAG="true"
        fi
        ;;
    h)
        # Set help flag
        HELP="help"
        ;;
    n)
        # check if number is power of 2 and if not display message and exit with error
        if isPowerOfTwo "${OPTARG}"; then
            PROCESSES=${OPTARG}
        else
            echo
            echo "Number of processes is not a power of 2!"
            echo
            echo
            usage
            exit 1
        fi
        ;;

    u)
        # Display help message and exit gracefully
        usage
        exit 0
        ;;
    :)
        # handles all the cases that an argument was needed and one was not provided
        echo "Error -${OPTARG} requires an argument."
        usage
        exit 1
        ;;
    ?)
        # Display help message and exit with error
        usage
        exit 1
        ;;
    esac
done

# If non of DIR or HELP were set exit with error
if [ "${DIR}" == "" ] && [ "${HELP}" == "" ]; then
    echo "Expected -d DIRECTORY or -h and got none."
    usage
    exit 1
fi

# DIR and HELP are mutually exclusive. So if they are both set DIR has the advantage and we unset the HELP
if [ "${DIR}" != "" ] && [ "${HELP}" != "" ]; then
    HELP=""
fi

# Write the hostfile file with the number of processes
echo "127.0.0.1 slots=${PROCESSES}" > hostfile

# Display welcome message
welcome

# Build the executable from scratch
make clean
make build_mpi
echo
echo

# Run all the benchmarks
for i in "${DIR}"/*.bin; do
  echo "Running for $i"
  make run_mpi "${i}" "${VERBOSITY}"
  echo
  echo
  echo
  echo
  if [ "${VERBOSITY_FLAG}" == "true" ]; then
      VERBOSITY="nv"
  fi
done




