MPICC = mpicc
CC = gcc -std=c99

BUILD_DIR := build
SRC_DIR := src

# Colors
#GREEN = \033[1;32m
#RED = \033[1;31m
#NC = \033[0m

GREEN =
RED =
NC =

# Directories
LIBRARIES_SRC := $(shell find $(SRC_DIR)/data $(SRC_DIR)/processes -name '*.c')
LIBRARIES_SRC := $(LIBRARIES_SRC:%=$(BUILD_DIR)/%.o)

MPI_SRC := $(SRC_DIR)/main.c
MPI_SRC := $(MPI_SRC:%=$(BUILD_DIR)/%.o)


# Every folder in ./src will need to be passed to GCC so that it can find header files
INC_DIRS := $(shell find $(SRC_DIR) -type d)

# Add a prefix to INC_DIRS. So moduleA would become -ImoduleA. GCC understands this -I flag
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CFLAGS := $(INC_FLAGS) -O3 -g
MPI_FLAGS := --mca opal_warn_on_missing_libcuda 0 -hostfile hostfile

# If the first argument is "run_mpi"...
ifeq (run_mpi,$(firstword $(MAKECMDGOALS)))
    # use the rest as arguments for "run_mpi"
    RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
    # ...and turn them into do-nothing targets
    $(eval $(RUN_ARGS):;@:)
else ifeq (debug_mpi,$(firstword $(MAKECMDGOALS)))
    # use the rest as arguments for "debug_mpi"
    RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
    # ...and turn them into do-nothing targets
    $(eval $(RUN_ARGS):;@:)
endif

pre_build:
	@echo
	@echo "$(GREEN)Building MPI binary ...$(NC)"
	@echo

build_mpi: pre_build $(MPI_SRC) $(LIBRARIES_SRC)
	@echo
	@$(MPICC) $(CFLAGS) -o $(BUILD_DIR)/mpi_main.out $(MPI_SRC) $(LIBRARIES_SRC) -lpthread -lm
	@echo "$(GREEN)Build finished successfully!$(NC)"
	@echo

$(BUILD_DIR)/%.c.o: %.c
	@echo "    Building" $@ "..."
	@mkdir -p $(dir $@)
	@$(MPICC) $(CFLAGS) -c $< -o $@

%.c:

run_mpi:
	@mpirun $(MPI_FLAGS) $(BUILD_DIR)/mpi_main.out $(RUN_ARGS)

debug_mpi:
	@mpirun $(MPI_FLAGS) valgrind $(BUILD_DIR)/mpi_main.out $(RUN_ARGS)


create_datafile:
	@mkdir -p datasets
	@echo
	@echo "$(GREEN)Creating data file with d=20 n=8388608 path=./datasets/data.bin$(NC)"
	@$(CC) -o $(BUILD_DIR)/create_sample_points.out $(SRC_DIR)/tools/create_sample_points.c
	@$(BUILD_DIR)/create_sample_points.out 20 8388608 ./datasets/data.bin
	@echo

.PHONY: clean
clean:
	@echo
	@echo "$(RED)Clearing build directory...$(NC)"
	@rm -rf $(BUILD_DIR)/*
	@echo "$(GREEN)All done!$(NC)"
	@echo