MPICC = mpicc


BUILD_DIR := ./build
SRC_DIR := ./src

# Colors
#GREEN = \033[1;32m
#RED = \033[1;31m
#NC = \033[0m

GREEN =
RED =
NC =

# Directories
LIBRARIES_SRC := $(shell find $(SRC_DIR)/data $(SRC_DIR)/logs -name '*.c')
LIBRARIES_SRC := $(LIBRARIES_SRC:%=$(BUILD_DIR)/%.o)

MPI_SRC := $(SRC_DIR)/main.c
MPI_SRC := $(MPI_SRC:%=$(BUILD_DIR)/%.o)


# Every folder in ./src will need to be passed to GCC so that it can find header files
INC_DIRS := $(shell find $(SRC_DIR) -type d)
# Add a prefix to INC_DIRS. So moduleA would become -ImoduleA. GCC understands this -I flag
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CFLAGS := $(INC_FLAGS) -O3
MPI_FLAGS := --mca opal_warn_on_missing_libcuda 0 -hostfile hostfile


all: build_mpi  # pre_build build_mpi post_build


build_mpi: $(MPI_SRC) $(LIBRARIES_SRC)
	@echo "$(GREEN)Building MPI binary...$(NC)"
	@$(MPICC) $(CFLAGS) -o $(BUILD_DIR)/mpi_main.out $(MPI_SRC) $(LIBRARIES_SRC)
	@echo "$(GREEN)Build finished successfully!$(NC)"
	@echo


$(BUILD_DIR)/%.c.o: %.c
	@mkdir -p $(dir $@)
	@$(MPICC) $(CFLAGS) -c $< -o $@

run_mpi: build_mpi
ifneq ("$(wildcard $(data_file))","")
	@mpirun $(MPI_FLAGS) $(BUILD_DIR)/mpi_main.out $(data_file)
else
	@echo "$(RED)Missing data_file={argument} or file does not exist$(NC)"
	@echo "Run with 'make data_file=path/to/file run_mpi'"
endif

.PHONY: clean
clean:
	@echo "$(RED)Clearing build directory...$(NC)"
	@rm -rf $(BUILD_DIR)
	@echo "$(GREEN)Done!$(NC)"