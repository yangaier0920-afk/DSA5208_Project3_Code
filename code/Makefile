# Build the CPU programs; run `make` inside the code/ directory.
# The GPU program is built with `bash openmp/build.sh nvidia` (NVHPC).
# Compile on a compute allocation, not on a login node.
CXX      = c++
MPICXX   = mpicxx
OMPFLAGS = -fopenmp
S_FLAGS  = -std=c++17 -O2 -Wall -Wextra -Wpedantic
M_FLAGS  = -std=c++17 -O2 -Wall -Wextra
BUILD    = build

.PHONY: all numerics scaling historical mpi openmp-cpu clean
all: numerics scaling historical mpi openmp-cpu

# Serial solver: convergence/accuracy study and problem-size scaling
numerics: $(BUILD)/jacobi_numerics
scaling:  $(BUILD)/jacobi_a_scaling

# Historical 2D programs used for the serial versus single-rank MPI comparison
historical: $(BUILD)/jacobi_serial $(BUILD)/jacobi_serial_local_update \
            $(BUILD)/jacobi_serial_row_update $(BUILD)/jacobi_mpi

# Cartesian MPI / MPI+OpenMP solver, instrumented variants and two-rank ping-pong
mpi: $(BUILD)/jacobi_cart_mpi $(BUILD)/jacobi_cart_hybrid \
     $(BUILD)/jacobi_cart_tuned $(BUILD)/jacobi_cart_axis $(BUILD)/mpi_link_probe

# OpenMP CPU build (GPU build: bash openmp/build.sh nvidia)
openmp-cpu:
	bash openmp/build.sh cpu

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/jacobi_numerics: serial/jacobi_numerics.cpp | $(BUILD)
	$(CXX) $(S_FLAGS) $< -o $@
$(BUILD)/jacobi_a_scaling: serial/scaling/jacobi_a_scaling.cpp | $(BUILD)
	$(CXX) $(S_FLAGS) $< -o $@

$(BUILD)/jacobi_serial: serial/historical_2d/jacobi.cpp | $(BUILD)
	$(CXX) $(S_FLAGS) $< -o $@
$(BUILD)/jacobi_serial_local_update: serial/historical_2d/jacobi_serial_local_update.cpp | $(BUILD)
	$(CXX) $(S_FLAGS) $< -o $@
$(BUILD)/jacobi_serial_row_update: serial/historical_2d/jacobi_serial_row_update.cpp | $(BUILD)
	$(CXX) $(S_FLAGS) $< -o $@
$(BUILD)/jacobi_mpi: serial/historical_2d/jacobi_mpi.cpp | $(BUILD)
	$(MPICXX) $(S_FLAGS) $< -o $@

$(BUILD)/jacobi_cart_mpi: mpi/jacobi_cart.cpp | $(BUILD)
	$(MPICXX) $(M_FLAGS) $< -o $@
$(BUILD)/jacobi_cart_hybrid: mpi/jacobi_cart.cpp | $(BUILD)
	$(MPICXX) $(M_FLAGS) $(OMPFLAGS) $< -o $@
$(BUILD)/jacobi_cart_tuned: mpi/supplementary/jacobi_cart_tuned.cpp | $(BUILD)
	$(MPICXX) $(M_FLAGS) $(OMPFLAGS) $< -o $@
$(BUILD)/jacobi_cart_axis: mpi/supplementary/jacobi_cart_tuned.cpp | $(BUILD)
	$(MPICXX) $(M_FLAGS) $< -o $@
$(BUILD)/mpi_link_probe: mpi/supplementary/mpi_link_probe.cpp | $(BUILD)
	$(MPICXX) $(M_FLAGS) $< -o $@

clean:
	rm -rf $(BUILD) openmp/build
