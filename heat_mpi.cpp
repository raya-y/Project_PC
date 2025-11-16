#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Heat equation parameters
#define ALPHA 0.01      // Thermal diffusivity
#define DX 0.01         // Spatial step in x
#define DY 0.01         // Spatial step in y
#define DT 0.0001       // Time step
#define MAX_N 4096      // Maximum grid size

// Global arrays
double grid_old[MAX_N][MAX_N];
double grid_new[MAX_N][MAX_N];

// Performance tracking structure
typedef struct {
	double total_time;
	double comm_time;
	double comp_time;
	double wait_time;
	int messages_sent;
} PerformanceMetrics;

PerformanceMetrics perf_metrics = { 0 };

// Initialize the temperature grid with boundary conditions
void initialize_grid(int N) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {

			// Boundary conditions: hot edges
			if (i == 0 || j == 0 || j == N - 1 || i == N-1) {
				grid_old[i][j] = 100.0;  // Hot boundary
				grid_new[i][j] = 100.0;
			}
			else {
				grid_old[i][j] = 0.0;    // Cold interior
				grid_new[i][j] = 0.0;
			}
		}
	}
}

// finite difference method
void heat_step(int rows, int cols) {
	double rx = ALPHA * DT / (DX * DX);
	double ry = ALPHA * DT / (DY * DY);

	// Update interior points (skip boundaries as if they generate heat)
	for (int i = 1; i < rows - 1; i++) {
		for (int j = 1; j < cols - 1; j++) {
			grid_new[i][j] = grid_old[i][j] +
				rx * (grid_old[i + 1][j] - 2 * grid_old[i][j] + grid_old[i - 1][j]) +
				ry * (grid_old[i][j + 1] - 2 * grid_old[i][j] + grid_old[i][j - 1]);
		}
	}
}

// Swap grids by copying values
void swap_grids(int rows, int cols) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			grid_old[i][j] = grid_new[i][j];
		}
	}
}

// Exchange boundary rows
void exchange_boundaries_nonblocking(int local_rows, int cols, int rank, int size,
	MPI_Comm comm, MPI_Request *requests) {
	double t_start = MPI_Wtime();
	int req_count = 0;

	if (rank > 0) {
		MPI_Isend(&grid_old[1][0], cols, MPI_DOUBLE, rank - 1, 0, comm, &requests[req_count++]);
		MPI_Irecv(&grid_old[0][0], cols, MPI_DOUBLE, rank - 1, 1, comm, &requests[req_count++]);
		perf_metrics.messages_sent++;
	}

	if (rank < size - 1) {
		MPI_Isend(&grid_old[local_rows - 2][0], cols, MPI_DOUBLE, rank + 1, 1, comm, &requests[req_count++]);
		MPI_Irecv(&grid_old[local_rows - 1][0], cols, MPI_DOUBLE, rank + 1, 0, comm, &requests[req_count++]);
		perf_metrics.messages_sent++;
	}

	double t_wait = MPI_Wtime();
	if (req_count > 0) {
		MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
	}
	double t_end = MPI_Wtime();

	perf_metrics.wait_time += (t_end - t_wait);
}

int main(int argc, char *argv[]) {
	// initialize grid size N and time steps   
	int N = 2048;
	int time_steps = 1000;

	// initialize MPI  
	MPI_Init(&argc, &argv);
	// initialize rank and size for each process
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// take input if there is
	if (argc > 1) N = atoi(argv[1]);
	if (argc > 2)time_steps = atoi(argv[2]);


	// Calculate local grid dimensions
	int rows_per_process = N / size;
	int remainder = N % size;
	int local_rows = rows_per_process + (rank < remainder ? 1 : 0); // if there is a remainder add it to some ps
	int global_row_offset = rank * rows_per_process + (rank < remainder ? rank : remainder); // which rows 
	int total_local_rows = local_rows + 2;

	//print main info once 
	if (rank == 0) {
		printf("       2D Heat Equation - MPI Performance Analysis\n");
		printf("-----------------------------------------------------------\n\n");
-		printf("  Grid: %d x %d / Processes: %d / Timesteps: %d\n", N, N, size, time_steps);
		printf("  Rows per process: ~%d / Halo cells: 2 per process\n", rows_per_process);
		// start off by implementing Dirichlet boundary conditions
		initialize_grid(N);
	}

	// wait until all reach here to start a timer
	MPI_Barrier(MPI_COMM_WORLD);
	double start_time = MPI_Wtime();

	// Time stepping loop with detailed performance tracking
	MPI_Request requests[4];

	for (int step = 0; step < time_steps; step++) {
		double t1, t2;

		// Computation phase
		t1 = MPI_Wtime();
		heat_step(total_local_rows, N);
		t2 = MPI_Wtime();
		perf_metrics.comp_time += (t2 - t1);

		// Communication phase
		t1 = MPI_Wtime();
		int req_count = 0;
		if (rank > 0) req_count += 2;
		if (rank < size - 1) req_count += 2;

		exchange_boundaries_nonblocking(total_local_rows, N, rank, size,MPI_COMM_WORLD, requests);
		t2 = MPI_Wtime();
		perf_metrics.comm_time += (t2 - t1);

		swap_grids(total_local_rows, N);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	double end_time = MPI_Wtime();
	perf_metrics.total_time = end_time - start_time;

	// detailed statistics by summing up
	double global_total, global_comm, global_comp, global_wait;
	double max_total, min_total, max_comm, min_comm, max_comp, min_comp;
	double total_bytes_sent, total_messages;
	//SUMS of all metrics 
	MPI_Reduce(&perf_metrics.total_time, &global_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&perf_metrics.comm_time, &global_comm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&perf_metrics.comp_time, &global_comp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&perf_metrics.wait_time, &global_wait, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MIN AND MAX of certain metrics 
	MPI_Reduce(&perf_metrics.total_time, &max_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&perf_metrics.total_time, &min_total, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&perf_metrics.comm_time, &max_comm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&perf_metrics.comm_time, &min_comm, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&perf_metrics.comp_time, &max_comp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&perf_metrics.comp_time, &min_comp, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

	double msg_count = (double)perf_metrics.messages_sent;
	MPI_Reduce(&msg_count, &total_messages, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		//Calculate averages
		double avg_total = global_total / size;
		double avg_comm = global_comm / size;
		double avg_comp = global_comp / size;
		double avg_wait = global_wait / size;
		double halo_time = avg_comm - avg_wait;
		// calculate imbalances
		double load_imbalance = ((max_total - min_total) / avg_total) * 100.0;
		double comm_imbalance = ((max_comm - min_comm) / avg_comm) * 100.0;
		double comp_imbalance = ((max_comp - min_comp) / avg_comp) * 100.0;


		printf("                    PERFORMANCE METRICS\n");
		printf("================================================================\n\n");

		printf("Timing Breakdown:\n");
		printf("----------------------------------------------------------------\n");
		printf("  Total execution time:      %.6f seconds\n", max_total);
		printf("  Average computation time:  %.6f seconds (%.2f%%)\n",
			avg_comp, (avg_comp / avg_total) * 100);
		printf("  Average communication time: %.6f seconds (%.2f%%)\n",
			avg_comm, (avg_comm / avg_total) * 100);
		printf("  Average wait time:         %.6f seconds (%.2f%%)\n",
			avg_wait, (avg_wait / avg_total) * 100);
		printf("\n");

		printf("Load Balance Analysis:\n");
		printf("----------------------------------------------------------------\n");
		printf("  Overall load imbalance:    %.2f%%\n", load_imbalance);
		printf("  Computation imbalance:     %.2f%%\n", comp_imbalance);
		printf("  Communication imbalance:   %.2f%%\n", comm_imbalance);
		printf("  Time range: [%.6f, %.6f] seconds\n", min_total, max_total);
		printf("\n");

		printf("Parallel Efficiency:\n");
		printf("----------------------------------------------------------------\n");
		printf("  Total messages sent:       %.0f messages\n", total_messages);
		printf("  Computation/Communication ratio: %.2f\n", avg_comp / avg_comm);
		printf("  Parallel efficiency estimate:    %.2f%%\n",
			(avg_comp / avg_total) * 100.0);
		printf("  Communication overhead:          %.2f%%\n",
			(avg_comm / avg_total) * 100.0);
		printf("\n");

		// Time breakdown - wait time is part of communication no double count
		double accounted_time = avg_comp + avg_comm;  // Wait is already in comm
		double other_time = avg_total - accounted_time;
		double other_percent = (other_time / avg_total) * 100.0;

		printf("Time Distribution:\n");
		printf("----------------------------------------------------------------\n");
		printf("  Computation:    %.1f%%\n", (avg_comp / avg_total) * 100);
		printf("  Communication:  %.1f%%\n", (avg_comm / avg_total) * 100);
		printf("(halo exchange):  %.1f%%\n", (halo_time / avg_total) * 100);
		printf("    (Wait/Sync):  %.1f%%\n", (avg_wait / avg_total) * 100);
		printf("  Other/Overhead: %.1f%%\n", other_percent);
		printf("\n");

	}

	MPI_Finalize();
	return 0;
}