#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

long parse_text_to_arrays(const char *chunk, int locstart, int locend, long fileLineNum, long *x, long *y, float *z);
void analyze_spatial_points(int rank, const long *x, const long *y, const float *z, long total_points, MPI_File *out);

void parprocess(MPI_File *in, MPI_File *out, const int rank, const int size, const int overlap);

int main(int argc, char *argv[]) {
    MPI_File in, out;
    int size, rank, ierr;
    const int overlap = 100; // Predefined padding byte size to safeguard broken lines at edges


    MPI_Init(&argc, &argv); // Initialize the global MPI execution framework
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Query and store total amount of available active tasks
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Query and store the unique numeric rank of current process

    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_filename> <output_filename>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    // Open a file across all processes in a given communicator group
    ierr = MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    if (ierr) {
        if (rank == 0) fprintf(stderr, "Error: Cannot open input file %s\n", argv[1]);
        MPI_Finalize();
        exit(2);
    }

    ierr = MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &out);
    if (ierr) {
        MPI_File_close(&in);
        if (rank == 0) fprintf(stderr, "Error: Cannot open output file %s\n", argv[2]);
        MPI_Finalize();
        exit(3);
    }

    // Execute Main Parallel Analysis
    parprocess(&in, &out, rank, size, overlap);

    MPI_File_close(&in);
    MPI_File_close(&out);

    // Exit the MPI framework
    MPI_Finalize();
    return 0;
}

void parprocess(MPI_File *in, MPI_File *out, const int rank, const int size, const int overlap) {
    MPI_Offset filesize = 0, globalstart = 0, globalend = 0;
    int mysize = 0;
    char *chunk = NULL;

    // Get current file size for each process
    MPI_File_get_size(*in, &filesize);
    if (filesize > 0) filesize--;

    // *** Important part : Calculate the start and end offsets for each process ***
    mysize = filesize / size;
    globalstart = rank * mysize;
    globalend = (rank == size - 1) ? (filesize - 1) : (globalstart + mysize - 1);

    // Read extra bytes to avoid incomplete records
    if (rank != size - 1) globalend += overlap;
    mysize = globalend - globalstart + 1;

    if (mysize <= 0) return;

    // Allocate buffer
    chunk = (char *)malloc((mysize + 1) * sizeof(char));
    if (!chunk) return;

    // All processes must call this collective function together to prevent deadlocks
    MPI_File_read_at_all(*in, globalstart, chunk, mysize, MPI_CHAR, MPI_STATUS_IGNORE);
    chunk[mysize] = '\0';

    // Refine boundary pointers to start and end
    int locstart = 0, locend = mysize - 1;
    if (rank != 0) {
        char prev_char = '\0';
        MPI_File_read_at(*in, globalstart-1, &prev_char, 1, MPI_CHAR, MPI_STATUS_IGNORE);

        if (prev_char != '\n') {
            while (locstart < mysize && chunk[locstart] != '\n') {
                locstart++;
            }
            locstart++;
        }
    }

    if (rank != size - 1) {
        locend -= overlap;
        if (locend < 0) locend = 0;
        while (locend < mysize && chunk[locend] != '\n') {
            locend++;
        }
    }

    if (locstart > locend) {
        free(chunk);
        return;
    }

    // Pre-calculate total line count in this processor's allocated chunk
    long fileLineNum = 0;
    for (int i = locstart; i < locend; i++) {
        if (chunk[i] == '\n') fileLineNum++;
    }
    fileLineNum++;

    // Dynamically allocate target arrays for storing point spatial coordinates
    long *x = (long *)malloc(fileLineNum * sizeof(long));
    long *y = (long *)malloc(fileLineNum * sizeof(long));
    float *z = (float *)malloc(fileLineNum * sizeof(float));

    if (!x || !y || !z) {
        if (x) free(x); if (y) free(y); if (z) free(z);
        free(chunk);
        return;
    }

    // Parse string data to arrays
    long valid_points = parse_text_to_arrays(chunk, locstart, locend, fileLineNum, x, y, z);
    free(chunk);
    chunk = NULL;

    /*printf("\n=== [DEBUG Rank %d] Array Date (Total Number : %ld) ===\n", rank, valid_points);
    for (long p = 0; p < valid_points; p++) {
        printf("Index %ld -> X: %ld, Y: %ld, Z: %.2f\n", p, x[p], y[p], z[p]);
    }*/

    // Run spatial proximity matrix filter and output statistics logs
    analyze_spatial_points(rank, x, y, z, valid_points, out);

    // Free dynamically allocated coordinate arrays
    free(x);
    free(y);
    free(z);
}

long parse_text_to_arrays(const char *chunk, int locstart, int locend,
    long fileLineNum, long *x, long *y, float *z) {

    char line[512] = {0};
    char delim[] = " ,-\r\n";
    int charcount = 0;
    long xyzcount = 0;

    for (int i = locstart; i <= locend; i++) {
        if (charcount < 511) {
            line[charcount++] = chunk[i];
        }

        if (chunk[i] == '\n' || i == locend - 1) {
            line[charcount] = '\0';
            charcount = 0;

            char *ptr = strtok(line, delim);
            if (ptr != NULL) {
                x[xyzcount] = atol(ptr);
                int idx = 1;
                while ((ptr = strtok(NULL, delim)) != NULL) {
                    if (idx == 1)      y[xyzcount] = atol(ptr);
                    else if (idx == 2) z[xyzcount] = (float)atof(ptr);
                    idx++;
                }
                xyzcount++;
                if (xyzcount >= fileLineNum) break;
            }
        }
    }
    return xyzcount;
}

void analyze_spatial_points(int rank, const long *x, const long *y, const float *z, long total_points, MPI_File *out) {

    char io_buffer[512];
    for (int i = 0; i < total_points; i++) {
        int count = 0;

        double mean = 0.0;
        double M2 = 0.0;

        for (int j = 0; j < total_points; j++) {
            float dist = sqrt(pow(x[j] - x[i], 2) + pow(y[j] - y[i], 2));
            if (dist > 0 && dist <= 1000) {
                count++;

                float val = z[j];
                double delta = val - mean;
                mean += delta/count;
                double delta2 = val - mean;
                M2 += delta * delta2;
            }
        }

        if (count > 0) {
            // Check if the dense neighborhood meets the 85% coverage threshold.
            // 7853 is the ideal total point count within a 1000m radius circle,
            // calculated as: Circle Area (pi * 1000^2) divided by resolution/grid area (400).
            //                -> (pi * 10^6) / 400

            float arr_percent  = (float)count / 7853.0f * 100.0f;
            float arr_mean = (float)mean;
            float arr_variance = (float)(M2 / count);
            if (arr_variance < 0) arr_variance = 0;
            float arr_SD = (float)sqrt(arr_variance);

            if ((arr_mean >= 5 && arr_mean <= 10) && (arr_percent >= 85.f) && (arr_SD < 1)) {
                printf("[Rank %d] | Coordination:(%ld, %ld) | Z:%.2f | Mean:%.2f | SD:%.2f\n",
                       rank, x[i], y[i], z[i], arr_mean, arr_SD);

                int msg_len = sprintf(io_buffer, "[Rank %d] | Coordination:(%ld, %ld) | Z:%.2f | Mean:%.2f | SD:%.2f\n",
                                      rank, x[i], y[i], z[i], arr_mean, arr_SD);

                MPI_File_write_shared(*out, io_buffer, msg_len, MPI_CHAR, MPI_STATUS_IGNORE);
            }
        }
    }

}