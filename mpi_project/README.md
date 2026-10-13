# Parallel Spatial Data Processing with MPI

## Overview

This project implements a parallel spatial data processing program in C using the **Message Passing Interface (MPI)**. 
It distributes a large grd-based spatial dataset across multiple processes, analyzes the local neighborhoods of spatial points, and identifies points that satisfy predefined statistical and density criteria.

## Processing Workflow

1. **Initialize MPI**
   Each process obtains its rank and the total number of processes in `MPI_COMM_WORLD`.

2. **Partition the Input File**
   The input file is divided into byte ranges, with each process assigned a portion of the data.

3. **Read Data in Parallel**
   Processes use collective MPI file reading to load their assigned chunks. Additional bytes are read near chunk boundaries to help handle records spanning multiple chunks.

4. **Parse Spatial Coordinates**
   Text records are parsed into arrays containing X, Y, and Z coordinates.

5. **Analyze Local Neighborhoods**
   For each point assigned to a process, the program calculates distances to the other points in that process's dataset and selects neighbors within a 1,000-meter radius.

6. **Apply Statistical Filters**
   The program calculates the neighborhood's mean and standard deviation and checks whether the point meets the predefined criteria.

7. **Write Results**
   Qualifying points are printed to the console and written to a shared output file using MPI file operations.

## Filtering Criteria

A point is selected when its analyzed neighborhood satisfies all the following conditions:

| Criterion                      | Threshold                   |
| ------------------------------ | --------------------------- |
| Neighborhood radius            | 1,000 meters                |
| Mean Z value                   | Between 5 and 10, inclusive |
| Point coverage                 | At least 85%                |
| Standard deviation of Z values | Less than 1                 |

The point coverage percentage is estimated using 7,853 as the reference number of points, based on the area of a circle with a 1,000-meter radius and an assumed grid-cell area of 400 square meters.

The neighborhood mean and standard deviation are calculated using Welford's algorithm to update the statistics incrementally as neighboring points are processed.

## Execution

Run the program with the desired number of MPI processes:

```bash
mpirun -np 4 ./spatial_analysis data.grd output.txt
```

Arguments:

* `data.grd`: Input file containing spatial coordinate records. (Data can access from <https://data.gov.tw/dataset/176927>)
* `output.txt`: Output file for points satisfying the filtering criteria.
* `-np 4`: Runs the program using four MPI processes; change this value as needed.

Each input record is expected to contain X, Y, and Z values in a grd format supported by the parser.

## Performance Evaluation

The program was tested with different numbers of MPI processes to evaluate the relationship between parallelism and execution time on 9,380,253 data points.


| Number of MPI Processes | Execution Time (min) |
| ----------------------: | -------------------: |
|                       1 |               37,520 |
|                      64 |                  375 |
|                     128 |                  149 |
|                     256 |                   67 |
|                     512 |                   21 |
|                     720 |                   13 |

Co-developed with Chi-Hung (wchiohme@gmail.com), utilizing resources from the Taiwan Central Weather Administration (CWA) to get this result.
