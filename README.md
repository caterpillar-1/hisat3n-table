# HISAT-3N-Table-NG

[repo](https://github.com/caterpillar-1/hisat3n-table)

This is a rewrite of non-performant `hisat-3n-table` component in  [HISAT-3N](https://github.com/DaehwanKimLab/hisat2/tree/hisat-3n). By employing optimized task partitioning algorithm and The rich ecosystem of the Rust programming language, we implemented the optimization below:

- better task partitioning algorithm for parallel base appending
- mmap-based I/O for economic memory footprint

All planned optimizations in our proposal have been implemented.

## Optimization Strategy

### Task Partitioning Algorithm

With profiling results from Intel VTune. We discovered that the handwritten `SafeQueue` MPMC queue and `class Positions` object pool 

Both the Alignment File and the Reference File are divided into chunks, and the chunks from the two files are paired to form a Task. This changes the task submission granularity from per-line to per-chunk and the result collection granularity from per-position to per-chunk. This reduces thread synchronization overhead and bus contention.

### Mmap-based I/O

In the original implementation, a `class Position` owns its referencing chromosome's name `String`. Consequently, it introduces dramatic unnecessary memory usage.

In our rewritten version, the alignment file and the reference file are Mmapped into the process's memory space as read-only static data (using Rust's `LazyLock` and `Box::leak`). Therefore, a `struct Position` contains simply a reference to the name in reference file in the memory.

## Benchmark result

<!-- TODO -->

## HOWTO

### Install

Clone the repository and use cargo to build `hisat-3n-table` and `dna_index`.

```sh
$ cargo build --bin hisat-3n-table --release
$ cargo build --bin dna_index --release
```

### Run

The command line arguments is almost the same as the original version. Run with `--help` for more details.

## Bug report for the original version

- Hand-written binary search

    ```cpp
    // https://github.com/DaehwanKimLab/hisat2/blob/hisat-3n/position_3n_table.h
    int searchReadNameID (unsigned long long&readNameID, int start, int end) {
        if (uniqueIDs.empty()) {
            return 0;
        }
        if (start <= end) {
            int middle = (start + end) / 2;
            if (uniqueIDs[middle].readNameID == readNameID) {
                return middle;
            }
            if (uniqueIDs[middle].readNameID > readNameID) {
                // ATTENTION HERE!
                // should be searchReadNameID(readNameID, start, middle)
                return searchReadNameID(readNameID, start, middle-1);
            }
            return searchReadNameID(readNameID, middle+1, end);
        }
        return start;
    }
    ```

    We **CANNOT** guarantee mathematical equivilance with the original code, for the bug in it.

    For other parts of the program, we are performing the same operation to the original implementation (appending an alignment to its reference position).

- Use of uninitialized memory

    The fast path to get the location of the alignment line is the 
    ```cpp
    // https://github.com/DaehwanKimLab/hisat2/blob/hisat-3n/hisat_3n_table.cpp#L211`
    bool getSAMChromosomePos(string* line, string& chr, long long int& pos);
    ```
    function. Despite the fast parsing technique, it ignores potential reference location change introduced by `adjustPos` function, which might lead to out-of-bound index in `Positions::refPositions`.

## Acknoledgement

We would like to extend our gratitude to chariri ([cqjjjzr](https://github.com/cqjjjzr)) for his effort in fixing the bugs and rewriting the optimized program.
