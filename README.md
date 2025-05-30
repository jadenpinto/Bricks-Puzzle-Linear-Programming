# Brick Wall Puzzle Solver

This project uses the GLPK (GNU Linear Programming Kit) library in C to solve the Brick Wall Puzzle by formulating it as a linear programming problem. It models the puzzle as a Maximum Flow application known as the Project Selection Problem (Maximum Weight Closure).

## Prerequisites

- [Cygwin](https://www.cygwin.com/) terminal
- GCC (GNU Compiler Collection), and GLPK library installed within Cygwin

## Building the Project

1. In the Cygwin terminal, navigate to the directory of C file:

   ```bash
   cd /path/to/brick_wall_puzzle.c
   ```
2. Compile the source code using gcc:
   ```bash
   gcc -o puzzle_out brick_wall_puzzle.c -lglpk
   ```

## Running the Code

To run the solver on a set of test input files:
```bash
./puzzle_out Tests/*
```

Use the `-d` flag to run the code in debug mode:
```bash
./puzzle_out -d Tests/in01.txt
```
