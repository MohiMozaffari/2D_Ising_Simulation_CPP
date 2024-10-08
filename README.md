# 2D Ising Simulation

This project simulates the 2D Ising model using Monte Carlo methods and calculates various thermodynamic properties such as energy, magnetization, specific heat, magnetic susceptibility, and Binder cumulant.

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Simulation Parameters](#simulation-parameters)

## Introduction

The Ising model is a mathematical model used in statistical mechanics to simulate the behavior of ferromagnetic materials. This simulation uses Monte Carlo methods to simulate the behavior of a 2D lattice of spins interacting with each other.

## Installation

To run the simulation, you will need a C++ compiler that supports C++11 or later. 

1. Clone the repository:

   ```bash
   git clone https://github.com/MohiMozaffari/2D_Ising_Simulation_CPP.git 
   ```

2. Compile the code:

    ```bash
    g++ -o Ising_2D/main.cpp
    ```

## Usage
To run the simulation, execute the compiled binary:
```bash
./main
```

The simulation will generate a CSV file containing the results of the simulation.

## Simulation Parameters
* Lattice Size (N): 16
* Number of Ensembles (nens): 25
* Number of Ensembles in Critical Region (nens_crit): 50
* Minimum Temperature (tmin): 2
* Maximum Temperature (tmax): 3.0
* Upper Critical Temperature (tcrit_up): 2.3
* Lower Critical Temperature (tcrit_down): 2.2
* Temperature Step Size (deltat): 0.2
* Temperature Step Size Near Critical Region (deltat_crit): 0.01
* Number of Monte Carlo Steps per Temperature (steps): N^3

  


