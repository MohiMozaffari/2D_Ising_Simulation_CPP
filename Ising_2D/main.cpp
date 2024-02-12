#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

// Function to initialize the lattice
void initializer(vector<vector<double>>& lattice, double agg, int N) {
    lattice.assign(N, vector<double>(N, agg));
}

// Function to calculate the Hamiltonian at a given position
double hamiltonian(int r, int c, const vector<vector<double>>& lattice, int N) {
    double E = -lattice[r][c] * (lattice[(r + 1) % N][c] + lattice[(r - 1 + N) % N][c] +
                                lattice[r][(c - 1 + N) % N] + lattice[r][(c + 1) % N]);
    return E;
}

// Function to flip the spin at a given position if the energy change is favorable
void check_flip(int r, int c, vector<vector<double>>& lattice, int N, double T,  mt19937& gen, uniform_real_distribution<double>& ran_u) {
    double dE = -2 * hamiltonian(r, c, lattice, N);
    if (dE <= 0 || ran_u(gen) < exp(-dE / T)) {
        lattice[r][c] *= -1;
    }
}


void MonteCarloSteps(const int numsteps, vector<vector<double>>& lattice, int N, double T, mt19937& gen, uniform_real_distribution<double>& ran_u, uniform_int_distribution<int>& ran_pos){
    for (int i = 0; i < numsteps; i++) {
        int row = ran_pos(gen);
        int col = ran_pos(gen);
        double dE = -2 * hamiltonian(row, col, lattice, N);
        if (dE <= 0 || ran_u(gen) < exp(-dE / T)) {
            lattice[row][col] *= -1;
        }
    }
}

// Function to calculate the magnetization of the lattice
double Magnetization(const vector<vector<double>>& lattice, int N) {
    double mag = 0.0;
    for (const auto& row : lattice) {
        for (double spin : row) {
            mag += spin;
        }
    }
    return abs(mag) / (N * N);
}

// Function to calculate the energy of the lattice
double Energy(const vector<vector<double>>& lattice, int N) {
    double energy = 0.0;
    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            energy += hamiltonian(r, c, lattice, N);
        }
    }
    return energy / (2 * N * N);
}


// Function to calculate the mean of a vector
double calculateMean(const std::vector<double>& data) {
    double sum = 0.0;
    for (double value : data) {
        sum += value;
    }
    return sum / data.size();
}


// Function to calculate the specific heat capacity
double calc_spec_heat(const std::vector<double>& E_array, double T) {
    double sum_E = 0.0;
    double sum_E2 = 0.0;
    
    // Calculate the sum of E and E^2
    for (double E : E_array) {
        sum_E += E;
        sum_E2 += E * E;
    }
    
    // Calculate the average of E and E^2
    double avE = sum_E / E_array.size();
    double avE2 = sum_E2 / E_array.size();
    
    // Calculate the specific heat capacity
    double avE_2 = avE * avE;
    double C = (1 / T) * (1 / T) * (avE2 - avE_2);
    
    return C;
}

// Function to calculate the magnetic susceptibility
double calc_mag_sus(const std::vector<double>& M_array, double T) {
    double sum_M = 0.0;
    double sum_M2 = 0.0;
    
    // Calculate the sum of M and M^2
    for (double M : M_array) {
        sum_M += M;
        sum_M2 += M * M;
    }
    
    // Calculate the average of M and M^2
    double avM = sum_M / M_array.size();
    double avM2 = sum_M2 / M_array.size();
    
    // Calculate the magnetic susceptibility
    double avM_2 = avM * avM;
    double X = (1 / T) * (avM2 - avM_2);
    
    return X;
}


// Function to calculate the Binder cumulant
double calculateBinderCumulant(const std::vector<double>& data) {
    double cum2 = 0.0;
    for (double value : data) {
        cum2 += (value ) * (value );
    }
    double variance = cum2 / data.size();

    double sum = 0.0;
    for (double value : data) {
        sum += pow(value , 4);
    }
    double fourthMoment = sum / data.size();
    return 1.0 - (fourthMoment / (3.0 * variance * variance));
}


// Function to write vectors to a CSV file
void writeVectorsToCSV(const vector<double>& temp, const vector<double>& Mag, const vector<double>& Energy, const vector<double>& C_v, const vector<double>& Chi, const vector<double>& binder, const string& filename) {
    ofstream file(filename);

    if (file.is_open()) {
        // Write column headers
        file << "Temperature,Magnetization,Energy,Specific Heat,Magnetic Susceptibility,Binder Cumulant\n";

        // Write data to the file
        size_t size = temp.size();
        for (size_t i = 0; i < size; ++i) {
            file << temp[i] << "," << Mag[i] << "," << Energy[i] << "," << C_v[i] << "," << Chi[i]<< "," << binder[i]<< "\n";
        }

        // Close the file
        file.close();
        cout << "File saved successfully!" << endl;
    }
    else {
        cout << "Error opening the file!" << endl;
    }
}



int main() {
    // Define simulation parameters
    const int N = 16; // Lattice size
    const int nens = 25; // Number of ensembles
    const int nens_crit = 50; // Number of ensembles in critical region
    const double tmin = 2; // Minimum temperature
    const double tmax = 3.0; // Maximum temperature
    const double tcrit_up = 2.3; // Upper critical temperature
    const double tcrit_down = 2.2; // Lower critical temperature
    const double deltat = 0.2; // Temperature step size
    const double deltat_crit = 0.01; // Temperature step size near critical region
    const int steps = pow(N,3); // Number of Monte Carlo steps per temperature

    mt19937 gen(random_device{}()); // Random number generator
    uniform_int_distribution<int> ran_pos(0, N-1); // Uniform distribution for selecting random lattice positions
    uniform_real_distribution<double> ran_u(0.0, 1.0); // Uniform distribution for generating random numbers

    // Define vectors to store simulation results
    vector<double> energy;
    vector<double> magnetization;
    vector<double> specificHeat;
    vector<double> magneticSusceptibility;
    vector<double> binderCumulant;
    vector<double> temperature;


    auto start = high_resolution_clock::now();

    vector<vector<double>> lattice(N, vector<double>(N));
    initializer(lattice, 1.0, N);

    //tmin --> tcrit_down
    for (double T = tmax ; T > tcrit_up; T -= deltat) {
        temperature.push_back(T);
        vector<double> energyTemp, magnetizationTemp;

        for (int ens = 0; ens < nens; ens++) {
            MonteCarloSteps(steps, lattice, N, T, gen, ran_u, ran_pos);

            energyTemp.push_back(Energy(lattice, N));
            magnetizationTemp.push_back(Magnetization(lattice, N));
        }
        energy.push_back(calculateMean(energyTemp));
        magnetization.push_back(calculateMean(magnetizationTemp));
        specificHeat.push_back(calc_spec_heat(energyTemp, T));
        magneticSusceptibility.push_back(calc_mag_sus(magnetizationTemp, T));
        binderCumulant.push_back(calculateBinderCumulant(magnetizationTemp));
    }

    //tcrit_down --> tcrit_up
    for (double T = tcrit_up ; T > tcrit_down ; T -= deltat_crit) {
        temperature.push_back(T);
        vector<double> energyTemp, magnetizationTemp;

        for (int ens = 0; ens < nens_crit; ens++) {

            MonteCarloSteps(steps, lattice, N, T, gen, ran_u, ran_pos);

            energyTemp.push_back(Energy(lattice, N));
            magnetizationTemp.push_back(Magnetization(lattice, N));
        }
        energy.push_back(calculateMean(energyTemp));
        magnetization.push_back(calculateMean(magnetizationTemp));
        specificHeat.push_back(calc_spec_heat(energyTemp, T));
        magneticSusceptibility.push_back(calc_mag_sus(magnetizationTemp, T));
        binderCumulant.push_back(calculateBinderCumulant(magnetizationTemp));
    }    
    
    //tcrit_up --> tmax
    for (double T = tcrit_down ; T > tmin; T -= deltat) {
        temperature.push_back(T);
        vector<double> energyTemp, magnetizationTemp;

        for (int ens = 0; ens < nens; ens++) {

            MonteCarloSteps(steps, lattice, N, T, gen, ran_u, ran_pos);

            energyTemp.push_back(Energy(lattice, N));
            magnetizationTemp.push_back(Magnetization(lattice, N));
        }

        energy.push_back(calculateMean(energyTemp));
        magnetization.push_back(calculateMean(magnetizationTemp));
        specificHeat.push_back(calc_spec_heat(energyTemp, T));
        magneticSusceptibility.push_back(calc_mag_sus(magnetizationTemp, T));
        binderCumulant.push_back(calculateBinderCumulant(magnetizationTemp));
    }
    
    string filePrefix = "C:\\Users\\Asus\\Documents\\Programming\\ISING_2D\\IsingII_";
    string fileExtension = ".csv";
    string fileName = filePrefix + to_string(N) + fileExtension;

    writeVectorsToCSV(temperature, magnetization, energy, specificHeat, magneticSusceptibility, binderCumulant, fileName);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);

    cout << "Simulation completed in " << duration.count() << " seconds." << endl;

    return 0;
}