# Radar-Based 3D Trajectory Estimation and Stability Analysis

This project contains MATLAB code for reconstructing the 3D trajectory of a moving person using time-of-flight data from radar sensors and for analyzing the sensitivity of the numerical method to noise.

## Overview

The script performs:

1. **3D position reconstruction** using LU decomposition.
2. **Trajectory visualization** in 3D and 2D over time.
3. **Error propagation analysis** with simulated multiplicative noise on distance measurements.
4. **Perturbation analysis** of radar sensor positions and its influence on solution stability, determinant, and condition number of the system matrix.

---

## Data

The file `MNUB_24L_P1_dane17.mat` contains:

- `D` — matrix of distances from radar sensors
- `R` — matrix of radar sensor positions
- `t` — time vector

Due to privacy concerns, the `.mat` file contains only synthetic or anonymized data.

---

## How to Use

1. Open `MATLAB`.
2. Load the script into the editor.
3. Run the file using:
   ```matlab
   run('nazwa_pliku.m')

Make sure the MNUB_24L_P1_dane17.mat file is in the same directory or update the path in the script:

load('MNUB_24L_P1_dane17.mat');

---

## Tasks
Task 1: Matrix Formulation
Construct matrix A and right-hand side b using distance data.

Set up a linear system Ax = b for position estimation.

Task 2: LU Decomposition
Solve the system using LU decomposition for 165 time points.

Visualize the trajectory in 3D and coordinate values over time.

Task 3: Noise Sensitivity
Apply multiplicative Gaussian noise to D.

Compare actual vs. estimated trajectory errors.

Use condition number cond(A) to estimate error bounds.

Task 4: Sensor Perturbation Analysis
Independently perturb the z-coordinate of each radar sensor.

Analyze how small errors in sensor positioning affect the solution, determinant, and condition number of matrix A.

---

## Plots Generated
3D trajectory of the moving person alongside radar sensor positions.

2D plots of x(t), y(t), and z(t).

Log-log plot of error vs. standard deviation sigma.

Plots showing how sensor position uncertainty affects:

Trajectory reconstruction error

Determinant of matrix A

Condition number of matrix A

---

## Requirements
MATLAB (tested with R2021 or newer)

No additional toolboxes required

---

## License
This project is intended for educational use only.
