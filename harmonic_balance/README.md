
# Harmonic Balance Module

This directory contains the implementation of the **Harmonic Balance** algorithm, a numerical technique widely used in analyzing and solving nonlinear systems, particularly in the context of radio frequency (RF) and microwave circuit simulations.

## Overview

The harmonic balance method is an efficient frequency-domain technique for solving nonlinear differential equations. It is particularly useful for:

- Analyzing steady-state responses of nonlinear systems.
- Simulating circuits with periodic or quasi-periodic excitations.
- Studying systems where time-domain approaches are computationally expensive.

This implementation is part of the larger **Nomad** project and leverages C++ for performance and flexibility.

## Directory Structure

- **Source Code**: The main source files implementing the harmonic balance technique. These files contain the core logic for solving nonlinear equations using frequency-domain methods.
- **Tests**: Unit tests to validate the correctness and robustness of the implementation.
- **Utilities**: Helper functions and additional utilities used by the harmonic balance solver.

## Key Features

- **Frequency-Domain Analysis**: Utilizes Fourier transformations for efficient computation in the frequency domain.
- **Nonlinear Solver**: Supports solving nonlinear systems using iterative techniques such as Newton-Raphson.
- **Scalability**: Designed to handle large systems with multiple harmonics efficiently.

## How to Use

1. **Building the Module**:
   - Ensure you have the necessary dependencies installed (e.g., C++ compiler, CMake).
   - Follow the build instructions in the [Nomad project README](../README.md).

2. **Running the Solver**:
   - Include the harmonic balance module in your simulation pipeline.
   - Provide the nonlinear system equations and specify the desired number of harmonics for analysis.

3. **Testing**:
   - Run the unit tests to ensure the module is functioning correctly:
     ```bash
     make test
     ```

## Dependencies

This module depends on the following libraries and tools:
- **C++ STL**: Standard Template Library for data structures and algorithms.
- **Eigen**: A C++ template library for linear algebra (if applicable).

## Contributing

We welcome contributions to the harmonic balance module! If you encounter issues or have suggestions for improvements, please open an issue or create a pull request.

## License

This module is part of the **Nomad** project and is licensed under the terms of the project's [LICENSE](../LICENSE).

## Contact

For questions or support, please contact the project maintainer or open an issue in the repository.

---

Thank you for using the Harmonic Balance module! We hope it helps you in your simulations and research.
