**README: Finite Element Method Implementation for Structural Analysis**

**Project Overview:**
This repository contains the implementation of a Finite Element Method (FEM) code in MATLAB for solving a 2D linear elastic structural problem. The project aims to familiarize users with the steps involved in commercial software packages like Abaqus and Ansys, offering a cost-effective alternative for simple structural problems.

**Problem Description (Problem 7.1):**
The structural problem involves a 9.99 mm by 5 mm rectangular structure represented by 5 quad elements. It includes fixed support at the left bottom edge, rollers at the left top and right bottom ends, and uniform forces applied to the top and right side edges.

**Key Features:**
1. **FEM Code Development:**
   - Implemented in MATLAB to replicate the algorithms used in commercial software.
   - Covers shape function interpolation, mapping, numerical integration, matrix assembly, and applying boundary conditions.

2. **Mesh Generation:**
   - Utilized four-noded quadrilateral plane-strain elements to create a mesh representing the given geometry.
   - Read nodal coordinates and element connectivity from input data.

3. **Boundary Conditions and Forces:**
   - Applied roller boundary conditions on the left side and bottom.
   - Implemented a uniform force of 10Kn/mm on the entire top side and 2KN/mm for 2.5mm of the right side edge.

4. **Material Properties:**
   - Incorporated an elastic modulus of 1 GPa and Poisson’s ratio of 0.3 for accurate structural simulation.

5. **Validation with Abaqus:**
   - Cross-validated results by independently solving the problem in Abaqus.
   - Addressed any discrepancies and provided explanations for variations.

6. **Comprehensive Analysis:**
   - Calculated nodal displacements, forces at imposed boundary nodes, and stress components (sig_11, sig_22) at specific integration points.

**Usage Instructions:**
1. Clone the repository.
   ```bash
   https://github.com/subhodeepbakshi/FEM_Structural_Analysis_MATLAB.git
   ```
2. Run the MATLAB script to execute the FEM code.
3. Input necessary parameters and data as per your problem requirements.

**Dependencies:**
- MATLAB (version X.X or higher)

**Acknowledgments:**
This project was completed as part of an assignment. Special thanks to the professor for the guidance and inspiration.

Feel free to explore the code, experiment with different parameters, and contribute to the improvement of the FEM code. If you encounter any issues or have suggestions, please open an issue, and I'll be happy to assist!
