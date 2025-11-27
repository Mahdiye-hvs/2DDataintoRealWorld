# **3D Printing of Two-dimensional Data for Chemical Education**

Mahdiye Hassanpoor, James J. Harynuk, and Ryland T. Giebelhaus*  
Department of Chemistry, University of Victoria, Victoria, Canada  
rgiebelhaus@uvic.ca  

---

## 1.0 About  
This repository contains a MATLAB workflow and application designed to convert two-dimensional datasets into 3D-printable STL models. Many analytical techniques like GC×GC and 2D NMR produce data that are typically visualized as 2D contour or color plots. Although helpful, traditional 2D plots can make it difficult to fully interpret the shape, depth, and spatial relationships present in the data.

To address this, we developed an easy-to-use workflow that transforms a 2D matrix into a physical, tactile 3D model. The process allows users to import raw data, apply optional smoothing, and generate an STL file suitable for most commercial 3D printers. The MATLAB app included in this repository also supports adding 3D text (such as titles and axis range labels) directly onto the model, so printed objects retain contextual information.

In addition to the MATLAB-based workflow, we provide a compiled standalone version of the app so that users can run the tool without installing MATLAB, using only the free MATLAB Runtime.

This tool is intended for educational use, accessibility support, and general data visualization in chemical research and teaching environments.

---

## 2.0 Use  
Users have two main options:

1. **MATLAB workflow and app** – for users who have access to MATLAB.  
2. **Standalone application** – for users without MATLAB, using the compiled installer and MATLAB Runtime.

For the MATLAB version, download or clone the entire repository and add the `src` folder to your MATLAB path. Users may either run the provided functions scripts in the `src/` directory or launch the interactive MATLAB app (`STLapp.mlapp`) located in `Release`.

The basic workflow consists of:  
1. Loading or importing a 2D matrix representing the chemical data  
2. (Optional) Applying smoothing to the surface  
3. Adding text labels for title and axis limits, if desired  
4. Generating a corresponding STL file for 3D printing  

The tool accepts any 2D matrix (Z) and, optionally, X and Y coordinate grids.

Currently in v1.0 for general use.

---

## 2.1 Inputs  

**Z:**  
A 2D matrix representing the surface or signal intensity (such as chromatographic or spectroscopic data).

**X, Y (optional):**  
Coordinate vectors or grids corresponding to the first and second dimensions. If not provided, indices of the matrix are used.

**smoothFactor:**  
User-defined scalar controlling the degree of smoothing applied to the data surface.

**titleText:**  
Text string to include as a 3D title on the model.

**xRangeText / yRangeText:**  
Numbers indicating the x- and y-axis range limits.

**outputFile:**  
Name of the STL file to be generated.

---

## 2.2 Outputs of the App  

The workflow outputs an STL file suitable for slicing and 3D printing. Depending on the user’s settings, the STL may include additional features such as:

- **3D surface model:** Representation of the input 2D matrix as a 3D object.  
- **Base:** An optional flat base under the model.  
- **3D text and numbers:** Title and axis range values attached directly to the model.  

---

## 3.0 Repository Structure  

```text
src/                Core MATLAB functions for STL generation
src/app/            Interactive MATLAB App (STLapp.mlapp)
examples/           Example scripts and demo data
figures/            Optional images for documentation
LICENSE             MIT License
CITATION.cff        Citation metadata
