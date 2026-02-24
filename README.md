# **TeraFlux**

**Thermodynamically Enabled and Reaction Attuned FLUXome estimation**

TeraFlux is a Python library designed to predict metabolic fluxes by integrating gene expression data (FPKM) with Genome-Scale Metabolic Models (GEMs). It utilizes the principle of maximum entropy production, formulated as a non-linear optimization problem and solved efficiently using CasADi and IPOPT. This repository contains the code and data to reproduce the analyses for the manuscript *"Thermodynamically Enabled and Reaction Attuned Estimation of Metabolic Fluxes"* submitted to *iScience*.

## **🚀 Key Features**

* **Entropy-Based Flux Estimation:** Employs a Shannon entropy term within the objective function to align metabolic fluxes with transcriptomic evidence.  
* **Objective Function:** Maximizes the entropy term ![][image1] while satisfying stoichiometric and thermodynamic constraints.  
* **High Performance:** Built with CasADi for symbolic differentiation and IPOPT for robust non-linear optimization.  
* **GPR Integration:** Includes a built-in Gene-Protein-Reaction (GPR) parser that handles isozymes (OR) and enzyme complexes (AND) to map gene expression to reactions.  
* **Standardized Workflow:** Automated model updating, medium application, and FBA-based initial guess (![][image2]) for faster solver convergence.

## **🛠 Installation**

To use the teraflux.py library, ensure you have the required dependencies installed:

pip install numpy pandas cobra casadi psutil

*Note: The solver requires IPOPT to be available within the CasADi installation.*

## **📂 Repository Structure**

The repository is organized to ensure the complete reproducibility of the manuscript figures (Figures 1 through 7\) and contains the necessary inputs for all evaluated organisms.

* teraflux.py: The core library and command-line execution script.  
* fig1/ through fig7/: Analysis directories containing Jupyter notebooks (.ipynb) and shell scripts (.sh). For computationally intensive tasks, .sh scripts generate the underlying data, which is subsequently read and processed by the notebooks to create the final visualizations.  
* figs2\_3/: Contains the specific core and genome-scale comparisons for *E. coli* (ecoli\_core and ecoli\_iJO1366).  
* **Organism Data Directories** (fig5/Ecoli/, fig5/Bsubtilis/, fig5/Scerevisiae/, fig5/Sstipitis/):  
  * gems/: SBML metabolic models (e.g., iJO1366.xml).  
  * transcriptomes/: Processed FPKM gene expression data files.  
  * mediums/: Culture medium definitions listing active exchange reactions.  
  * knownFluxes/: Experimental fluxes used for constraint or validation.  
  * results/: Output directories containing the computed fluxes and Lagrange multipliers, separated by method (e.g., teraflux/, pheflux/, teraflux\_allReversible/).

## **💻 Usage**

### **Command Line Execution**

You can run TeraFlux directly from the terminal by providing a tab-separated input file that lists the experimental conditions you want to process.

python teraflux.py inputData\_Ecoli\_glc.csv \-o results/teraflux \-p ecoli\_run

### **1\. Main Input File (inputData\_\*.csv)**

This configuration file coordinates the paths to the necessary biological and computational data. It must be a tab-separated file with the following columns:

| Organism | Condition | GeneExpFile | Medium | Network | KnownFluxes |
| :---- | :---- | :---- | :---- | :---- | :---- |
| Ecoli | glc | ./transcriptomes/Ecoli/Ecoli\_Expfile\_glc.csv | ./mediums/Ecoli\_Medium\_glc.csv | ./gems/iJO1366.xml | ./knownFluxes/Ecoli/Ecoli\_knownFluxes\_glc.csv |

### **2\. Medium Specification File (\*\_Medium\_\*.csv)**

The medium file identifies which exchange reactions allow the uptake of metabolites. TeraFlux reads this file and sets the lower bounds of these specified exchange reactions to a negative value to allow influx.

| Metabolite\_Name | Reaction\_ID |
| :---- | :---- |
| ca2 | EX\_ca2\_e |
| cbl1 | EX\_cbl1\_e |
| cl | EX\_cl\_e |
| cobalt2 | EX\_cobalt2\_e |
| cu2 | EX\_cu2\_e |

## **🧪 Example Analysis Pipeline**

The Jupyter notebooks provided in the repository (such as runTeraflux\_iJO1366.ipynb) demonstrate the complete end-to-end pipeline used for the manuscript:

1. **Input Generation:** Programmatically creating the inputData configuration files for various carbon sources (e.g., glucose, acetate, galactose).  
2. **Model Customization:** Re-writing SBML models to open bounds (e.g., setting all internal reactions to be reversible between \-1000 and 1000\) to explore the unconstrained metabolic solution space.  
3. **Execution:** Batch processing all conditions using the getFluxes function from the teraflux.py library.  
4. **Output Parsing:** Reading the resulting .fluxes.csv, .forward\_reverse\_fluxes.csv, and .Lagrange.csv files to generate summary tables and performance metrics.

## **🧬 Biological Context**

TeraFlux bridges the gap between static metabolic networks and dynamic transcriptomic states. By applying an entropy principle, it provides a "Reaction Attuned" estimation of the fluxome, offering a thermodynamically consistent way to study metabolic shifts under varying environmental and genetic conditions.

*Developed by Marcelo Rivas Astroza and the Metabolic Modeling team at Universidad Tecnológica Metropolitana (UTEM).*

[image1]: <data:image/png;base64,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>

[image2]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABMAAAAYCAYAAAAYl8YPAAABfElEQVR4XqVTPUsEMRC9rVQQsVlsdjMJt5XNgZWFWAgigvb2h6WlFiLWFjZ2WtvY2Oqv0EIQe8HC4kAsLTzfZBOTdT+SxQdhZ+a9NzuZvRsMKkiqaReC0qDAISBto/16EqVqLxlopp0OoGasFQza6kFEGCMkJaKFjfiHW0p5IYTY55jbID8jopM/sm4opZakkk8co9kBGnzhTDknQQ9E4rzq6IA1MtI0nedcSTWSUq1yjAn3un5jGpaEYcXWMNmR3xyNZm1s8i3wbziToihmfK4GiD79Zj5gXgD3YvMWnbsCC3CuK7QB6tgfjb2cVyB/r494EUue5nk+tPvCVZedQTzrICnN4HYdRxOcQ5tzsysWZVk2h+ejmUyxG195BPONZ2Zu28vfwV/anJEYEb91M8+zoXD5qScrJyOxYytUTnbsNLFI9Jd+perOvnE2fFkUeMm49jrMd7bGk3qS/kCDW5x7rOIDk671/P/2Ehv09QT1AYGmA5pGNHuaqxY/b41Pgt35vJcAAAAASUVORK5CYII=>
