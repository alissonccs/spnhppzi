# **spnhppzi: Bayesian Modeling of Recurrent Event Data with Zero Inflation and Spatial Correlation**  

The **spnhppzi** package is designed for **modeling recurrent event data** that exhibit **zero inflation and spatial correlation**. It follows a **Bayesian approach**, implementing **hierarchical models** to capture **complex structures** associated with event repetition and spatial dependencies.  

Spatial modeling is based on the **Intrinsic Conditional Autoregressive (ICAR) model**, which efficiently incorporates spatial correlation. Additionally, the package provides **both parametric and semiparametric models**, where **Bernstein polynomials** are used for **modeling the baseline intensity function**. This approach enhances the model's flexibility, making it applicable in scenarios where **purely parametric models may struggle to adequately capture data complexity**.  

## **Data Simulation**  

The **spnhppzi** package also includes **functions for simulating recurrent event data**, using an adaptation of the **SIMREC** package ([Farrington et al., 2014](https://doi.org/10.18637/jss.v058.i02)). This enables the generation of **spatially correlated recurrent event data**, allowing for **model performance evaluation** and experimentation with different recurrence and spatial dependence scenarios.  

## **References and Application**  

This repository contains the functions used to generate the results presented in the following article:  

📄 **"The Analysis of Criminal Recidivism: A Hierarchical Model-Based Approach for the Analysis of Zero-Inflated, Spatially Correlated Recurrent Events Data"**, available at [arXiv:2405.02666](https://arxiv.org/abs/2405.02666).  

Additional methodological details can be found in the dissertation:  

📖 **"Hierarchical Models for the Analysis of Recurrent Event Data with Zero Inflation and Spatial Correlation"**, available upon request from the author via **alisson.ccs2@gmail.com**.  

## **Features**  

- **Bayesian hierarchical modeling** for recurrent event data;  
- **Spatial correlation modeling** using the **ICAR model**;  
- **Parametric and semiparametric models**, incorporating **Bernstein polynomials** for greater flexibility in modeling the baseline intensity function;  
- **Data simulation** for spatially correlated recurrent events (adaptation of **SIMREC**);  
- **Tools for data preprocessing and model estimation**.  

## **Contributions and Contact**  

Thank you for your interest in **spnhppzi**! If you encounter any issues or have suggestions, feel free to **open an issue** or contact the author via **alisson.ccs2@gmail.com**. 
