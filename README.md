# SuccessiveBarriers
This repository contains code for synthesizing Successive Control Barrier Functions. For details refer to our HSCC 2025 [paper](https://dl.acm.org/doi/pdf/10.1145/3716863.3718043)
## Getting started
1. Download and install [Julia](https://julialang.org/). Installation instructions [here](https://docs.julialang.org/en/v1/manual/installation/).
2. Install [Jupyter](https://jupyter.org/).
3. Add required Jilia packages ... ([how?](https://docs.julialang.org/en/v1/stdlib/Pkg/))
   using Pkg;
   Pkg.add(["IJulia", "DynamicPolynomials", "Plots", "SumOfSquares", "CSDP", "LinearAlgebra", "DifferentialEquations", "PyCall", "JuMP", "LaTeXStrings"]
   
4. Now clone this repo
5. Run jupyter notebooks for each case study
6. All plots are saved in the 'figures' folder

```bibtex
    @inproceedings{10.1145/3716863.3718043,
    author = {Wajid, Rameez and Sankaranarayanan, Sriram},
    title = {Successive Control Barrier Functions for Nonlinear Systems},
    year = {2025},
    isbn = {9798400715044},
    publisher = {Association for Computing Machinery},
    address = {New York, NY, USA},
    url = {https://doi.org/10.1145/3716863.3718043},
    doi = {10.1145/3716863.3718043},
    booktitle = {Proceedings of the 28th ACM International Conference on Hybrid Systems: Computation and Control},
    articleno = {21},
    numpages = {11},
    keywords = {Control Barrier Functions, Nonlinear Hybrid Systems, Sum Of Squares Programming},
    location = {Irvine, CA, USA},
    series = {HSCC '25}
    }
```
