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
@inproceedings{wajid2025successive,
  title={Successive Control Barrier Functions for Nonlinear Systems},
  author={Wajid, Rameez and Sankaranarayanan, Sriram},
  booktitle={Proceedings of the 28th ACM International Conference on Hybrid Systems: Computation and Control},
  pages={1--11},
  year={2025}
  location = {Irvine, CA, USA},
  series = {HSCC '25},
  doi = {10.1145/3716863.3718043},
}
```
