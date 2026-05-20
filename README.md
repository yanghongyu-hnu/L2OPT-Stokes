# L2OPT-Stokes

This repository contains the official MATLAB implementation of the paper  
**"A non-intrusive model order reduction method based on nonlinear optimization for parameterized Stokes problems"**, published in *Journal of Computational and Applied Mathematics* (JCAM), 2026.

---

## Paper Information

- **Title:**  
  A non-intrusive model order reduction method based on nonlinear optimization for parameterized Stokes problems

- **Authors:**  
  Liang Chen, Qiuqi Li, Hongyu Yang

- **Journal:**  
  Journal of Computational and Applied Mathematics 482 (2026) 117283

- **DOI:**  
  https://doi.org/10.1016/j.cam.2025.117283

---

## Citation

If you use this code in your research, please cite our paper:

```bibtex
@article{chen2026nonintrusive,
  title={A non-intrusive model order reduction method based on nonlinear optimization for parameterized Stokes problems},
  author={Chen, Liang and Li, Qiuqi and Yang, Hongyu},
  journal={Journal of Computational and Applied Mathematics},
  volume={482},
  pages={117283},
  year={2026},
  publisher={Elsevier}
}
```

---

## Method Overview

We propose a **non-intrusive, data-driven reduced-order modeling (ROM) framework** termed **$\mathcal{L}_2$-Opt-ROM** for steady parameterized Stokes problems.

Key features:

1. **Non-intrusive:**  
   Only requires output samples of the full-order model (FOM), no access to internal operators or state variables.

2. **Weighted loss function:**  
   Introduces a weight parameter $\lambda$ to balance approximation errors between velocity and pressure outputs.

3. **Closed-form gradients:**  
   Derives analytical gradients of the objective function with respect to reduced-order matrices using parameter-separable forms.

4. **Block coordinate descent (BCD) algorithm:**  
   Decomposes the non-convex optimization problem into multiple convex subproblems for efficient solution.

5. **Offline-online decomposition:**  
   All computationally expensive optimization is performed offline, enabling rapid online evaluation for new parameters.
