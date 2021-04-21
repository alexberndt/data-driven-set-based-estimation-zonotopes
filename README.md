# Data-Driven Set-Based Estimation using Matrix Zonotopes with Set Containment Guarantees

by Alexander Berndt, Amr Alanwar, Karl Henrik Johansson and Henrik Sandberg

Manuscript: [Data-Driven Set-Based Estimation using Matrix Zonotopes with Set Containment Guarantees](https://arxiv.org/abs/2101.10784) 

## Abstract

We propose a method to perform set-based state estimation of an unknown dynamical system using a data-driven set propagation function. Our method comes with set-containment guarantees, making it applicable to the estimation of safety-critical systems. The method consists of two phases: (1) an offline learning phase where we collect noisy state-input data to determine a function to propagate the state-set ahead in time; and (2) an online estimation phase consisting of a time update and a measurement update. It is assumed that sets bound measurement noise and disturbances, but we assume no knowledge of their statistical properties. These sets are described using zonotopes, allowing efficient propagation and intersection operations. We propose two approaches to perform the measurement update. The method is extended to constrained zonotopes. Simulations show that the proposed estimator yields state sets comparable in volume to the confidence bounds obtained by a Kalman filter approach, but with the addition of state set-containment guarantees. We observe that using constrained zonotopes yields smaller sets, but with higher computational cost compared to unconstrained zonotopes. 

## Gettting Started

### Dependencies

The code requires the installation of the [CORA MATLAB toolbox](https://tumcps.github.io/CORA/)

### File Descriptions

The main code used throughout the [manuscript](https://arxiv.org/abs/2101.10784) can be found in
```
estimator_simulation.m
```
The comparison of singular values as in Fig. 4 of the [manuscript](https://arxiv.org/abs/2101.10784) can be found in 
```
singular_value_comparison.m
```

## Reference 

If you find this code and/or our paper useful, consider citing it using the following bibtex entry

```
@article{berndt2021setbasedestimation,
  title={{Data-Driven Set-Based Estimation using Matrix Zonotopes with Set Containment Guarantees}},
  author={Berndt, Alexander and Alanwar, Amr and Johansson, Karl Henrik and Sandberg, Henrik}
  journal={arXiv preprint arXiv:2101.10784},
  year={2021}
}
```
