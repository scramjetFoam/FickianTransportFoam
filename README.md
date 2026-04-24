# FickianTransportFoam 
**Stable** mixture-averaged and specie-specific constant Lewis transport models that are missing from the current OpenFOAM release. Supported OpenFOAM versions are: 13, 12, 11 and 10.

Developed by Aleksi Rintanen & Ilya Morev, Aalto University, Finland 

The related advanced tutorials are aviable in [DLBFoam-Hydrogen-Tutorials.](https://github.com/Aalto-CFD/DLBFoam-Hydrogen-Tutorials)
## Theory
### Mass diffusion fluxes
FickianTransportFoam implements a transport model based on Fick's law using the formulation, where diffusion fluxes are evaluated with respect to the mass fraction gradient as
```math
\mathbf{j}_k = -\rho D_k \nabla Y_k
```
where $D_k$ is a model specific diffusion coefficient.   
To ensure mass conservation, a correction velocity is introduced to the diffusion fluxes [2] as
```math
 \mathbf{j}_k = -\rho D_k \nabla Y_k + \rho Y_k \sum_{j=1}^{N} D_j \nabla Y_j 
```


Diffusion coefficients are evaluted based on two assumptions:
1. **Mixture-averaged diffusion coefficients**  
   Diffusion coefficients are given by Eq. 12.178 in [2]  
   $$D_k = \left(\sum_{j\neq k} \frac{X_j}{D_{jk}} + \frac{X_k}{1-Y_k}\sum_{j\neq k} \frac{Y_j}{D_{jk}}\right)^{-1}$$  
3. **Constant Lewis number diffusion coefficients**  
   $$D_k = \frac{\lambda}{c_p}\frac{1}{\mathrm{Le}}$$

In LES-models, additional subgrid scale mixing is added using eddy diffusivity concept similarly as in the OpenFOAM's native transport models. Turbulent mass diffusion coefficients are calculated from the turbulent thermal diffusivity using constant turbulent Prantl and Schmidt numbers defined in the model dictionary. 

### Heat fluxes
In the enthalpy equation, heat flux $\mathbf{q}$ is given assuming Fourier's law as
```math
\mathbf{q} = -\lambda \nabla T + \sum_{k=1}^N h_k \mathbf{j}_k
```
Evaluating term $\nabla\cdot\bf q$ is problematic, since the term $-\nabla\cdot \lambda \nabla T$ cannot added directly explicitly for stability reasons.
This library implements two models for evaluating the term $\nabla\cdot\bf q$ in entalphy equation:
1. **Reformulate and add implicitly (recommended)**  
   Reformulate in terms of entalphy using the perfect gas assumption and relation $dh_{s,k}=c_{p,k}dT$  
   $$-\lambda \nabla T = -\frac{\lambda}{c_p} \nabla h_s  + \frac{\lambda}{c_p}\sum_{k=1}^N h_k \nabla Y_k$$
2. **Add $-\nabla\cdot \lambda \nabla T$ explicitly**  
    An additional correction term is included that stabilizes the equation.  
    This approach is used in the native implementation of FickianFourier in OpenFOAM 

> [!NOTE]
> The implicit formulation is not compatible with the default implementation of a coupled temperature boundary condition in OpenFOAM, which balances the heat fluxes using temperature gradients [3]. Thus, for conjugate heat transfer applications with the default coupled temperature boundary condition, the native explicit formulation should be used.
    

## Why OpenFOAM's native implementation "FickianFourier" can be unstable? 
#### 1. Correction velocity is missing
The disparities in diffusion velocities are dumped to the inert specie to ensure mass conservation. This works fine, when the mixture is diluted greatly (for example the inlet is premixed mixture of air and fuel and N2 is thus the dominant specie), but not when there is a separate fuel inlet.

#### 2. Mixture-averaged diffusion coefficients do not behave well
The formula for mixture-averaged diffusion coefficients is undefined when $Y_k\rightarrow 1$ (leads to $\frac{0}{0}$).   
```math
D_k = \frac{1-X_k}{\sum_{j\neq k} X_j/D_{jk}}
```
In FickianFourier this is "fixed" by adding $\epsilon$ to the denominator, which is defined by default as "small" ($\approx10^{-16}$)
```math
D_k = \frac{1-X_k}{\sum_{j\neq k} X_j/D_{jk} + \epsilon}
```
This causes the diffusion coefficient to go zero when $X_k\rightarrow 1$ and negative when $X_k\gt1$. This is both unphysical and numerically unstable leading to divergent simulations.  
(Side note: OpenFOAM uses wrong formulation (Eq. 12.176 in [2]), which is meant for evaluating diffusion flux respect to the molar
average velocity)

In this implementation, this issue is handled by setting $D_k = D_{kk}$, when the mass fraction approaches to 1.

## Usage
Include `libFickianTransport.so` in controlDict and select the desired model in `thermophysicalTransport`.  Available models:

```
Laminar:

mixtureAveraged
constantLewis

LES/RAS:

mixtureAveragedEddyDiffusivity
constantLewisEddyDiffusivity
```

For constant Lewis number models, Lewis numbers are included in the subdictionary as
```
Le {
  H2 1.0;
  O2 1.0;
  N2 1.0;
}
```
A file with Lewis numbers can be generated using [generateLewisNumbersDict.py](utilities/generateLewisNumbersDict.py) utility (use `--help` option to learn more). Note that Lewis numbers have to be selected carefully and default parameters in the script are given for simplicity only.

For binary diffusion coefficients, we included the [implementation of log-polynomial fit](src/ThermophysicalTransportModels/ThermophysicalTransportFunctions/binaryDiffusionCoefficientsPolynomial/binaryDiffusionCoefficientsPolynomial.H), used by Cantera [4]
```math
 D_{jk} = \frac{1}{p} T^{3/2} \sum_{i=0}^{4} a_{jki} \log^i{T}
```
These coefficients can be obtained using [generateBinaryDiffusionCoefficientsPolynomialDict.py](utilities/generateBinaryDiffusionCoefficientsPolynomialDict.py) utility (use `--help` option to learn more).

Example of the dictionary entry in thermophysicalTransport:  

```
laminar
{
    model           mixtureAveraged;
    D
    {
      O-H2 {
        type binaryDiffusionCoefficientsPolynomial;
        polynomialCoeffs (-9.65147826e-03  5.71246016e-03 -1.09656624e-03  9.84367103e-05 -3.22064414e-06);
        }
      O-H {
        type binaryDiffusionCoefficientsPolynomial;
        polynomialCoeffs (-2.64601309e-02  1.37960208e-02 -2.44482561e-03  1.99897138e-04 -5.97357847e-06);
        }
    // ....     
    }
    implicitHeatFlux true; //Default true
    selfDiffusionLimit 1e-6; // Only for OF versions < 13
}
```
Note that since OpenFOAM 13, this implementation [was included](https://github.com/OpenFOAM/OpenFOAM-dev/commit/f9139b7a2cc053f0196fa569f9018405a1dd810b) in the release with a slightly different name and input dictionary format.
 
> [!IMPORTANT]
> On OpenFOAM-dev, both the sensible enthalpy (`sensibleEnthalpy`) and the sensible internal energy (`sensibleInternalEnergy`) formulations are supported for the energy equation. The diffusive heat flux $\sum_k h_{s,k}\mathbf{j}_k$ always uses sensible enthalpy $h_{s,k}$, regardless of the solved variable (this mirrors OpenFOAM-dev's native `Fickian`).
> However, the *implicit* heat-flux formulation is derived from $dh_s = c_p\,dT$ and is only strictly valid when the solved variable is sensible enthalpy. When `sensibleInternalEnergy` is used, the library automatically switches to the explicit formulation (equivalent to `implicitHeatFlux false`) and emits a warning; it is recommended to set `implicitHeatFlux false` explicitly in that case.
> On older releases that only provide `hsi` (e.g. OpenFOAM 10, 11) only the sensible enthalpy formulation is supported.


## Citation

If you use our model, please cite the publication describing its implementation [1].

## References
<a id="1">[1]</a>
A. Haider, I. Morev, A. Rintanen, Z. Shahin, P. Tamadonfar, S. Karimkashi, A. Wehrfritz, V. Vuorinen, Accelerated numerical simulations of hydrogen flames: Open-source implementation of an advanced diffusion model library in OpenFOAM, International Journal of Hydrogen Energy, Volume 189, 152115, [10.1016/j.ijhydene.2025.152115](https://doi.org/10.1016/j.ijhydene.2025.152115) (2025)
<details>
<summary>BibTex</summary>
<p>
 
```
@article{haider2025accelerated,
 author = {Ali Haider and Ilya Morev and Aleksi Rintanen and Zin Shahin and Parsa Tamadonfar and Shervin Karimkashi and Armin Wehrfritz and Ville Vuorinen},
 title = {{Accelerated numerical simulations of hydrogen flames: Open-source implementation of an advanced diffusion model library in OpenFOAM}},
 journal = {International Journal of Hydrogen Energy},
 volume = {189},
 pages = {152115},
 year = {2025},
 issn = {0360-3199},
 doi = {10.1016/j.ijhydene.2025.152115}
}
```
</p>
</details>

[2] Kee, R. J., Coltrin, M. E. & Glarborg, P. Chemically Reacting Flow: Theory and Practice (John Wiley & Sons, 2003).

[3] R. Tuominen, Coupling Serpent and OpenFOAM for neutronics - CFD multi-physics calculations. Master's thesis, Aalto university, Espoo, Helsinki, Aug. 2015

[4] David G. Goodwin, Harry K. Moffat, Ingmar Schoegl, Raymond L. Speth, and Bryan W. Weber. Cantera: An object-oriented software toolkit for chemical kinetics, thermodynamics, and transport processes. https://www.cantera.org, 2023. Version 3.0.0.
