# FEniCSx_Kamarei_Sozio_Lopez-Pamies

! This repository contains FEniCSx example codes for solving boundary value problems of a family of internal-variable-based constitutive models for the viscoelastic response of elastomers of any compressibility undergoing finite deformations, put forth in [1]. A large variety of non-Gaussian elasticities and nonlinear viscosities can be taken into account by using the models introduced by [2]. The 'Uniaxial test' and the 'Single edge notch tension test' [1] are the examples used. These examples serve to illustrate the capability of the model to describe the behavior of viscoelastic elastomers. Google Colab, which provides an easy-to-use interface for running ipynb files, is one way to run these examples.

!**********************************************************************

! Usage: ! ! This code takes as inputs the 16 material properties listed below.

! For Elasticity:

! 1. mu1 = $$\mu_1$$ ! 2. mu2 = $$\mu_2$$ ! 3. mu_pr = $$\kappa$$ ! 4. alph1 = $$\alpha_1$$ ! 5. alph2 = $$\alpha_2$$ ! 6. m1 = $$\nu_1$$ ! 7. m2 = $$\nu_2$$ ! 8. a1 = $$\beta_1$$ ! 9. a2 = $$\beta_2$$

! For Viscosity:

! 10. K1 = $$K_1$$ ! 11. K2 = $$K_2$$ ! 12. bta1 = $$\gamma_1$$ ! 13. bta1 = $$\gamma_2$$ ! 14. eta0 = $$\eta_0$$ ! 15. etaInf = $$\eta_{\inf}$$ ! 14. N0 = $$N_0$$

! In addition, the user must specify the geometry or dicretization of interest.

!**********************************************************************

For any inquiry, please contact me at kamarei2@illinois.edu

Alternatively, you may also reach out to my Ph.D. advisor at pamies@illinois.edu

!**********************************************************************

! References:

! [1] Kamarei, F., Sozio, F., Lopez-Pamies, O. 2024. The single edge notch fracture test for viscoelastic elastomers. https://arxiv.org/pdf/2410.15380

! [2] Kumar, A., Lopez-Pamies, O., 2016. On the two-potential constitutive modeling of rubber viscoelastic materials. Comptes Rendus Mécanique 344, 102–112. http://pamies.cee.illinois.edu/Publications_files/CRM_2016.pdf
