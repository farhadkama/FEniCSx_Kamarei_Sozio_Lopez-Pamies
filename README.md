# FEniCSx_Kamarei_Sozio_Lopez-Pamies

! This repository contains FEniCSx example codes for solving boundary value problems of a family of internal-variable-based constitutive models for the viscoelastic response of elastomers of any compressibility undergoing finite deformations, put forth in [1]. A large variety of non-Gaussian elasticities and nonlinear viscosities can be taken into account by using the models introduced by [2]. The 'Uniaxial test' and the 'Single edge notch tension test' [1] are the examples used. These examples serve to illustrate the capability of the model to describe the behavior of viscoelastic elastomers. Google Colab, which provides an easy-to-use interface for running ipynb files, is one way to run these examples.

!**********************************************************************

! Usage: ! ! This code takes as inputs the 5 material properties listed below: ! 1. E = Young's modulus! 2. nu = Poisson's ratio ! 3. Gc = Critical energy release rate ! 4. sts = Tensile strength ! 5. scs = Compressive strength

! In addition, the user must specify the regularization length for these boundary value problems. Typically, this length should be chosen so that it is smaller than the smallest size of the structure, as well as the material characteristic length scale (3Gc)/(16 Wts).

!**********************************************************************

For any inquiry, please contact me at kamarei2@illinois.edu

Alternatively, you may also reach out to my Ph.D. advisor at pamies@illinois.edu

!**********************************************************************

! References:

! [1] Kamarei, F., Sozio, F., Lopez-Pamies, O. 2024. The single edge notch fracture test for viscoelastic elastomers. https://arxiv.org/pdf/2410.15380

! [2] Kumar, A., Lopez-Pamies, O., 2016. On the two-potential constitutive modeling of rubber viscoelastic materials. Comptes Rendus Mécanique 344, 102–112. http://pamies.cee.illinois.edu/Publications_files/CRM_2016.pdf
