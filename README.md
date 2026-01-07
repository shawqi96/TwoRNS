# TwoRNS

TwoRNS - Two fluid version of RNS, created by Andreas Konstantinou.
Based on Nik Stergioulas' RNS code

Solves for rotating stars in relativity with two non-interacting fluids.

Changes by Shafayat Shawqi
- Dark fluid rotates differentially with zero angular momentum and angular velocity equal
  to the frame-dragging velocity of spacetime caused by the rigidly rotating baryonic fluid.
- Fermionic DM pressure and energy density is calculated with user provided dark matter
  particle mass and self-interaction strength.
- Total mass within the baryonic radius is calculated using the global definition, along with dark cloud mass.

run make,
then 
./trns -e e_cBM -c e_cDM -b eosBM -d f -m mchi -y selfint -r r_ratio

Comments:
- e_cBM, e_cDM in MeV/fm^3
- eosBM is the BM EOS file name
- 'f' means fermionic DM
- mchi is the DM particle mass in MeV
- y is the DM self-interaction strength which is unitless
- r_ratio is the ratio of the polar baryonic radius to the equatorial baryonic radius (0, 1]

