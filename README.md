# Pinching-Antenna Systems (PASS): Power Radiation Model and Optimal Beamforming Design
This is code repository for our recent paper "Pinching-Antenna Systems (PASS): Power Radiation Model and Optimal Beamforming Design" accepted by IEEE Transactions on Communications (TCOM) [(https://ieeexplore.ieee.org/document/11263923)](https://ieeexplore.ieee.org/document/11263923).

A pinching antenna system (PASS) relying on practical discrete activation is considered. We propose a novel adjustable power radiation model. By tuning the waveguide-antenna spacing to adjust coupling strength, it can achieve flexible or equal power radiation, without the need to modify the fabricated lengths of pinching antennas (PAs).

Under the commonly assumed equal power radiation, the transmit beamforming, pinching beamforming (determined by activated locations of PAs), and the numbers of activated PAs along each waveguide are jointly optimzed. 

We propose **the first globally optimal** joint pinching beamforming and digital beamforming algorithm. Then, a **low-complexity** matching-based algorithm is also developed, which has demonstrated **near-optimal** performance compared to the globally optimum.

## Algorithm Implementation
- **Globally optimal beamforming algorithms for PASS**: The proposed BnB algorithms are implemented in [BnB_SU.m](BnB_SU.m) for single-user scenario and [BnB.m](BnB.m) for multi-user scenario.
- **Near-optimal beamforming algorithms for PASS based on many-to-many matching theory**:The proposed welfare-driven many-to-many matching algorithms are implemented in [Matching.m](Matching.m) for multi-user scenario and [Matching_SU.m](Matching_SU.m) for single-user scenario.

## Code Reproduction Guideline
- Test transmit power versus number of pre-mounted PAs  $L$ per-waveguide
  - Run [main_L.m](main_L.m)
- Test transmit power under different spacial ranges $S_{x}$
  - Run [main_S.m](main_S.m)
  - For comparison with continuous architecture, further run [main_continuous.m](main_continuous.m)
- Test transmit power under differnet minimum SINR requirements $\gamma_{\min}$
  - Single-user scenario: [main_SU_gamma_min.m](main_SU_gamma_min.m)
  - Multi-user scenario: [main_MU_gamma_min.m](main_MU_gamma_min.m)
  > Since these scripts perform the BnB algorithm, it takes about 5-8 hours to complete the calculations.
- Plot figures for the obtained results using [plot_MU.m](plot_MU.m) for multi-user scenarios and [plot_SU.m](plot_SU.m) for single-user scenarios

## Future Works
- Power radiation optimization
- Discretization of power radiation

## Reference
If you find this code useful for your research, please consider citing 
> X. Xu, X. Mu, Z. Wang, Y. Liu, and A. Nallanathan, ``Pinching-Antenna Systems (PASS): Power Radiation Model and Optimal Beamforming Design'', *IEEE Trans. Commun.*, early access, 2026.
