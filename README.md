# Tutorial: Statistical System Identification Methods for Subsynchronous Forced Oscillations

This repository contains the MATLAB code and data accompanying the paper, _A Tutorial on Identification of Subsynchronous Mode Frequencies in Power Transmission Systems using Parametric and Non-parametric Methods_, presented at the [2023 Texas Power and Energy Conference](https://tpec.engr.tamu.edu/).

The tutorial encompasses two cases: a single mode (unimodal) forced oscillation (case 1), and a multi-modal forced oscillation, labeled as case 2. Each of the `*.m` files uses a different observable among current, frequency, and voltage. The two different statistical techniques are:

- *Non-parametric method*: Welch-based power spectral density (PSD) estimation;
- *Parametric method*: output-based autoregressive (AR) model.

## Compatibility

The code in this repository was developed using MATLAB R2021b in Windows.

## Citing

If you find our work useful, we kindly ask you cite the following paper:

```
@INPROCEEDINGS{Dorado-Rojas2023-ds,
  title     = "{A Tutorial on Identification of Subsynchronous Mode Frequencies
               in Power Transmission Systems using Parametric and
               Non-parametric Methods}",
  booktitle = "{2023 Texas Power and Energy Conference}",
  author    = "Dorado-Rojas, Sergio A and Wang, Zongjie",
  month     =  feb,
  year      =  2023
}
```

## License

This repository is [licensed](./LICENSE.txt) under the [GNU Affero General Public License (GNU AGPL)](https://www.gnu.org/licenses/agpl-3.0.en.html).

## Acknowledgements

This work was funded in part by the Colombian Ministry of Science, Technology, and Innovation (Minciencias) under grant agreement number 885-2020, and in part by the Eversource Energy Center project entitled _A Pathway to Enable Sustainable Modern Power Systems: Optimal System Dispatch_.
