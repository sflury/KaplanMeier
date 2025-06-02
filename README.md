# KaplanMeier

Statistical assessments with the Kaplan-Meier survival function (lower/upper limits) to test whether a measured value `x0` (typically the mean of a distribution) is associated with some population `x`, accounting for lower limits in `x`. If the Kaplan-Meier survival function at `x0` is outside the limits of 0.01 to 0.99, one can confidently reject the null hypothesis that the measurement `x0` is associated with the measurements `x`.

I originally developed and implemented this script for Flury et al. 2024
for tests involving a sample of 89 Lyman continuum measurements,
39 of which were upper limits requiring the censoring treatment of
the Kaplan-Meier survival curve.

## Examples
``` python
import matplotlib.pyplot as plt
from numpy.random import seed,rand,randn
from KaplanMeier import *

seed(123)
x = randn(100)
c = rand(100)<0.3
x0 = array([1.65])
x0_err = array([[0.3],[0.5]])

km_x,km_y = km_curve(x,c)
p_x,p_e = km_eval(x0,x,c,x0_err=x0_err)
```
which gives the results below

![image of Kaplan-Meier curve with test measurement](km_examp.png "example Kaplan-Meier test")


## BibTeX reference

While this code is provided publicly, I request that any use 
thereof be cited in any publications in which this code is used.
BibTeX formatted reference provided below.

``` bibtex
@ARTICLE{Flury2024,
       author = {{Flury}, Sophia R. and {Jaskot}, Anne E. and {the LzLCS Collaboration}},
        title = "{The Low-Redshift Lyman Continuum Survey: The Roles of Stellar Feedback and ISM Geometry in LyC Escape}",
      journal = {\apjs},
     keywords = {Reionization, Galactic and extragalactic astronomy, Ultraviolet astronomy, Hubble Space Telescope, 1383, 563, 1736, 761, Astrophysics - Astrophysics of Galaxies, Astrophysics - Cosmology and Nongalactic Astrophysics},
         year = 2024,
        month = {},
       volume = {},
       number = {},
          eid = {},
        pages = {},
          doi = {10.48550/arXiv.2409.12118},
archivePrefix = {arXiv},
       eprint = {2409.12118},
          url = {https://github.com/sflury/KaplanMeier},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System} }
```

An additional reference to consider is the original Kaplan & Meier (1958) paper 
which first presented the Kaplan-Meier statistic. The BibTeX entry for their
paper is listed below.

``` bibtex
@article{KaplanMeier1958,
author = {E. L. Kaplan and Paul Meier},
title = {Nonparametric Estimation from Incomplete Observations},
journal = {Journal of the American Statistical Association},
volume = {53},
number = {282},
pages = {457-481},
year = {1958},
publisher = {Taylor & Francis},
doi = {10.1080/01621459.1958.10501452},
}
```

## DOI

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11406486.svg)](https://doi.org/10.5281/zenodo.11406486)

## Licensing
<a href="https://github.com/sflury/KaplanMeier">KaplanMeier</a> © 2024 by <a href="https://sflury.github.io">Sophia Flury</a> is licensed under <a href="https://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International</a>

<img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/nc.svg" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/nd.svg" style="max-width: 1em;max-height:1em;margin-left: .2em;">

This license enables reusers to copy and distribute KaplanMeier in any medium or format in unadapted form only, for noncommercial purposes only, and only so long as attribution is given to the creator. CC BY-NC-ND 4.0 includes the following elements:

BY: credit must be given to the creator.
NC: Only noncommercial uses of the work are permitted.
ND: No derivatives or adaptations of the work are permitted.

You should have received a copy of the CC BY-NC-ND 4.0 along with this program. If not, see <https://creativecommons.org/licenses/by-nc-nd/4.0/>.
