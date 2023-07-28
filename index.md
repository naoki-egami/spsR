## sps: Synthetic Purposive Sampling

------------------------------------------------------------------------

**Description:**

R package `sps` implements the SPS algorithm for site selection and the
SPS estimator.

**Authors:**

-   [Naoki Egami](https://naokiegami.com) (Maintainer)
-   [Diana Da In Lee](https://www.dianadainlee.com)

**Reference:**

-   Egami and Lee. (2023+). Designing Multi-Context Studies for External
    Validity: Site Selection via Synthetic Purposive Sampling.

### Installation Instructions

You can install the most recent development version using the `devtools`
package. First you have to install `devtools` using the following code.
Note that you only have to do this once:

``` r
if(!require(devtools)) install.packages("spsR")
```

Then, load `devtools` and use the function `install_github()` to install
`spsR`:

``` r
library(devtools)
install_github("naoki-egami/spsR", dependencies = TRUE)
```
