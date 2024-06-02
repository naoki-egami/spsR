## spsR: Synthetic Purposive Sampling

------------------------------------------------------------------------

**Description:**

R package `spsR` implements Synthetic Purposive Sampling (SPS) that
allows users to select diverse study sites for external validity, as
developed in [Egami and Lee
(2023+)](https://naokiegami.com/paper/sps.pdf).

**Authors:**

-   [Naoki Egami](https://naokiegami.com) (Maintainer)
-   [Diana Da In Lee](https://www.dianadainlee.com)

**Reference:**

-   [Egami and Lee. (2023+)](https://naokiegami.com/paper/sps.pdf).
    Designing Multi-Context Studies for External Validity: Site
    Selection via Synthetic Purposive Sampling.

### <span style="color:#3366CC">**Overview of the Website**</span>

<p class="lh-lg">
</p>

-   **Please start with [Get Started
    Page](http://naokiegami.com/spsR/articles/intro.html).**

    -   This page should be sufficient for understanding the basic use
        of the package.

<p class="lh-lg">
</p>

-   **We provide more detailed practical guides on [Practical
    Guides](http://naokiegami.com/spsR/articles/intro.html).**

    -   Include (1) step-by-step instructions, (2) answers to frequently
        asked questions, (3) suggestions about data collection.

<p class="lh-lg">
</p>

-   **For methodological details, please read [Egami and Lee
    (2023+)](https://naokiegami.com/paper/sps.pdf).**

<p class="lh-lg">
</p>

-   **We are happy to receive any questions and suggestions about R
    package `spsR`.** <br> Please reach out to Naoki
    (<naoki.egami@columbia.edu>) and Diana (<dl2860@columbia.edu>).

### <span style="color:#3366CC">**Installation Instructions**</span>

<p class="lh-lg">
</p>

You can install the most recent development version using the `devtools`
package. First you have to install `devtools` using the following code.
Note that you only have to do this once:

``` r
if(!require(devtools)) install.packages("devtools")
```

Then, load `devtools` and use the function `install_github()` to install
`spsR`:

``` r
library(devtools)
install_github("naoki-egami/spsR", dependencies = TRUE)
```

### <span style="color:#3366CC">**Acknowledgement**</span>

This research is supported by the National Science Foundation
(SES–2318659).
