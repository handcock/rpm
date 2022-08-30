<img src="man/figures/rpm_hl.png" align="left" width="250" height="250" alt="RDS network"/>

# The rpm package
This is an R package to estimate revealed preferences based on observed bipartite matchings.

It was originally developed by Ryan Admiraal and Mark S. Handcock.

## Installation
To install it, you can also use:
```
# If devtools is not installed:
# install.packages("devtools")

devtools::install_github("eelawson/rpm", auth_token="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")
```
## Resources

To run an example use:
```
library(rpm)
data(fauxmatching)
fit <- rpm(~match("edu") + WtoM_diff("edu",3),
          Xdata=fauxmatching$Xdata, Zdata=fauxmatching$Zdata,
          X_w="X_w", Z_w="Z_w",
          pair_w="pair_w", pair_id="pair_id", Xid="pid", Zid="pid",
          sampled="sampled")
summary(fit)
```

For details on how to construct data for input to `rpm()` see the documentation:
```
help(fauxmatching)
```

For information on the current terms that can be used in formulas for `rpm()` see the documentation:
```
help("rpm-terms")
```

<!-- badges: start -->
[![R-CMD-check](https://github.com/eelawson/rpm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/eelawson/rpm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

