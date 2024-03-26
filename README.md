# vcPB: Varying Coefficient Peters-Belson Method for Longitudinal Data

[![R-CMD-check](https://github.com/SangkyuStat/vcPB/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SangkyuStat/vcPB/actions/workflows/R-CMD-check.yaml)

A package for estimating the disparity between a majority group and minority group based on the extended model of the Peters-Belson method. Our model is the first extension of Peters-Belson method to the longitudinal data. 

Furthermore, our method can set a modifiable variable to find the complicated association between other variables and the modifiable variable.

### Installation

The current version can be installed from source using the package `devtools`
```R
devtools::install_github("SangkyuStat/vcPB")
```

It can also be found on CRAN
```R
install.packages("vcPB")
```


### Usage Examples

#### - `vc.pb` function 
`vc.pb` function provides three different types of models based on the different input arguments: `modifier` and time varying coefficients. 

If `modifier` is `NULL` (the default setting is `NULL`) and at least a time-varying variable exists, then the simple varying-coefficient Peters-Belson method using a gaussian kernel regression can be performed as below:
```R
vc.pb(formula = response ~ (time varying variable | time variable) + 
                variable, 
                id,
                data = input_data, 
                group = disparity_group)
```
If `modifier` is not `NULL` and is a discrete variable, and at least a time-varying variable exists, then the modifiable varying-coefficient Peters-Belson method using a gaussian kernel regression can be performed as below:
```R
vc.pb(formula = response ~ (time varying variable | time variable) + 
                variable + discrete modifier, 
                id,
                data = input_data, 
                group = disparity_group, 
                modifier = "discrete modifier")
```
If `modifier` is not `NULL` and is a continuous variable, and at least a time-varying variable exists, then the simple varying-coefficient Peters-Belson method using a gaussian kernel regression can be performed as below:
```R
vc.pb(formula = response ~ (time varying variable | time variable) + 
                variable + continuous modifier, 
                id,
                data = input_data, 
                group = disparity_group, 
                modifier = "continuous modifier")
```
The type of modifier returns the different results. If there are more than one time-varying variables, the user can perform the function as below:
```R
vc.pb(formula = response ~ (time varying variable1 | time variable) + 
                (time varying variable2 | time variable) + variable,
                id,
                data = input_data, 
                group = disparity_group)
```
<!--the user has to indicate whether the variable is time-varying or not. If there is no time-varying variable, then user can perform the function as below:
```R
vc.pb(formula = response ~ variable + 
                any modifier, 
                id,
                data = input_data, 
                group = disparity_group, 
                modifier = "any modifier")
```-->
If there is no modifier and time-varying variable, then the model is just the naive PB model. For this case, the user can use `pb` function instead.

The user needs to define `group` properly to measure the disparity between two groups in `group` variable, there should be 2 levels for this variable. 

The user needs to define `id` properly to have the exact identification on observations whether they are measured repeatedly across the time. 

The selection of bandwidths is essential and important for the kernel regression. If there is nothing given as initial values, we get and use the default marginal bandwidth from the function `KernSmooth::dpill`. For all models, `bandwidth_M`, `bandwidth_m`, `bandwidth_xM` and `bandwidth_xm` are essential. If `modifier` is not `NULL` and is a continuous variable, then `bandwidth_Z_M`, `bandwidth_Z_m`, `bandwidth_Z_xM` and `bandwidth_Z_xm` are needed more.

Also, use needs to specify local time points (`local_time`) for the time-varying kernel regression. The function will automatically give the time points if there is nothing given. The local time points will be returned in the fitted object.

#### - `pb` function
`pb` function provides the original Peters-Belson method of Peters (1941) and Belson (1956). The usage is as similar as the `vc.pb` but the user should not put the time varying coefficients and a modifier variable.

### Developing

- The conditional version will be uploaded very soon.
- The cross-validation function for choosing the bandwidths will be developed.
- We are trying to develop other methods as well.

### References

Peters, C. C. (1941) A method of matching groups for experiment with no loss of population. Journal of Educational Research, 34, 606-612.

Belson, W. A. (1956) A Technique for Studying the Effects of a Television 
Broadcast.  JRSSC, 5(3), 195-202.
