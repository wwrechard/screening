# screening

An R package for efficient variable **screening** for linear and generalized linear models

## Description
This package implements four popular variable screening algorithms, the **sure independence screening**<sup>[1](#myfootnote1)</sup> (SIS), **high-dimensional ordinary least squares projection**<sup>[2](#myfootnote2)</sup> (HOLP), **rank-robust correlation screening**<sup>[3](#myfootnote3)</sup> (RRCS) and the **Forward regression**<sup>[4](#myfootnote4)</sup> (FR) for efficient screening operaion over large p small n models. **SIS**, **HOLP** and **FR** are implemented to work with both linear models and generalized linear models while **RRCS** works only for linear models. We refer the reader to the references and the R-help docs for detailed information.

## Dependency
For linear model screening, the package does not rely on any external packages. For generalized linear model, **HOLP** relies on the package `glmnet` to compute the ridge-regression estimator. This part can be replaced by any user function or package that achieves a similar functionality.

## Installation
To install the package in a painless way, you will need to first install `devtools` package and then execute the following lines in your favorate R IDE
```
> library(devtools)
> install_github('wwrechard/screening')
> library(screening)
```
The package is now ready for use. Simply call the helper by typing
```
> ?screening
```
will provide further information on the usage of the function. If there is any issue installing the `devtools` package, you may instead downlowd the R files from directory `/R` directly. The file `screening.r` contains all working functions required for using the algorithm.

## Quick start
A simple working example. Assume we observe 200 samples with 500 dimension, stored in `x` and `y`. Following code illustrates screening using a linear model
```
> library(screening)
> output = screening(x, y, method = 'holp', num.select = 50)
> print (output$screen)
```
`output$screen` contains the indexes of the screened variables, which, in this case, is a 50-dimensional vector. If you opt to use extended BIC for determining the screening set, then simply replace the second line in the above code with
```
> output = screening(x, y, method = 'holp', ebic = TRUE)
```
The package contains other two testing functions `linearModelTest` and `logisticTest`. The two functions provide integrated tests on the four screening methods. The full functionality includes simulating synthetic datasets and apply screening methods and report the result. Their codes also serve as good examples of using `screening`. For more detailed information, please refer to their R-help docs.

### Author
[Xiangyu (Samuel) Wang](https://github.com/wwrechard). I appreciate your feedback under [*issues*](https://github.com/wwrechard/screening/issues) of this repository. 

### References
<a name="myfootnote1">1. </a>Fan, Jianqing, and Jinchi Lv. "Sure independence screening for ultrahigh dimensional feature space." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70.5 (2008): 849-911.

<a name="myfootnote2">2. </a>Wang, Xiangyu, and Chenlei Leng. "High-dimensional ordinary least-squares projection for screening variables." Journal of the Royal Statistical Society: Series B (Statistical Methodology) (To appear)

<a name="myfootnote3">3. </a>Li, Gaorong, et al. "Robust rank correlation based screening." The Annals of Statistics 40.3 (2012): 1846-1877.

<a name="myfootnote4">4. </a>Wang, Hansheng. "Forward regression for ultra-high dimensional variable screening." Journal of the American Statistical Association 104.488 (2009): 1512-1524.
