# SURFTest

* This branch allows for the remote accessing of some R code that conducts a MSE-based hypothesis test which is computationally efficient, statistically valid, and powerful. The accompanying paper can be found at https://arxiv.org/pdf/1904.07830.pdf.

* Contents:

    * `MSE_Test_File.R` contains all the relevant functions for conducting the test, including plotting and summary methods. Moreover, it contains information about the proposed permutation test importance measures with additional plotting methods. `source` this file first!
    
    * `MSE_Test_Exs.R` contains all examples of implementations of the code for some toy datasets. This is essentially an abridged version of the package documentation. 
    
    
* Required Packages:

  * `ggplot2`
  * `dplyr`
  * `rpart`
  * `party`
  * `randomForest`
  * `ranger`
  * `glmnet`
  
* Contact: Tim Coleman; tsc35@pitt.edu
