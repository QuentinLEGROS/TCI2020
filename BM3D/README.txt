------------------------------------------------------------------

  demo software for iterative Poisson image denoising 
                   Public release v1.0 (16 March 2016) 

------------------------------------------------------------------

Copyright (c) 2016 Tampere University of Technology. 
All rights reserved.
This work should be used only for nonprofit purposes.

Authors:                     Lucio Azzari
                             Alessandro Foi

Web page:                    http://www.cs.tut.fi/~foi/invansc/


------------------------------------------------------------------
Installation and requirements
------------------------------------------------------------------
This software is designed to run on
*) MS Windows, Linux, or Mac OSX (32-bit or 64-bit CPU)
*) Matlab v.7 or later

It requires:
 * Exact unbiased inverse package "invansc" 
    http://www.cs.tut.fi/~foi/invansc/
 * BM3D denoising filter
    http://www.cs.tut.fi/~foi/GCF-BM3D/

The demo script uses also:
 * Statistics Toolbox (to generate Poisson data "poissrnd"),
 * Image Processing Toolbox (only for visualization with "imshow").



------------------------------------------------------------------
Content
------------------------------------------------------------------
demo_iterVSTpoisson.m      demo script
iterVSTpoisson.m           main denoising function
bin_B_h.m                  binning function
debin_Binv_h.m             debinning function
Anscombe_lambda.mat        precomputed expectation vectors and matrix
                             for defining the exact unbiased inverse
                             of noisy+estimate combinations
paramsFromQfun.mat         function handle to set parameters

images_for_table_1 (folder)   test images used for Table 1 of [1]
images_for_table_2 (folder)   test images used for Table 2 of [1]


------------------------------------------------------------------
Reference
------------------------------------------------------------------

[1] L. Azzari and A. Foi, "Variance Stabilization for Noisy+Estimate
Combination in Iterative Poisson Denoising", submitted, March 2016


------------------------------------------------------------------
Disclaimer
------------------------------------------------------------------

Any unauthorized use of these routines for industrial or profit-
oriented activities is expressively prohibited. By downloading 
and/or using any of these files, you implicitly agree to all the 
terms of the TUT limited license:
http://www.cs.tut.fi/~foi/invansc/legal_notice.html

