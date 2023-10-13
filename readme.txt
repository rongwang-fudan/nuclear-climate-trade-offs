All files are for a paper that is under review.

Please use these files with caution, as the paper is still under review.

Data, code or files here may be modified.

A notification will be released when there is a decision on the paper.

Notes about data and code
1)code in in the folder of FDUM is our central case, and the results are displayed in fig 2. And code in in the folder of FDUM incorporates sensitivity tests, and the results are shown in fig 3.

2） the main code for both two parts are run in fumccm.m, which produces files in the folder of '../fdum/'. This is the main code of the model, which calls ModelCN.m to run simulations from 1971 to 2300 and calculate the utility based on the welfare dynamics.

3)  after running  fumccm.m, you can run A1plot_dashboard_fig2.m to produce fig 2 in our paper.

4)  or after running  fumccm.m, you can run A1plot_nscc_fig3_readdata.m which execute the post-simulation processing based on the files generated in '..fdum'. The output of  post-simulation processing are stored in the folder of '..nuclear'.

5) after producing the files in '..nuclear' by running  A1plot_nscc_fig3_readdata.m , you can run A1plot_nscc_fig3_plotfig.m to produce fig 3 in our paper.

Rong Wang and Ruipu Yang

13 October 2023

Contact rongwang@fudan.edu.cn