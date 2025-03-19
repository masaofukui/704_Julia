# EC 704 Codes (Boston University)

This repository has codes for the second part of EC704, a part of the first-year macro sequence at Boston University, as taught in spring 2025.

The course covers labor market and financial market frictions.
* Topic 1: Unemployment facts.
  * [lecture1.do](./Topic1/lecture1.do): Top-level code to produce all the figures in lecture note 1. You need to obtain the FRED API key and "set fredkey" in Stata.
  * [time_aggreation.R](./Topic1/time_aggregation.R): This needs to be run before line 383 of the above file.
  * The underlying data is too large to be uploaded on github. You can download the data from [here](https://www.dropbox.com/scl/fo/aofw98nppaey0pjoi1aql/AL0pFltslJinYj_tKfWq4IY?rlkey=a9v0i1ugkf30ik5zyiiam77eb&dl=1). Save the folder as "./Topic1/oriignal_data"

* Topic 6: Financial Friction and Capital Misallocation
  * [Toplevel.jl](./Topic6/Toplevel.jl): Top-level code to solve the steady state and the transition dynamics of the discrete-time version of [Moll (2014)](https://benjaminmoll.com/wp-content/uploads/2019/07/TFPFF.pdf).

* Topic 7: Canonical Incomplete Market Models (co-written with Aruzhan Nurlankul)
  * [Toplevel.jl](./Topic7/Toplevel.jl): Top-level code to solve the steady state and the transition dynamics of Bewley-Hugget-Aiyagari model. Transition dynamics is obtained by sequence space Jacobian method by [Auclret, Bardóczy, Rognlie and Straub (2021)](https://web.stanford.edu/~aauclert/sequence_space_jacobian.pdf).


