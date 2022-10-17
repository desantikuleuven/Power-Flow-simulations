# Power Flow simulations

This repo contains different scripts for computing Power Flows by exploiting PowerModels.jl

File named My_ref contains all the necessary PowerModel references for building the PF models. 

Most relevant files:

1) Basic_pf.jl --> performs AC PF as defined in PowerModels
2) Flex_pf.jl and Lazy_flex_pf.jl --> the latter script adds voltage and congestions constraints iteratively only when violated. Done just to improve speed.
It considers upwards and downwards flexibility offered by all nodes. DGs are fixed. 
3) DG_curtail_and_flex_lazy,jl --> allows to perform AC PF by including as well DG curtailment. 
4) My_functions.jl --> contains all the relevant functions used in all scripts 
