# Epidemic tracing through Statistical Inference


## Table of contents
* [Message-Passing strategies for epidemic tracing](#message-passing-strategies-for-epidemic-tracing)
* [Epidemic Control](#epidemic-control)
* [Inference of epidemic parameters](#inference-of-epidemic-parameters)

## Message-Passing strategies for epidemic tracing

Description:

- Tracing
- BP
- MF

## ROCS
ROC curves and precisions at differents epidemic sizes. 

![](./figs/roc_.png)

Averaged ROC area at different epidemic size, changing app adoptions (100\%,66\%, 62\%, 55\%) :

![](./figs/auc.gif)

## Epidemic control

OpenABM model API: brief description and link to for

Intervention results: 10K - single

Intervention results: 50K - multiple

![intervention_multiple_50K](https://github.com/sibyl-team/sib/blob/master/examples/figs/N50K_o400_linear_and_log.svg)


## Inference of epidemic parameters

Temporal graph with 10K nodes and contacts up to day T= 30.
After 5initial days, we select 10 % of the nodes are observed uniformly at random.
Their states are observed on a daily basis.

![inference_auc_parameters_10K](https://github.com/sibyl-team/sib/blob/master/examples/figs/inference_parameters_openABM_gamma.png)
