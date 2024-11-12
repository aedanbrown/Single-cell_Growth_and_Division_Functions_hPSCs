Determination of single-cell growth rate and division intensity of human pluripotent stem cells
==================================================
This repository contains the code used to analyze experimental and synthetic data to determine physiological state functions (PSFs) of human pluripotent stem cells.

# Running the code
The code is split into files that analyze the synthetic results and those that analyze the experimental results. With the synthetic results, the code to generate synthetic data is included as well. The experimental data is not included in this repository, but the code related to experiments is provided as an example of the analysis.

## Synthetic Results

| File name | Comment |
|-|-|
|Synthetic/inverseProblem_synthetic_diff.m| Compute the PSFs of differentiating cells with a single instance of data.|
|Synthetic/inverseProblem_synthetic_diff_multi.m| Compute the PSFs of differentiating cells with multiple instances of data.|
|Synthetic/inverseProblem_synthetic_growth.m| Compute the PSFs of exclusively growing cells with a single instance of data.|
|Synthetic/inverseProblem_synthetic_growth_multi.m| Compute the PSFs of exclusively growing cells with multiple instances of data.|
|Synthetic/diffSimIQ.m| Simulate the PBE with specified PSFs, including differentiation. Generates files that can be analyzed with the inverseProblem_synthetic_diff files.|
|Synthetic/growthSimIQ.m |Simulate the PBE with specified PSFs, without differentiation (only growth). Generates files that can be analyzed with the inverseProblem_synthetic_growth files.|

## Experimental Results

The main file of interest is inverseProblem_experimental_main.m, which would generate the main figures of the paper if given the experimental data.

| File name | Comment |
|-|-|
|Experimental/inverseProblem_experimental_main.m|Uses the other functions in the Experimental folder to estimate the PSFs from each experimental condition and plot the results.|
|Experimental/inverseProblem_experimental_H9_func.m| Analyze the experimental H9 data and compute the PSFs. |
|Experimental/inverseProblem_experimental_H9_lactate_func.m| Analyze the experimental H9+lactate data and compute the PSFs. |
|Experimental/inverseProblem_experimental_IMR90_func.m| Analyze the experimental IMR90 data and compute the PSFs. |
|Experimental/inverseProblem_experimental_IMR90_lactate_func.m| Analyze the experimental IMR90+lactate data and compute the PSFs. |
|Experimental/estimatePSFs.m|Helper function to calculate PSFs given the various distributions.|


