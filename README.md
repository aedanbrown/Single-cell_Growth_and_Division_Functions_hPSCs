Single-cell Growth and Division Functions in hPSCs
==================================================
This repository contains the code used to analyze experimental and synthetic data to determine physiological state functions (PSFs) of human pluripotent stem cells.

# Running the code
Each file needs appropriate data to be run.

| File name | Comment |
|-|-|
|inverseProblem_synthetic_diff.m| Compute the PSFs of differentiating cells with a single instance of data.|
|inverseProblem_synthetic_diff_multi.m| Compute the PSFs of differentiating cells with multiple instances of data.|
|inverseProblem_synthetic_growth.m| Compute the PSFs of exclusively growing cells with a single instance of data.|
|inverseProblem_synthetic_growth_multi.m| Compute the PSFs of exclusively growing cells with multiple instances of data.|
|inverseProblem_experimental_H9.m| Analyze the experimental H9 data and compute the PSFs. The code is the same as inverseProblem_experimental_IMR90.m, with the initial inputs adjusted. |
|inverseProblem_experimental_IMR90.m| Analyze the experimental H9 data and compute the PSFs. The code is the same as inverseProblem_experimental_H9.m, with the initial inputs adjusted. |


