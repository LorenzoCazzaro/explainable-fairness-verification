# Explainable Global Fairness Verification

This repository contains the implementation of the synthesiser of sufficient conditions for fairness for decision-tree ensembles proposed by **Calzavara *et. al.*** in their research paper titled [<em>Explainable Global Fairness Verification of Tree-Based Classifiers</em>](https://openreview.net/forum?id=HOu7LgqCTqd). At the moment, this repository contains the code of the synthesiser. The scripts to scripts to reproduce the experiments described in the paper will be added soon.

## Installation
Download the repository. Remember to compile using the flags <em>-Iinclude</em> and <em>-lpthread</em> to use the synthesiser.

### Requirments
The following requirements are necessary to use the analyzer. For reproducibility reasons, specific versions of the libraries are required.
- Python3
- C++14 compiler
- sklearn (version 0.23.2 used for the experiments)
- pandas (version 1.1.3 used for the experiments)
- numpy (version 1.19.1 used for the experiments)
- boost C++ libraries
  
### Data-Independent Stability Analysis
The resilience-verification repository at [https://github.com/FedericoMarcuzzi/resilience-verification](https://github.com/FedericoMarcuzzi/resilience-verification) contains the instructions to execute the data-independent stability analysis on the tree-based classifier.
At the moment, this repo contains the essential code to perform the data-independent stability analysis. 

## Usage
Compile the synthesiser using the following string:
```bash
g++ g++ -o synthesiser ./src/exec_analyser_synthesiser.cpp -Iinclude -lpthread
```

## Run the experiments presented in the paper
Coming soon...
