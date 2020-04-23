# AQI-deploy
Optimizing sensor deployment and maintenance costs in large-scale environmental monitoring.

## Getting Started

Test environment: MATLAB R2019b/R2020a.

**Note:** Need the Curve Fitting Toolbox for the `fit` function.  Need the Bioinformatics Toolbox for the `graphminspantree` function.

The `tutorial.m` will walk you through all of our algorithms on the small dataset.

## File Structure

```
.
├── LICENSE
├── README.md     // this file
├── SFO           // the SFO toolbox by Krause et al. from FileExchange
├── alg           // our algorithms
├── data-large    // large dataset around LA after preprocess
├── data-small    // small dataset around SD after preprocess
├── exp           // scripts to run experiments
├── gp            // sensing quality library based on Gaussian Process
├── libs          // general library
├── lldistkm      // distance calculation library from FileExchange
├── mlibs         // maintenance cost library
└── tutorial.m
```

## Available Algorithms

* Distance-Weighted Greedy (DWG) [[Awerbuch et al. 1999]](https://www.cs.tau.ac.il/~azar/kmst.pdf).
* Information-Driven Sensor Querying (IDSQ) [[Zhao and Guibas 2004]](https://books.google.com/books?hl=en&lr=&id=BkaQkhkWGfoC&oi=fnd&pg=PP2&dq=wireless+sensor+networks:+an+information+processing+approach&ots=HVaEhKt6Q9&sig=Cst8PgUBGemOE0nt4yd6ujTyTro#v=onepage&q=wireless%20sensor%20networks%3A%20an%20information%20processing%20approach&f=false).
* Padded Sensor Placements at Informative and cost-Effective Locations (pSPIEL) [[Krause et al. 2011]](https://dl.acm.org/doi/10.1145/1921621.1921625). We download their implementation of Submodular Function Optimization toolbox from FileExchange [here](https://www.mathworks.com/matlabcentral/fileexchange/20504-submodular-function-optimization).
* Particle Swarm Optimization (PSO) [[Clerc and Jennedy 2002]](https://ieeexplore.ieee.org/abstract/document/985692).
* Artificial Bee Colony (ABC) Optimization [[Karaboga and Basturk 2007]](https://link.springer.com/article/10.1007/s10898-007-9149-x).

For PSO and ABC, we devise the cost function for our sensor deployment problem in `./alg/CostFunction.m`.

## License

MIT

