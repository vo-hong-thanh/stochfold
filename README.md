Stochfold implements different kinetic algorithms for predicting the folding of RNA. It includes
  1) StochFold algorithm: predicting folding kinetics of RNA. The folding is defined by element actions on a base pair, i.e., adding a base pair, deleting a base-pair and shifting a base pair, following the Watson-Crick rule. The transition between RNA conformations defines a continuous-time Markov process where its rate is measured by the Metropolis rule. The algorithm implements the kinetic Monte Carlo (also known as the Gillespie algorithm) to realize a folding path.  
  2) CostochFold algorithm: predicting cotranscriptional folding of RNA by taking into account its transcription process. RNA folds simultaneously after new nucleotides are synthesized during the transcription. The transcription process is explicitly considered by including the action of elongation. In elongation, the current RNA chain increases in length and a newly synthesized nucleotide is added to its 3' end.

The extension of these folding algorithms to include pseudoknots was developed by Dani Korpela.  
## Main reference:
Vo Hong Thanh, Dani Korpela and Pekka Orponen. Cotranscriptional kinetic folding of RNA secondary structures including pseudoknots. _Journal of Computational Biology_ 28(9):892-908. doi: 10.1089/cmb.2020.0606, 2021.

## Further readings and references:
The theoretical background in RNA folding is in [1]. The Zuker algorithm is developed by Zuker and Stiegler [2]. The stochastic folding of RNA is developed by Flamm *et al.* [3]. The book by Marchetti et al. [5] gives a comprehensive review and recent developments in stochastic simulation. 

[1] Jörg Fallmann, Sebastian Will, Jan Engelhardt, Björn Grüning, Rolf Backofenc, and Peter F. Stadler. Recent advances in RNA folding. J. Biotechnol., 261:97–104, 2017.

[2] Michael Zuker and Patrick Stiegler. Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information,	Nucleic Acids Res. 9(1):133-148, 1981.

[3] Christoph Flamm, Walter Fontana, Ivo L. Hofacker, and Peter Schuster. RNA folding at elementary step resolution. RNA, 6:325–338, 2000.

[4] Marchetti, L., Priami, C. and Thanh, V.H. Simulation Algorithms for Computational Systems Biology, Springer. 2017. https://www.springer.com/gp/book/9783319631110.
