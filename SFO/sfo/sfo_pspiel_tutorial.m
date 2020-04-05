%% Matlab Toolbox for Submodular Function Optimization (v 2.0)
% Tutorial script and all implementations Andreas Krause (krausea@gmail.com)
%
% Tested in MATLAB 7.0.1 (R14), 7.2.0 (R2006a), 7.4.0 (R2007a, MAC), 7.9.0 (MAC)
% 
% First some information on conventions:
% --------------------------------------
% All algorithms will use function objects.
% For example, to measure variance reduction in a Gaussian model, call 
%   F = sfo_fn_varred(sigma,V)
% where sigma is the covariance matrix and V is the ground set, e.g., 1:size(sigma,1)
%
% A note on Octave compatibility:
% ------------------------------
% This toolbox also works under Octave; however, since Octave handles
% function objects differently from Matlab. Use the function sfo_octavize
% to make a submodular function object Octave ready; 
% type 'help sfo_octavize' for more information.
% The script sfo_tutorial_octave has been tested under Octave 3.2.3
%
% Implemented algorithms:
% ----------------------- 
%
% Minimization:
%
% sfo_min_norm_point: Fujishige's minimum norm point algorithm for minimizing
%       general submodular functions 
% sfo_queyranne: Queyranne's algorithm for minimizing symmetric submodular
%       functions
% sfo_ssp: Submodular-supermodular procedure of Narasimhan & Bilmes for 
%       minimizing the difference of two submodular functions
% sfo_s_t_min_cut: For solving min F(A) s.t. s in A, t not in A
% sfo_minbound: Return an online bound on the minimum solution
% sfo_greedy_splitting: Greedy splitting algorithm for clustering of 
%       Zhao et al
%
% Maximization:
%
% sfo_polyhedrongreedy: For solving an LP over the submodular polytope
% sfo_greedy_lazy: The greedy algorithm for constrained maximization / coverage
% sfo_greedy_welfare: The greedy algorithm for solving allocation problems
% sfo_cover: Greedy coverage algorithm
% sfo_celf: The CELF algorithm for budgeted maximization
% sfo_ls_lazy: Local search algorithm for maximizing nonnegative submodular functions
% sfo_saturate: The Saturate algorithm of Krause et al. for robust optimization of submodular
%       functions
% sfo_max_dca_lazy: The Data Correcting algorithm of Goldengorin et al. for 
%       maximizing general (not necessarily nondecreasing) submodular functions
% sfo_maxbound: Return an online bound on the maximum solution
% sfo_pspiel: pSPIEL algorithm for trading off information and
%       communication cost
% sfo_pspiel_orienteering: pSPIEL algorithm for submodular orienteering
% sfo_balance: eSPASS algorithm for simultaneous placement and balanced scheduling
%
% Miscellaneous
%
% sfo_lovaszext: Computes the Lovasz extension for a submodular function
% sfo_mi_cluster: Clustering algorithm using both maximization and
%       minimization
% sfo_pspiel_get_path: Convert a tree into a path using the MST heuristic
%       algorithm
% sfo_pspiel_get_cost: Compute the Steiner cost of a tree / path
%
% Submodular functions (also try 'help sfo_fn')
% 
% sfo_fn_cutfun: Cut function
% sfo_fn_detect: Outbreak detection / facility location
% sfo_fn_entropy: Entropy of Gaussian random variables
% sfo_fn_infogain: Information gain about gaussian random variables
% sfo_fn_mi: Gaussian mutual information I(A; V\A)
% sfo_fn_varred: Gaussian variance reduction (orthogonal matching pursuit)
% sfo_fn_example: 2 element example from tutorial
% sfo_fn_iwata: Iwata's test function for testing minimization
% sfo_fn_ising: Energy function for Ising model for image denoising
% sfo_fn_residual: For defining residual submodular functions
% sfo_fn_invert: For defining F(A) = F'(V\A)-'F(V)
% sfo_fn_lincomb: For defining linear combinations of submodular functions
%
% Example: sfo_tutorial
%
% If you use the toolbox for your research, please cite
% A. Krause. "SFO: A Toolbox for Submodular Function Optimization". Journal
%   of Machine Learning Research (2010). 
%
% Here is an overview reference for submodularity in AI
% A. Krause, C. Guestrin. "Near-optimal Observation Selection Using Submodular Functions". 
%   Survey paper, Proc. of 22nd Conference on Artificial Intelligence (AAAI) 2007 -- Nectar Track
%
% Change log
% -----------
%
% Version 2.0
% * Modified specification of optional parameters (using sfo_opt)
% * Added sfo_ls_lazy for maximizing nonnegative submodular functions
% * Added sfo_fn_infogain, sfo_fn_lincomb, sfo_fn_invert, ...
% * Added additional documentation and more examples
% * Now Octave ready
%
% Version 1.1
% * added pSPIEL for informative path planning
% * added eSPASS for simultaneous placement and scheduling
% * new convention for submodular functions (incremental computations,
%   etc.) Much faster!
%

do_plot = 1; %switch visualization on or off

disp(' ');
disp('---------------------------------------------------------');
disp('---------------------------------------------------------');
disp('Welcome to the submodular function optimization tutorial');
disp('---------------------------------------------------------');
disp('---------------------------------------------------------');
disp(' ')
if ~do_plot
    disp('Visualization disabled');
end
pause

%% INITIALIZATION

% We first define several examples of submodular functions
% In order to obtain more information about submodular function objects,
% type 'help sfo_fn'
%
% First some functions for experimental design on a Gaussian Process
% trained on pH data from a lake in Merced, California

% load the data: Contains covariance matrix merced_data.sigma, 
% and locations (coordinates) merced_data.coords
load merced_data;
% the ground set for the submodular functions
V_sigma = 1:size(merced_data.sigma,1); 

% Mutual information: F_mi(A) = H(V\A) - H(V\A | A)
F_mi = sfo_fn_mi(merced_data.sigma,V_sigma);
% Mutual information: F_mi(A) = H(V\A) - H(V\A | A)
F_ig = sfo_fn_infogain(merced_data.sigma,V_sigma,.01);
% Variance reduction: F_var(A) = Var(V)-Var(V | A)
F_var = sfo_fn_varred(merced_data.sigma,V_sigma);

% Helper function for evaluating the maximum remaining variance,
% eval_maxvar(A) = max_s Var(s | A);
eval_maxvar = @(A) sfo_eval_maxvar(merced_data.sigma,A); 

% *Note:* mutual information is submodular, but not monotonic! But for when
% selecting sets of size k << n (where n is the total # of elements) it's
% approximately monotonic, and that's enough (see Krause et al., JMLR '08)
% Also, variance reduction is not always submodular (see Das & Kempe, STOC '08)

% -----------------------------------------------------
% Now the 2 element example function from the tutorial slides at www.submodularity.org
% 'F{[]} = 0, F({a}) = -1, F({b}) = 2, F({a,b}) = 0';
F_ex = sfo_fn_example;
V_ex = 1:2;

% -----------------------------------------------------
% now some cut functions on directed and undirected graphs
% first define adjacency matrices
G_dir=[0 1 1.2 0 0 0; 1.3 0 1.4 0 0 0; 1.5 1.6 0 1.7 0 0; 0 0 1.8 0 1.9 2; 0 0 0 2.1 0 2.2; 0 0 0 2.3 2.4 0];
G_un=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
% here's the cut functions:
F_cut_dir = sfo_fn_cutfun(G_dir);
F_cut_un = sfo_fn_cutfun(G_un);
% and the ground set V_G
V_G = 1:6;

% -----------------------------------------------------
% now an objective function for modeling sensor detections
n_sensors = 100;
n_targets = 50;
detprob = rand(n_sensors,n_targets);
V_det = 1:n_sensors;
F_det = sfo_fn_detect(detprob,V_det);

% -----------------------------------------------------
% a test function by Iwata for testing submodular minimization algorithms
% described, e.g., in Fujishige '06
V_iw = 1:100;
F_iw = sfo_fn_iwata(length(V_iw));



%% The pSPIEL Algorithm for trading off accuracy and communication cost

disp(' ')
disp('---------------------------------------------------------');
disp('Now let''s use the pSPIEL algorithm to trade off informativeness and')
disp('communication cost in GP regression on Merced Lake')
disp(' ');
V = V_sigma;
D = merced_data.dists + ones(length(V)); %cost of links + cost of nodes
Q = 0.6*F_var(V); %want 80% of optimal variance reduction

disp('Run greedy algorithm: ')
disp('AG = sfo_cover(F_var,V,Q)');
AG = sfo_cover(F_var,V,Q)

disp('AP = sfo_pspiel(F_var,V,Q,D)');
AP = sfo_pspiel(F_var,V,Q,D)

utility_greedy = F_var(AG);
utility_pspiel = F_var(AP);
[cost_greedy edges_greedy steiner_greedy] = sfo_pspiel_get_cost(AG,D);
[cost_pspiel edges_pspiel steiner_pspiel]= sfo_pspiel_get_cost(AP,D);

if do_plot
    subplot(211)
    sfo_plot_subgraph(merced_data.coords,edges_greedy,steiner_greedy)
    title('Greedy-connect');
    subplot(212)
    sfo_plot_subgraph(merced_data.coords,edges_pspiel,steiner_pspiel)
    title('pSPIEL');
end

disp(sprintf('Greedy-connect: Utility = %f, Cost = %f. pSPIEL: Utility = %f, Cost = %f.',utility_greedy,cost_greedy,utility_pspiel,cost_pspiel));
pause

