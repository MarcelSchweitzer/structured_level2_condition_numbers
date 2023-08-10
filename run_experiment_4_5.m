%% Run experiment for matrix exponential of triangular matrices generated as Schur factors of benchmark matrices (Experiment 4.5 in manuscript)
% Warning: This code is not optimized for speed, computations take quite
% some time

% As the matrices in this example come from various different sources,
% toolboxes etc., we collected them in the file
% data/triangular_benchmark_matrices.mat to not overcrowd this repository
% with too many things from other sources. Details on the origin of the
% all test matrices are given in Appendix A of [B. Arslan, S. D. Relton, M.
% Schweitzer, Structured level-2 condition numbers of matrix functions]

addpath('aux_functions')
addpath('main_algorithms')

rng(1);
[x,l_2,sl_2,slb_2,lb_2, CN] = struc_level2_benchmark();

[CN_sorted, sort_idx] = sort(CN);
l_2 = l_2(sort_idx);
sl_2 = sl_2(sort_idx);
slb_2 = slb_2(sort_idx);
lb_2 = lb_2(sort_idx);

loglog(CN_sorted,l_2,'b-*','markersize',8)
hold all
loglog(CN_sorted,lb_2,'k:x','markersize',8)
loglog(CN_sorted,sl_2,'r--d','markersize',8)
loglog(CN_sorted,slb_2,'g-o','markersize',8)
hold off
legend('ub\_uscond2','lb\_uscond2','ub\_scond2','lb\_scond2','Location','NorthWest')
xlabel('\kappa_2(U)')
xlim([1e0,1e22])
exportgraphics(gcf,'figures/expm_triu_benchmark.pdf')