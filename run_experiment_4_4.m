%% Run experiment for matrix exponential of random triangular matrices (Experiment 4.4 in manuscript)
% Warning: This code is not optimized for speed, computations take quite
% some time

addpath('aux_functions')
addpath('main_algorithms')

rng(1);
n = 8;
[x,r,l_2,sl_2,slb_2,lb_2, CN] = struc_level2(n,'triangular',@expm,'expm','T',[]);

loglog(CN,l_2,'b-*','markersize',8)
hold all
loglog(CN,lb_2,'k:x','markersize',8)
loglog(CN,sl_2,'r--d','markersize',8)
loglog(CN,slb_2,'g-o','markersize',8)
hold off
legend('ub\_uscond2','lb\_uscond2','ub\_scond2','lb\_scond2','Location','NorthWest')
xlabel('\kappa_2(U)')
exportgraphics(gcf,'figures/expm_triu.pdf')