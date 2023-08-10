%% Run experiment for matrix logarithm (Experiment 4.1 in manuscript)
% Warning: This code is not optimized for speed, computations take quite
% some time

addpath('toolboxes/jagger')
addpath('aux_functions')
addpath('main_algorithms')
rng(1)

disp('--- Running experiment for logm with orthogonal test matrices ---')
[x_orth,r_orth,l_2_orth,sl_2_orth,slb_2_orth,lb_2_orth, CN_orth] = struc_level2(4,'orth',@logm,'logm','I','auto');

disp('--- Running experiment for logm with symplectic test matrices ---')
[x_symp,r_symp,l_2_symp,sl_2_symp,slb_2_symp,lb_2_symp, CN_symp] = struc_level2(4,'symp',@logm,'logm','J','auto');

disp('--- Running experiment for logm with perplectic test matrices ---')
[x_perp,r_perp,l_2_perp,sl_2_perp,slb_2_perp,lb_2_perp, CN_perp] = struc_level2(4,'perp',@logm,'logm','R','auto');


% Plot results for orthogonal matrices
figure(1)
[l_2_orth,sort_idx]= sort(l_2_orth);
sl_2_orth = sl_2_orth(sort_idx);
lb_2_orth = lb_2_orth(sort_idx);
slb_2_orth = slb_2_orth(sort_idx);
r_orth = r_orth(sort_idx);

semilogy(x_orth,l_2_orth,'b-*','markersize',8)
hold on
semilogy(x_orth,lb_2_orth,'k:x','markersize',8)
semilogy(x_orth,sl_2_orth,'r--d','markersize',8)
semilogy(x_orth,slb_2_orth,'g-.o','markersize',8)
hold off
legend('ub\_uscond2','lb\_uscond2','ub\_scond2', 'lb\_scond2','Location','north west')
xlabel('Number of test matrices')
exportgraphics(gcf,'figures/logm_orth.pdf')

% Plot results for symplectic matrices 
figure(2)
loglog(CN_symp,l_2_symp,'b-*','markersize',8)
hold on
loglog(CN_symp,lb_2_symp,'k:x','markersize',8)
loglog(CN_symp,sl_2_symp,'r--d','markersize',8)
loglog(CN_symp,slb_2_symp,'g-.o','markersize',8)
legend('ub\_uscond2','lb\_uscond2','ub\_scond2', 'lb\_scond2','Location','north west')
xlabel('\kappa_2(A)')
hold off

set(gca,'YminorTick','off')
set(gca,'XminorTick','off')
exportgraphics(gcf,'figures/logm_symp.pdf')


% Plot results for perplectic matrices 
figure(3)
loglog(CN_perp,l_2_perp,'b-*','markersize',8)
hold on
loglog(CN_perp,lb_2_perp,'k:x','markersize',8)
loglog(CN_perp,sl_2_perp,'r--d','markersize',8)
loglog(CN_perp,slb_2_perp,'g-.o','markersize',8)
legend('ub\_uscond2','lb\_uscond2','ub\_scond2', 'lb\_scond2','Location','north west')
xlabel('\kappa_2(A)')
hold off

set(gca,'YminorTick','off')
set(gca,'XminorTick','off')
exportgraphics(gcf,'figures/logm_perp.pdf')