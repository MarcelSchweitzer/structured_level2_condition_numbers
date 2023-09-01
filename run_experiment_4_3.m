%% Run experiment for matrix exponential (Experiment 4.3 in manuscript)
% Warning: This code is not optimized for speed, computations take quite
% some time

addpath('toolboxes/jagger')
addpath('aux_functions')
addpath('main_algorithms')
rng(1)

disp('--- Running experiment for expm with skew-symmetric test matrices ---')
[x_skew,r_skew,l_2_skew,sl_2_skew,slb_2_skew,lb_2_skew, CN_skew] = struc_level2(4,'skew-sym',@expm,'expm','I','lie');

disp('--- Running experiment for expm with Hamiltonian test matrices ---')
[x_hamil,r_hamil,l_2_hamil,sl_2_hamil,slb_2_hamil,lb_2_hamil, CN_hamil] = struc_level2(4,'Hamilton',@expm,'expm','J','lie');

% Plot results for orthogonal matrices
figure(1)
[l_2_skew,sort_idx]= sort(l_2_skew);
sl_2_skew = sl_2_skew(sort_idx);
lb_2_skew = lb_2_skew(sort_idx);
slb_2_skew = slb_2_skew(sort_idx);
r_skew = r_skew(sort_idx);

semilogy(x_skew,l_2_skew,'b-*','markersize',8)
hold on
semilogy(x_skew,sl_2_skew,'r--d','markersize',8)
semilogy(x_skew,lb_2_skew,'k:x','markersize',8)
semilogy(x_skew,slb_2_skew,'g-.o','markersize',8)
hold off
legend('ub\_uscond2','ub\_scond2','lb\_uscond2', 'lb\_scond2','Location','west')
ylim([1e-14, 1e1])
xlabel('Number of test matrices')
exportgraphics(gcf,'figures/expm_skew.pdf')

% Plot results for symplectic matrices 
figure(2)
[l_2_hamil,sort_idx]= sort(l_2_hamil);
sl_2_hamil = sl_2_hamil(sort_idx);
lb_2_hamil = lb_2_hamil(sort_idx);
slb_2_hamil = slb_2_hamil(sort_idx);
r_hamil = r_hamil(sort_idx);

semilogy(x_hamil,l_2_hamil,'b-*','markersize',8)
hold on
semilogy(x_hamil,lb_2_hamil,'k:x','markersize',8)
semilogy(x_hamil,sl_2_hamil,'r--d','markersize',8)
semilogy(x_hamil,slb_2_hamil,'g-.o','markersize',8)
xlabel('Number of test matrices')
legend('ub\_uscond2','lb\_uscond2','ub\_scond2','lb\_scond2','Location','north west')
ylim([9e-1, 1e1])
hold off

set(gca,'YminorTick','off')
set(gca,'XminorTick','off')
exportgraphics(gcf,'figures/expm_hamil.pdf')
