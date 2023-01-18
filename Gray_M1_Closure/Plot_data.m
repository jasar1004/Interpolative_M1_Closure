clear all
close all
clc

%%
Lagrange = importdata('Lag_M1_1D_000000.dat');
fclose all;

i = 50;
npts_E = 100;
npts_f = 99;
Index_l = i*npts_f + 1;
Index_u = (i+1)*npts_f;
Lagrange = Lagrange(Index_l:Index_u,:);

figure
plot(Lagrange(:,2), Lagrange(:,3))
xlabel('N^{1}')
ylabel('\lambda_{0}')

figure
plot(Lagrange(:,2), Lagrange(:,4))
xlabel('N^{1}')
ylabel('\lambda_{1}')

%%
Partial_Mom_file_Id = fopen('Partial_Moms_M1_1D.dat');
Partial_Mom = fread(Partial_Mom_file_Id, [2 20000],'double');
fclose all;
% n_Part_Mom = numel(Partial_Mom);
Partial_Mom_E_plus = Partial_Mom(1,:)';
Partial_Mom_F_plus = Partial_Mom(2,:)';

Partial_Mom_E_plus = reshape(Partial_Mom_E_plus, [200, 100]);
Partial_Mom_F_plus = reshape(Partial_Mom_F_plus, [200, 100]);

Partial_Mom_E_plus = sort(Partial_Mom_E_plus, 1);
Partial_Mom_F_plus = sort(Partial_Mom_F_plus, 1);

% Partial_Mom_E_plus

figure
plot(Partial_Mom_E_plus)
ylabel('E^{+}')

figure
plot(Partial_Mom_F_plus)
ylabel('F^{+}')