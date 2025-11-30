string_params = struct();
string_params.n = 3;
string_params.M = 5;
string_params.Tf = 5;
string_params.L = 10;
string_params.c = 0.001;
string_params.dx = string_params.L/(string_params.n+1);

[M_mat, K_mat] = construct_2nd_order_matrices(string_params);

% Generalized eigenvalue
[Phi, Lambda] = eig(K_mat, M_mat);
% frequencies
omega_vals = sqrt(diag(Lambda));

[omega_vals, sort_idx] = sort(omega_vals);

Phi = Phi(:, sort_idx);

% normalize
for k = 1:length(omega_vals)
    Phi(:,k) = Phi(:,k) / max(abs(Phi(:,k)));
end

disp('Natural frequencies (rad/s):');
disp(omega_vals);

disp('Mode shapes (each column is a mode):');
disp(Phi);


xlist = linspace(0, string_params.L, string_params.n+2);
figure; 
hold on;

for k = 1:length(omega_vals)
    y = [0; Phi(:,k); 0];
    plot(xlist, y, 'LineWidth', 2);
end

title('Mode Shapes of Discrete String');
xlabel('X'); 
ylabel('Displacement');
legend("Mode 1","Mode 2","Mode 3");
