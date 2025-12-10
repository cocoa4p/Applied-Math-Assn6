
function modal_continuous_vs_discrete()

    total_mass = 5;
    tension_force = 5;
    string_length = 10;

    % continuous mode

    c = sqrt(tension_force / (total_mass / string_length)); % wave speed

    mode_index = 3;  

    % continuous x-grid
    x_cont = linspace(0, string_length, 1000);

    % continuous mode shape
    X_cont = sin(mode_index * pi * x_cont / string_length);

    omega_cont = mode_index * pi * c / string_length;
    fprintf("Continuous natural frequency = %.5f\n", mode_index, omega_cont);

    % discrete modes for varying N
 
    N_list = [4, 8, 20, 80]; % different numbers of masses
    
    figure; 
    for i = 1:length(N_list)
        N = N_list(i);

        dx = string_length/(N+1);

        % build param struct
        string_params = struct();
        string_params.n = N;
        string_params.M = total_mass;
        string_params.Tf = tension_force;
        string_params.L = string_length;
        string_params.c = 0;
        string_params.dx = dx;

        % get matrices
        [M_mat, K_mat] = construct_2nd_order_matrices(string_params);

        % solve eigenproblem
        [V, D] = eig(K_mat, M_mat);

        % discrete nat frequencies
        omega_disc_all = sqrt(diag(-D));   

        v_disc = V(:, mode_index);
        
        % discrete x-grid
        x_disc = linspace(0, string_length, N+2);
        v_disc_full = [0; v_disc; 0];

        % matches continuous amplitude vector size
        scale = max(abs(X_cont)) / max(abs(v_disc_full));
        v_disc_full = scale * v_disc_full;

        % plot comparison
%         subplot(2,2,i);
%         plot(x_cont, X_cont, 'b-', 'LineWidth', 2); hold on;
%         plot(x_disc, v_disc_full, 'r.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red');
%         title(sprintf('N = %d Masses', N));
%         xlabel('x'); ylabel('Mode shape');
%         legend('Continuous', 'Siscrete');

        % Semilog Error Plot for Frequencies
        % frequency approximation error
        N_vals = 4:4:200;
        freq_error = zeros(size(N_vals));
    
        for j = 1:length(N_vals)
            N = N_vals(j);
            dx = string_length/(N+1);
    
            string_params = struct();
            string_params.n = N;
            string_params.M = total_mass;
            string_params.Tf = tension_force;
            string_params.L = string_length;
            string_params.c = 0;
            string_params.dx = dx;
    
            [M_mat,K_mat] = construct_2nd_order_matrices(string_params);
            [~,D] = eig(K_mat,M_mat);
            omega_disc = sqrt(-D(mode_index,mode_index));
    
            freq_error(j) = abs(omega_disc - omega_cont);
        end
    
        figure;
        semilogy(N_vals, freq_error, 'LineWidth', 2);
        xlabel('Number of Masses N');
        ylabel('|Discrete Freq - Continuous Freq|');
        title(sprintf('Frequency Approximation Error (Mode %d)', mode_index));

    end

    % harmonic comparison discrete vs continuous
    N_low = 4;
    N_high = 80;

    figure;
    tiledlayout(1,2);

    % ----- low N -----
    nexttile;
    plot_harmonics(N_low, total_mass, tension_force, string_length);
    title('Harmonics: N = 8');

    % ----- high N -----
    nexttile;
    plot_harmonics(N_high, total_mass, tension_force, string_length);
    title('Harmonics: N = 80');

end


function plot_harmonics(N, M_total, Tf, L)
    dx = L/(N+1);
    c = sqrt(Tf/(M_total/L));

    string_params = struct();
    string_params.n = N;
    string_params.M = M_total;
    string_params.Tf = Tf;
    string_params.L = L;
    string_params.c = 0;
    string_params.dx = dx;

    [M_mat, K_mat] = construct_2nd_order_matrices(string_params);
    [~, D] = eig(K_mat, M_mat);

    omega_disc = sqrt(abs(diag(D)));

    % continuous frequencies
    n_vals = 1:N;
    omega_cont = n_vals * pi * c / L;

    plot(n_vals, omega_disc(1:N), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
    plot(n_vals, omega_cont, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 5);
    xlabel('Mode Index n');
    ylabel('Frequency ω');
    legend('Discrete','Continuous');
end

