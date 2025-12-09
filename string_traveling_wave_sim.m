function string_traveling_wave_sim()

    num_masses = 300;
    total_mass = 5;
    tension_force = 5;
    string_length = 10;
    damping_coeff = 0.0001;  

    dx = string_length/(num_masses+1);

    % wave speed
    c = sqrt(tension_force / (total_mass/string_length));
    fprintf("Wave speed c = %.3f\n", c);

    xlist = linspace(0, string_length, num_masses+2);  

    % Triangle pulse parameters
    pulse_width = string_length/4;
    pulse_height = 1.0;

    % Right-end pulse function and derivative using your triangle_pulse functions
    Uf_func = @(t) triangle_pulse(t, pulse_width/c, pulse_height);
    dUfdt_func = @(t) triangle_pulse_derivative(t, pulse_width/c, pulse_height);

    % Store string parameters
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;

    % ODE function
    my_rate_func = @(t, V) string_rate_func01(t, V, string_params);

    % Initial conditions
    U0 = zeros(num_masses,1);
    dU0 = zeros(num_masses,1);
    V0 = [U0; dU0];

    T_end = 4 * string_length / c;  % several passes across the string
    tspan = [0, T_end];

    [tlist, Vlist] = ode45(my_rate_func, tspan, V0);

    % extract displacement only
    Ulist = Vlist(:, 1:num_masses); 

    figure; hold on;

    hString = plot(xlist, zeros(size(xlist)), 'k-', 'LineWidth', 2);
    hMarker = plot(0, 0, 'r', 'MarkerSize', 10, 'LineWidth', 2); 

    ylim([-5, 5]);
    xlim([-5 string_length]);
    xlabel('x');
    ylabel('displacement u(x,t)');
    title('Traveling Wave on a Discretized String (Triangle Pulse)');

    dT = 0.002;    
    next_t = 0;

    for k = 1:length(tlist)

        if tlist(k) < next_t
            continue
        end
        next_t = next_t + dT;

        Uf_val = Uf_func(tlist(k));
        if numel(Uf_val) > 1
            Uf_val = Uf_val(1);
        end

        y_now = [0; Ulist(k,:)'; Uf_val];

        set(hString, 'YData', y_now);

        % move marker at speed c (just a visualization reference)
        marker_x = string_length - c*tlist(k) + pulse_width/2;
        %short blurb showing how to find x-coord of tracking line
        %x = x-coord of tracking line, t = current time
        %c = wave speed, w = pulse width (in time), L = string length
        marker_x = mod(marker_x,2*string_length);
        if marker_x > string_length
        marker_x = 2*string_length - marker_x;
        end
        set(hMarker, 'XData', [marker_x,marker_x], 'YData', pulse_height*[-1,1]);
        % if marker_x >= 0 && marker_x <= string_length
            % set(hMarker, 'XData', marker_x, 'YData', 0);
        % end

        drawnow;
    end

end


% Triangle pulse function

function res = triangle_pulse(t,w,h)
    t = t*(2/w);
    res = 1 - min(abs(t-1),1);
    res = h*res;
end

function res = triangle_pulse_derivative(t,w,h)
    t = t*(2/w);
    res = -sign(t-1) .* (abs(t-1)<1);
    res = (2*h/w)*res;
end
