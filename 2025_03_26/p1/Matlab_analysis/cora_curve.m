clear; clc;

% Parametros

p.g = 0.01;
p.gU = 0.05;
p.gW = 0.0001;
p.e0 = 0.0001;
p.eM = 0.5;

p.mU = 1.3261942096116295;
p.mW = 1.830965707922456;
p.eP = 23.259110949433275;

p.gY = 0.1;
p.mY = 1;
p.mY_p = p.mY*1.05;
%p.mY_p = 2;

p.mV = p.mY;
p.mV_p = p.mY_p;
p.gV = p.gY;

p.tmax = 5000;
points = 100000;

%% simulation of ATF
% Simulate ODE dynamics:
n = 20;
coras = zeros(1,n);
values = logspace(log10(0.001), log10(10), n);
y0 = [0, 0, 0, 0]; % Initial conditions 

for i = 1:length(values)
    p.mY = values(i);
    p.mY_p = p.mY*1.05;

    tspan = linspace(0, p.tmax, points);
    % Usa una función anónima para pasar la estructura de parámetros a ode45
    [t, y] = ode45(@(t, y) odefun(t, y, p), tspan, y0);
    
    
    % % Plot ODE results
    figure()
    plot(t, y(:, 1), '-', 'LineWidth', 2.5); hold on;
    plot(t, y(:, 2), '-', 'LineWidth', 2.5);
    plot(t, y(:, 3), '-', 'LineWidth', 2.5);
    plot(t, y(:, 4), '-', 'LineWidth', 2.5);
    hold off;

    ylim([0 30])
    legend({"W", "Y", "U", "C"})
    xlabel("time")
    ylabel("molecules")
    title('Plot 1')
    
    SS = y((points/2)-1,:);     % array of steady states before perturbation
    SS_p = y(points,:);   % array of steady states after perturbation 
    Yss = SS(2);
    
    % calculate cora
    Y_f = SS(2);
    Y_f_p = SS_p(2);
    
    % Constant imput for analogous system
    p.mUs = p.mU*Yss;
    
    %% Simulate ODE dynamics for analog system:
    ciA = SS; % Initial conditions 
    y0 = SS; 

    tspan = linspace(0, p.tmax, points);
    [t, y] = ode45(@(t, y) odefun2(t, y, p), tspan, ciA);
    
    
    %SS = y((points/2)-1,:);     % array of steady states before perturbation
    SS_p = y(points,:);   % array of steady states after perturbation 
    
    Y_nf = SS(2);
    Y_nf_p = SS_p(2);
    
    
    % % % Plot ODE results
    figure()
    plot(t, y(:, 1), '-', 'LineWidth', 2.5); hold on;
    plot(t, y(:, 2), '-', 'LineWidth', 2.5);
    plot(t, y(:, 3), '-', 'LineWidth', 2.5);
    plot(t, y(:, 4), '-', 'LineWidth', 2.5);hold off;
    legend({"W","Y","U","C"})
    xlabel("time")
    ylabel("molecules")
    title('Plot 1')
    
    
    coras(i) = log10(Y_f_p/Y_f) / log10(Y_nf_p/Y_nf);
end
figure()
plot(coras)
ylim([0 1])



figure()
plot(values, coras)
set(gca, 'XScale', 'log');
ylim([0 1])

function dydt = odefun(t, y, p)
    W = y(1);
    Y = y(2);
    U = y(3);
    C = y(4);

    if (t > p.tmax/2)
       p.mY = p.mY_p;
    end

     dY = (p.mY * W) - ((p.g + p.gY) * Y);
    dU = (p.mU * Y) - ((p.g + p.gU) * U) - (p.eP * U * W) + ((p.e0 + p.gW + p.eM) * C);
    dW =    p.mW    - ((p.g + p.gW) * W) - (p.eP * U * W) + ((p.e0 + p.gU) * C);
    dC =              - ((p.g + p.eM) * C) + (p.eP * U * W) - ((p.gU + p.gW + p.e0) * C);

    dydt = [dW; dY; dU; dC];
end


function dydt = odefun2(t, y, p)
    W = y(1);
    Y = y(2);
    U = y(3);
    C = y(4);

   
    p.mY = p.mY_p;
    

    % ODEs:   

    dY = (p.mY * W) - ((p.g + p.gY) * Y);
    dU = p.mUs - ((p.g + p.gU) * U) - (p.eP * U * W) + ((p.e0 + p.gW + p.eM) * C);
    dW =    p.mW   - ((p.g + p.gW) * W) - (p.eP * U * W) + ((p.e0 + p.gU) * C);
    dC =            - ((p.g + p.eM) * C) + (p.eP * U * W) - ((p.gU + p.gW + p.e0) * C);


    dydt = [dW; dY; dU; dC];
end