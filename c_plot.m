close all
clc,clear
t0=20000;%minutos
v_CO2=5.7e-5;%[m3] V=57+70 mL=127mL
n_CO2_0=0.0949;% moles in 57 mL at 23ºC and 4.1 Mpa 
Dia=0.0127;%[m]
L=v_CO2/(pi*Dia^2/4);%[m]
T=23+273;
R=8.3145;
    P0=3.23962;%[Pa]
    A1=0.04464;
    A2=0.32977; 
    A3=0.42264;
    tE0=8.75335;
    t1=13.41189;
    t2=137.49629;
    t3=901.37248;
t_E=0:t0/100:t0;      
P_lab=(P0+A1.*exp(-(t_E-tE0)./t1)+A2.*exp(-(t_E-tE0)./t2)+A3.*exp(-(t_E-tE0)./t3))*10e5; %plot(t_E(1:100),P_lab(1:100)) TE IN MINUTES
P_eq=3.24e6;%by looking at the plot
c_sat=6.1896e-4.*P_eq-0.255;% is it valid????
% Define ranges for KD and D
%KD_range = [100, 200];
% D_range = [5.4807e-8, 15e-9];% AT 40ºC 2.1559e-8

% Initialize variables
pressure_difference_min = 40000;%40kPa

c_over_csat=zeros(100,101);

%%
% Load the CSV data
data = readmatrix('dataset.csv');

% Extract time and diffusion coefficient columns
time = 1200.*data(:, 1);              % First column is time IN MINUTES           
diff=6e-8*data(:, 2); % Second column is diffusion coefficient

t_interp = linspace(min(time), max(time), 101);

% Interpolar el coeficiente utilizando 'interp1'
diff_interp = interp1(time, diff, t_interp, 'linear');

% plot(data(:, 1), data(:, 2), 'k-', 'LineWidth', 2)

% Create the matrix D_T
D_T = [t_interp; diff_interp];

%
% figure;
% plot(t_interp, diff_interp, 'k-', 'LineWidth', 2);
% xlabel('Time (min)','Interpreter', 'latex', 'FontSize', 14);
% ylabel('Diffusion Coefficient D (m^2/s)','Interpreter', 'latex', 'FontSize', 14);
% title('Diffusion Coefficient vs Time at 23°C','Interpreter', 'latex', 'FontSize', 16);
% grid on;
% grid minor; 
% set(gca, 'FontSize', 12); 
% set(gca, 'LineWidth', 0.5);
% set(gca, 'TickLabelInterpreter', 'latex');


%% Iterate over KD and D values
%KD = KD_range(1);
D = diff_interp;%[AT 23ºC]AT 40ºC= 2.1559e-8
        X1=1;
        for x=0:0.0045:L
            tt1=1;
            for t=0:t0*0.6:t0*60
                c_over_csat(X1,tt1) = solve_equation(x, t, L, D(tt1)); 
                tt1=tt1+1;
            end
            X1=X1+1;
        end
        c=c_over_csat.*c_sat;
 time2=0:t0/100:t0;       
 figure (1)
 plot (time2,c(30,:),time2,c(50,:),time2,c(75,:),time2,c(90,:),'k-', 'LineWidth', 2) 
 xlabel('Time (min)','Interpreter', 'latex', 'FontSize', 14);
 ylabel('Concentration (mol/m3)','Interpreter', 'latex', 'FontSize', 14)
 title('Concentration evolution at different positions','Interpreter', 'latex', 'FontSize', 16);
 legend ('Node 30','Node 50','Node 75','Node 90')
 grid on;
 grid minor; 
 set(gca, 'FontSize', 12); 
 set(gca, 'LineWidth', 0.5);
 set(gca, 'TickLabelInterpreter', 'latex');
 % Calculate  n, and Z
 %n= solve_equation_n(c, L, t0, n_CO2_0); %something is wrong

 % Get the dimensions of the concentration matrix
 [num_nodes, num_timesteps] = size(c);

 % Create time and node vectors
    time_values = linspace(0, 2e4, num_timesteps);
    node_values = linspace(0, 100, num_nodes);

 % Create a meshgrid for contour plotting
    [X, Y] = meshgrid(time_values, node_values);

 % Create the contour plot
 figure (2)
    contourf(X, Y, c, 1200,'LineStyle', 'none');

 % Set the colormap and colorbar
    colormap(jet);
    colorbar;
% Specify the desired constant values for the lines
    constant_values = [275.0, 550.0, 825.0, 1100, 1375, 1650];

    % Plot the constant-concentration lines
    hold on;
    contour(X, Y, c, constant_values, 'LineColor', 'black');
    hold off;

    % Set labels and title
    xlabel('Time(min)','Interpreter', 'latex', 'FontSize', 14);
    ylabel('Nodes','Interpreter', 'latex', 'FontSize', 14);
    title('CO2 Concentration @ 23°C (mol/m³)','Interpreter', 'latex', 'FontSize', 16);
% grid on;
% grid minor; 
set(gca, 'FontSize', 12); 
set(gca, 'LineWidth', 0.5);
set(gca, 'TickLabelInterpreter', 'latex');


max(max(c));
