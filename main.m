% Coded by Panqi Chen, 2023.
% This file test VSD-fort for the simulated data and plot figure 1 in the
% paper entitled: "Estimating Channels With Hundreds of Sub-Paths for
% MU-MIMO Uplink: A Structured High-Rank Tensor Approach";

clc;clear all;close all; 
addpath(genpath(pwd));
% generate simulated tensor
% determine the size and the rank of the simulated tensor.
size_tens = [16 32 4 2]; R = 420; 
U = cpd_rnd(size_tens,R,'Imag',@randn);
U{1} = struct_vander(exp(120*rand(R,1)*1j),[],size_tens(1)-1,true);
U{2} = struct_vander(exp(500*rand(R,1)*1j),[],size_tens(2)-1,true);
U{3} = struct_vander(exp(15*rand(R,1)*1j),[],size_tens(3)-1,true);
T = cpdgen(U);

% specify the hyperparameters
K1=11;K2=21;K3=2;
L1=size_tens(1)+1-K1;
L2=size_tens(2)+1-K2;
L3=size_tens(3)+1-K3;
% maximum paths that can be uniquely identified
R_max = min((K1-1)*K2*K3,L1*L2*L3*2);



% run VSDfort
S=[K1,K2,K3];
[Utmp,d,err] = VSDfort(T,R,S);
      
            
        
         
        
%% plot figure
b1=U{1}(2,:);b2=U{2}(2,:);b3=U{3}(2,:);d1=Utmp{1}(2,:);d2=Utmp{2}(2,:);d3=Utmp{3}(2,:);
figure(1);

markersize=9;
fontsize=33;

set(gcf,'Position',[100 0 1000 700]);
plot(real(b1), imag(b1), 'ko', 'LineWidth', 0.5,'MarkerSize', markersize);
hold on;
plot(real(d1), imag(d1), 'r+', 'LineWidth', 0.7,'MarkerSize', markersize);
hold on;
axis([-1.2 1.2 -1.5 1.2]);
legend('True','Estimated','fontname','Times New Roman','fontsize',fontsize);
xlabel('Real'), ylabel('Imag');
set(gca,'fontsize',fontsize,'fontname','Times');
% Define the area of the plot you want to magnify
x_zoom = [-0.98,-0.83];
y_zoom = [0.36,0.49];
% Highlight the zoomed area 5ith a rectangle on the main plot
rectangle('Position', [x_zoom(1), y_zoom(1), diff(x_zoom), diff(y_zoom)], ...
          'EdgeColor', 'b', 'LineWidth', 2);

% Create the inset plot in the right bottom corner
inset_axes = axes('Position', [.554 .17 .35 .35]);
box on; % Put a box around the inset
plot(real(b1), imag(b1), 'ko', 'LineWidth', 0.7,'MarkerSize', markersize+5);
hold on;
plot(real(d1), imag(d1), 'r+', 'LineWidth', 1,'MarkerSize', markersize+5);
hold on;
xlim(inset_axes, x_zoom);
ylim(inset_axes, y_zoom);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
ax = gca; % Gets the current axes
ax.XColor = 'blue'; % Changes the x-axis border color to red
ax.YColor = 'blue'; % Changes the y-axis border color to red
ax.LineWidth = 2;
set(gca,'fontsize',fontsize,'fontname','Times');



 
 