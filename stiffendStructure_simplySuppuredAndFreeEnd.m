clear all
close all
clc

AREA_BASELINE = 766.6227;

%% Material properties and design choices

% properties in MPa
E1 = 62046;
E2 = 62046;
G12 = 4826;

Xt= 1517;
Xc= 1379;
Yt= 1450;
Yc= 1379;
S = 99;

nu12 = .05;
t = .1905; % in mm

% dimensions of the panel in mm
a = 1000;   % length
b = 480;    % width


stiffners_on_edge = true;
num_beams = 7;
if stiffners_on_edge
    num_skin_sections = num_beams - 1;
else
    num_skin_sections = num_beams + 1;
end

theta_skin = [45, 45, 45, 45];    % -> [(0/90), (+-45)]s
theta_beam = [0, 0, 0];      % -> [(0/90), (0/90)]s

h_skin = t * length(theta_skin);
h_beam = t * length(theta_beam);

%% Geometry of single beam-skin section

% defining dimensions of the cross-setion of the beam
% alpha_beam = acosd(2.5/15);    % angle between angled beam part and skin (in deg)
% a1 = 15;            % outer legs of the beam
% a2 = 20;            % upper horizontal part
% a3 = sqrt(15^2-2.5^2);             % height of beam in terms of the mid-planes
% L = 2*a1 + 25;

alpha_beam = 85;
a1 = 5;
a2 = 10;
a3 = 17;
L = 2*a1 + 2*a3/tand(alpha_beam) + a2;
d_s = (b-L)/num_skin_sections;

% geometric difinition of skin and beam cross-section
% row indeces for geometric definitions:
%       1: skin left
%       2: skin right
%       3: flange skin left
%       4: flange skin center
%       5: flange skin right
%       6: web left
%       7: web right
%       8: flange up
OEF = [false, false, true, false, true, false, false, false];

y_geom = zeros(8,4);
z_geom = zeros(8,4);

% skin left
y_geom(1,:) = [-d_s/2, 0, 0, -d_s/2];
z_geom(1,:) = [0, 0, h_skin, h_skin];

% skin right
y_geom(2,:) = [0, d_s/2, d_s/2, 0];
z_geom(2,:) = z_geom(1,:);

% flange skin left
z_geom(3,:) = h_skin + [0, 0, h_beam, h_beam];
y_geom(3,:) = - L/2 + [0, a1, a1, 0];

% flange skin center
z_geom(4,:) = h_skin + [0, 0, h_beam, h_beam];
y_geom(4,:) = (-L/2 + a1) * [1, -1, -1, 1];

% flange skin right
z_geom(5,:) = h_skin + [0, 0, h_beam, h_beam];
y_geom(5,:) = L/2 - [0, a1, a1, 0];

% web left
z_geom(6,:) = [h_skin+h_beam,...
               h_skin+h_beam+a3,...
               h_skin+a3,...
               h_skin+h_beam];
y_geom(6,:) = - L/2 + a1 +... 
             [h_beam/2 * cosd(alpha_beam+90) + (z_geom(6,1) - (h_skin + h_beam/2) - h_beam/2*sind(alpha_beam+90))/tand(alpha_beam),...
              h_beam/2 * cosd(alpha_beam+90) + (z_geom(6,2) - (h_skin + h_beam/2) - h_beam/2*sind(alpha_beam+90))/tand(alpha_beam),...
              h_beam/2 * cosd(alpha_beam-90) + (z_geom(6,3) - (h_skin + h_beam/2) - h_beam/2*sind(alpha_beam-90))/tand(alpha_beam),...
              h_beam/2 * cosd(alpha_beam-90) + (z_geom(6,4) - (h_skin + h_beam/2) - h_beam/2*sind(alpha_beam-90))/tand(alpha_beam)];

% web right
y_geom(7,:) =- y_geom(6,:);
z_geom(7,:) = z_geom(6,:);

% flange up
z_geom(8,:) = z_geom(6, [2,2,3,3]);
y_geom(8,:) = [1, -1, -1, 1] .* (y_geom(6, [2, 2, 3, 3]));




% Areas of different structural components in mm²
areas = polyarea(y_geom, z_geom, 2);
%areas(4) = 0;

% Centroids of the sections
yc = zeros(8,1);
zc = zeros(8,1);
for i = 1:8
    [yc(i), zc(i)] = centroid(polyshape(y_geom(i,:), z_geom(i,:)));
end


% plot beam and skin cross-section and top view
figure(1)
subplot(2,1,1)
hold on
for i = 1:8
    if areas(i) ~= 0
        plot(y_geom(i,[1,2,3,4,1]), z_geom(i,[1,2,3,4,1]), 'k')
    end
end

% plot(d_s/2 + [-L/2, L/2], (h_skin + h_beam/2) * ones(1,2), 'r:')
% plot(d_s/2 - L/2 + a1 + [0, a3/tand(alpha_beam)], h_skin + h_beam/2 + [0, a3], 'r:')
% plot(d_s/2 + L/2 - a1 - [0, a3/tand(alpha_beam)], h_skin + h_beam/2 + [0, a3], 'r:')
% plot(d_s/2 + (L/2 - a1 - a3/tand(alpha_beam)) * [1, -1], (h_skin + h_beam/2 + a3) * ones(1,2), 'r:')

plot(yc, zc, 'bh')

axis equal
title(sprintf("Cross-section of skin and beam structure for one of %d beams", num_beams))
xlabel('y-axis [mm]')
ylabel('z-direction [mm]')

subplot(2,1,2)
hold on
plot([0, a, a, 0, 0], [0, 0, b, b, 0], 'k')
for i=1:num_beams
    if stiffners_on_edge
        yc_beam = L/2 + (i-1) * (b-L)/(num_beams-1);
        fill([0, a, a, 0, 0], yc_beam+L/2*[-1, -1, 1, 1, -1], 'r')
    else
        yc_beam = ((i-1) + .5) * b/num_beams;
        fill([0, a, a, 0, 0], yc_beam+L/2*[-1, -1, 1, 1, -1], 'r')
    end
end
axis equal
title(sprintf("Top view of the whole panel for %d beams", num_beams))
legend('panel', 'stiffeners', 'location', 'northeastoutside')
xlabel('x-axis [mm]')
ylabel('y-axis [mm]')


%% Elastic properties of skin and beam

Q_beam = Q_th(E1, E2, nu12, G12, theta_beam);
D_beam = D_Qt(Q_beam,t);
D_beam_prime = inv(D_beam);
A_beam = A_Qt(Q_beam,t);
A_beam_prime = inv(A_beam);


Q_skin = Q_th(E1, E2, nu12, G12, theta_skin);
D_skin = D_Qt(Q_skin,t);
D_skin_prime = inv(D_skin);
A_skin = A_Qt(Q_skin,t);
A_skin_prime = inv(A_skin);




%% Properties of skin and beam components
% neutral axis membrane
Em_skin = 1/(h_skin * A_skin_prime(1,1));
Em_beam = 1/(h_beam * A_beam_prime(1,1));

Ems = [Em_skin * ones(1,2), Em_beam * ones(1,6)]';

z_neutral_membrane = sum(Ems.*areas.*zc)/sum(Ems.*areas);
figure(1)
subplot(2,1,1)
plot(d_s/2 * [-1, 1], z_neutral_membrane*ones(1,2), 'k--')


% equivalent EA
EAs = Ems.*areas;
EA_eq = sum(EAs);
EA_beam_eq = sum(Ems(3:end).*areas(3:end));


% neutral axis bending
Eb_skin = 12 / (h_skin^3 * D_skin_prime(1,1));
Eb_beam = 12 / (h_beam^3 * D_beam_prime(1,1));

Ebs = [Eb_skin * ones(1,2), Eb_beam * ones(1,6)]';

z_neutral_bending = sum(Ebs(3:end).*areas(3:end).*zc(3:end))/sum(Ebs(3:end).*areas(3:end));

% equivalent EI
factor_EI = [a1*h_beam^3; ...
             (L-2*a1)*h_beam^3; ...
             a1*h_beam^3; ...
             areas(6)*a3^2; ...
             areas(7)*a3^2; ...
             a2*h_beam^3]/12 ...
            + areas(3:end) .* (zc(3:end) - z_neutral_bending).^2;
EI_beam_eq = sum(Ebs(3:end) .* factor_EI);


% equivalent properties of the whole structure
A_eq = A_skin;
A_eq(1,1) = A_eq(1,1) + EA_beam_eq/d_s;

GJ = 0;
D_eq = D_skin;
D_eq(1,1) = D_eq(1,1) + EI_beam_eq * num_beams/b;
E_eq(3,3) = D_eq(3,3) + GJ/2/d_s;

%% Failure analysis
F_tot = 7e3;


% buckling of skin alone between stiffeners
AR = a/d_s;
%N0_skin = @(m) pi^2*(D_skin(1,1)*m.^4 + 2*(D_skin(1,2)+2*D_skin(3,3)) * m.^2 * AR^2 + D_skin(2,2) * AR^4) / (a^2 * m.^2);
% m = 1:3;
% N0_m = N0_skin(m);
% Nx_skin_buckling = N0_skin(m(N0_m == min(N0_m)));

A_param = A_skin(1,1) - A_skin(1,2)^2/A_skin(2,2);

F_tot_skin = @(k) (A_param + EA_beam_eq * num_beams/b)/A_param * pi^2*b/a^2 * (D_skin(1,1)*k.^4 + 2*(D_skin(1,2)+2*D_skin(3,3)) * k.^2 * AR^2 + D_skin(2,2) * AR^4) ./ k.^2;

if F_tot >= min(F_tot_skin(1:10))
    fprintf("Skin buckels! Applied force: %.2f, failure force: %.2f\n", F_tot, min(F_tot_skin(1:10)))
else
    fprintf("Skin buckling failure index: %.3f\n", F_tot/min(F_tot_skin(1:10)))
end

% bucking of the whole panel
AR = a/b;
F_tot_panel = @(m) pi^2*b/a^2 * (D_eq(1,1)*m.^4 + 2*(D_eq(1,2)+2*D_eq(3,3)) * m.^2 * AR^2 + D_eq(2,2) * AR^4) ./ m.^2;

if F_tot >= min(F_tot_panel(1:10))
    %disp('Whole panel buckels!')
    fprintf("Whole panel buckels! Applied force: %.2f, failure force: %.2f\n", F_tot, min(F_tot_panel(1:10)))
else
    fprintf("Panel buckling failure index: %.3f\n", F_tot/min(F_tot_panel(1:10)))
end


% cripling of beams
F_skin = A_param/(A_param + EA_beam_eq * num_beams/b) * F_tot;
F_beam = (F_tot - F_skin) / num_beams;

F_i_beam = EAs(3:end)/EA_beam_eq * F_beam;
sigma_i_beam = F_i_beam ./ areas(3:end);

b_t = areas(3:end) / h_beam^2;
sig_crip_cu = zeros(length(b_t),1);
sig_crip_cu(OEF(3:end) == true) = 1.63 ./ (b_t(OEF(3:end) == true)).^0.717;
sig_crip_cu(OEF(3:end) == false) = 11 ./ (b_t(OEF(3:end) == false)).^1.124;

% calculate first-ply failure stress for sig_cu
eps0_unitLoadNx = A_beam \ [-1; 0; 0];

sig_local = zeros(3,length(theta_beam));
failure_index = zeros(5,length(theta_beam));
for i = 1:length(theta_beam)
    sig_local(:,i) = rotate_stress(Q_beam(:,:,i) * eps0_unitLoadNx, theta_beam(i));
    failure_index(1:3,i) = sig_local(:,i) ./ [-Xc, -Yc, sign(sig_local(3,i))*S]';
    failure_index(4:5,i) = sig_local(1:2,i) ./ [Xt, Yt]';
end

% Compressive load of 1 N/mm² achieves this max. failure index
% This is used to calculate the load for which FFP occurs
max_FI = max(failure_index, [], 'all'); 
sig_cu = 1/max_FI/h_beam;


% failure stress
sig_fail = sig_crip_cu * sig_cu;

if any(sigma_i_beam > sig_fail)
    disp('Cripling')
else
    fprintf("Cripling failure index: %.3f\n", max(sigma_i_beam./sig_fail))
end


% FPF of skin
eps0_skin = A_skin\[-F_skin/b; 0; 0];

sig_local = zeros(3,length(theta_beam));
failure_index = zeros(5,length(theta_beam));
for i = 1:length(theta_beam)
    sig_local(:,i) = rotate_stress(Q_beam(:,:,i) * eps0_unitLoadNx, theta_beam(i));
    failure_index(1:3,i) = sig_local(:,i) ./ [-Xc, -Yc, sign(sig_local(3,i))*S]';
    failure_index(4:5,i) = sig_local(1:2,i) ./ [Xt, Yt]';
end

if max(failure_index, [], 'all') > 1
    disp('FPF in the skin')
else
    fprintf("Skin FPF failure index: %.3f\n", max(failure_index, [], 'all'))
end

%%
fprintf('Area ratio to baseline: %.3f\n', (sum(areas(3:end))*num_beams + b*h_beam) / AREA_BASELINE)
