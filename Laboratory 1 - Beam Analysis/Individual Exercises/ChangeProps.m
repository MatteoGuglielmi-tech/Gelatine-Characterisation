% Include k-Wave folder and subfolders to the current path
addpath(genpath('k-Wave'));

close all 
clear zaxis
clc

%Parameters initialization

% Source parameters
f0 = 2.0e6;                                                                   % center frequency (of the signal transmitted by the source) [Hz]
Bandwidth = 0.5e6;                                                            % -6 dB half-bandwidth (of the signal transmitted by the source) [Hz]

Ps = 1e6;                                                                     % source peak pressure (max pressure level of the pressure wave field at the source) [Pa]

% Acoustic properties of background medium
c0 = 1500;                                                                    % speed of sound [m/s]
rho0 = 1000;                                                                     % density of mass [kg/m^3]
lambda = c0/f0;                                                               % wave length [m]

% Define the grid size (map the physical domain into a numerical domain)
% We need to define the dimension of a pixel
dx = lambda/8;                                                                % grid size along x-dimension -> lambda/8 is a good approximation
dz = dx;                                                                      % grid size along z-dimension

Width = 0.02;                                                                 % computational domain physical size along x-dimension [m]
Lenght = 0.04;                                                                % computational domain physical size along z-dimension [m]

Nx = round(Width./dx);                                                        % computational domain size in pixels along x-dimension
Nz = round(Lenght./dz);                                                       % computational domain size in pixels along z-dimension

% grid initialization
kgrid = makeGrid(Nz, dx, Nx, dz);

% acoustic properties within the domain
medium.sound_speed              = c0*ones(Nz, Nx);                          % sound speed for each pixel[m/s]
medium.density                  = rho0*ones(Nz, Nx);                        % [kg/m^3]

% grid size check
cfl=0.1;

% simulation time [physical time]
% here we take a margin for simulation time in order to be sure
% the wave reaches the end of the space
t_end = (Width/10+sqrt(Lenght^2+(Width/2)^2))/c0;   % [s]

% grid (makeGrid creates the time array needed to run the simulation)
kgrid.t_array = makeTime(kgrid, max(medium.sound_speed(:)), cfl, t_end);


% %% ====================================================================================
% %% Plane Wave =========================================================================
% %% ====================================================================================
% % define source mask for a linear transducer with an odd number of elements   
% % transducers transmit the waves with no time delay.
Nelements = 40;                         
Pitch = 4;                                                                    % Distance between centers of two sensors[grid points]
kerf = 1;                                                                     % Distance between sensors[grid points]
num_elements = Nelements*Pitch;                                               % [grid points] Total aperture size
x_offset = 15;                                                                % The transducer is not initialized exactly on top because the transducer when excited propagates also backwards
                                                                              % [grid points] Needed for PML (Perfectly matched layer)

% Create the source mask  
source.p_mask = zeros(Nz, Nx);
start_index = round(Nx/2) - round(num_elements/2) + 1;                        % Middle of numerical domain along x-axis minus half of the physical aperture
source.p_mask(x_offset, start_index:start_index + num_elements - 1) = 1;
% Create the kerf 
for indexElement1=1:Nelements
    source.p_mask(x_offset,start_index + (indexElement1*Pitch)-(kerf-1) : start_index + indexElement1*Pitch) = 0;
end
source.p_mask(x_offset,start_index) = 0;

% The input signal to each element is then created using the function toneBurst with a geometrically steered temporal offset that varies across the source.
% Define the properties of the tone burst used to drive the transducer
dt = (kgrid.t_array(end)-kgrid.t_array(end-1));                                % Sampling interval [s]
sampling_freq = 1/dt;                                                          % Sampling frequency [Hz]
steering_angle = 0;                                                            % [deg]
element_spacing = dx;                                                          % [m]
tone_burst_freq = f0;                                                          % Center frequency of the pulse [Hz]
tone_burst_cycles = round(tone_burst_freq/Bandwidth); 

% create an element index relative to the centre element of the transducer
element_index = -(num_elements - 1)/2:(num_elements - 1)/2;

% use geometric beam forming to calculate the tone burst offsets for each transducer element based on the element index
tone_burst_offset = element_spacing*element_index*(sin(steering_angle*pi/180)./cos(steering_angle*pi/180));
tone_burst_offset = tone_burst_offset-min(tone_burst_offset);
tone_burst_offset = tone_burst_offset./(c0*dt); 

% create the tone burst signals
source.p = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles, 'SignalOffset', tone_burst_offset); % by default is a Gaussian Pulse
SourceSize = size(source.p);
source.p = (source.p./max(max(abs(source.p))));

% plots
figure(1)
subplot(2,1,1)
plot((0:dt:(SourceSize(2)-1)*dt)*1e6,source.p(1,:))
hold on
plot((0:dt:(SourceSize(2)-1)*dt)*1e6,abs(hilbert(source.p(1,:))))
axis([0 (SourceSize(2)-1)*dt*1e6 -1.2 1.2])
xlabel('Time [\mus]')
ylabel('Normalized amplitude')
title('simulated signal')
subplot(2,1,2)
plot((0:1/((SourceSize(2)-1)*dt):1/dt)./1e6,20*log10(abs(fft((source.p(1,:) ))))-max(20*log10(abs(fft((source.p(1,:) ))))))
axis([0 5 -20 5])
xlabel('Frequency [MHz]')
ylabel('Normalized amplitude [dB]')
title('simulated signal')

% compute the pulselength by counting the number of values of abs(hilbert(source.p(1,:)) greater than 0.5 
pulseLength = find(abs(hilbert(source.p(1,:)))>0.5);
% axial resolution
ar = pulseLength(end)-pulseLength(1);
ar = ar*dt*c0/2; % convert the axial resolution in meters
ar = ar*1e3;     % convert the axial resolution in millimeters

% creates kerf 
for indexElement1 = 0:Nelements-1
source2.p(1+(Pitch-kerf)*indexElement1:Pitch-kerf+(Pitch-kerf)*indexElement1,:) = source.p(1+(Pitch)*indexElement1:Pitch-kerf+(Pitch)*indexElement1,:);
end

source.p = source2.p;
clear source2


%% =========================================================================
% create a sensor mask covering the entire computational domain using the opposing corners of a rectangle
sensor.mask = [1, 1, Nz, Nx].';

% set the record mode to capture the final wave-field and the statistics at each sensor point 
sensor.record = {'p_final', 'p_max', 'p_rms'};

% create a display mask to display the transducer
display_mask = source.p_mask;

% assign the input options
input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false,'RecordMovie', true, 'MovieName', 'PlaneWave'};
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% generate axis
xaxis = (0:Nx-1)*dx;
xaxis = xaxis-mean(xaxis);
zaxis = (0:Nz-24)*dz;

% clear data
sensor_data.p_max = sensor_data.p_max(x_offset:end,:);
sensor_data.p_max = sensor_data.p_max(1:end-x_offset,:);

% plots
figure(2)
imagesc(xaxis*1e3,zaxis*1e3,20*log10(abs(sensor_data.p_max))-max(max(20*log10(abs(sensor_data.p_max)))))
clim([-20 0])
axis image
xlabel('x-axis [mm]')
ylabel('z-axis [mm]')
title('Plane Wave')

PlaneWaveField = sensor_data.p_max*Ps;

%% =================================================================================
%% Focused =========================================================================
%% =================================================================================
% define source mask for a linear transducer with an odd number of elements   
Nelements = 40;                         
Pitch = 4;                                                                    % [grid points]
kerf = 1;                                                                     % [grid points]
num_elements = Nelements*Pitch;                                               % [grid points] Total aperture size
x_offset = 15;                                                                % [grid points] Needed for PML

source.p_mask = zeros(Nz, Nx);
start_index = round(Nx/2) - round(num_elements/2) + 1;
source.p_mask(x_offset, start_index:start_index + num_elements - 1) = 1;
for indexElement1=1:Nelements
source.p_mask(x_offset,start_index + (indexElement1*Pitch)-(kerf-1) :start_index + indexElement1*Pitch) = 0;
end
source.p_mask(x_offset,start_index) = 0;

% The input signal to each element is then created using the function toneBurst with a geometrically steered temporal offset that varies across the source.
% define the properties of the tone burst used to drive the transducer

dt = (kgrid.t_array(end)-kgrid.t_array(end-1));                                % [s]
sampling_freq = 1/dt;                                                          % [Hz]
steering_angle = 0;                                                            % [deg]
element_spacing = dx;                                                          % [m]
tone_burst_freq = f0;                                                          % [Hz]
tone_burst_cycles = round(tone_burst_freq/Bandwidth); 

% create an element index relative to the centre element of the transducer
element_index = -(num_elements - 1)/2:(num_elements - 1)/2;

% use geometric beam forming to calculate the tone burst offsets for each transducer element based on the element index
tone_burst_offset = element_spacing*element_index*(sin(steering_angle*pi/180)./cos(steering_angle*pi/180));


%% Focusing
FocalDepth = 0.020;                                                            % [m]
Radius = sqrt(FocalDepth^2+(dx*num_elements/2)^2);
alpha = asin(((0:(num_elements)/2)*dx)./Radius);
A = Radius*cos(alpha)-FocalDepth;
dR = A./(c0*dt);
tone_burst_offset((num_elements)/2:end) = dR;
tone_burst_offset(1:(num_elements)/2) = fliplr(dR(1:end-1));

% create the tone burst signals
source.p = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles, 'SignalOffset', tone_burst_offset);
SourceSize = size(source.p);
source.p = (source.p./max(max(abs(source.p))));

% creates kerf 
for indexElement1=0:Nelements-1
source2.p(1+(Pitch-kerf)*indexElement1:Pitch-kerf+(Pitch-kerf)*indexElement1,:) = source.p(1+(Pitch)*indexElement1:Pitch-kerf+(Pitch)*indexElement1,:);
end

source.p=source2.p;
clear source2


%% =========================================================================
% create a sensor mask covering the entire computational domain using the opposing corners of a rectangle
sensor.mask = [1, 1, Nz, Nx].';

% set the record mode to capture the final wave-field and the statistics at each sensor point 
sensor.record = {'p_final', 'p_max', 'p_rms'};

% create a display mask to display the transducer
display_mask = source.p_mask;

% assign the input options
input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false,'RecordMovie', true, 'MovieName', 'Focused'};
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% generate axis
xaxis = (0:Nx-1)*dx;
xaxis = xaxis-mean(xaxis);
zaxis = (0:Nz-24)*dz;

% clear data
sensor_data.p_max = sensor_data.p_max(x_offset:end,:);
sensor_data.p_max = sensor_data.p_max(1:end-x_offset,:);

% plots
figure(3)
imagesc(xaxis*1e3,zaxis*1e3,20*log10(abs(sensor_data.p_max))-max(max(20*log10(abs(sensor_data.p_max)))))
clim([-20 0])
axis image
xlabel('x-axis [mm]')
ylabel('z-axis [mm]')
title('Focused')

FocusedField = sensor_data.p_max*Ps;
