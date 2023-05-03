clear; close all; clc;

addpath(genpath('~/repos/pulseq/'))
addpath(genpath('thirdparty/'));
addpath(genpath('utils/'));


%% Scanner Params
filename = "RealTimeSpiralSSFP055_FOV24_n13";

[RF_phase, RF_amp, Gx, Gy, Gz, param] = loadFromSpinBench(filename + ".spv", filename + ".wfm");

% plot_waveforms_publish(RF_amp, Gx, Gy, Gz, param);

% HACK for view order, only for this speech scan. It is bitreversed order. 
view_order = [0, 8, 4, 12, 2, 10, 6, 1, 9, 5, 3, 11, 7] + 1;
Gx = Gx(:, view_order);
Gy = Gy(:, view_order);
Gz = Gz(:, view_order);
RF_amp = RF_amp(:, view_order);
RF_phase = RF_phase(:, view_order);


Gmax = param.gradLimit * 10 + 0.1; % convert from G/cm to mT/m. add arbitrary 0.1 for passing checks. 
Smax = param.slewRateLimit * 10; % same conversion as above.
B0 = param.B0;

if param.gradSamplingRate ~= 100
    error("grad sampling rate in spv file is incorrect! pulseq file will not convert.");
end

rfRasterTime = 10e-6 / param.rfRateMultiplier;
gradRasterTime = 10e-6;
RFDeadTime = 100e-6;
dt_adc = 1e-6;

sys = mr.opts('MaxGrad',Gmax,'GradUnit','mT/m',...
    'MaxSlew',Smax,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6, 'rfDeadTime', RFDeadTime, ...
    'adcDeadTime', 10e-6, 'gradRasterTime', gradRasterTime, 'rfRasterTime', rfRasterTime, 'B0', param.B0);

seq=mr.Sequence(sys);           % Create a new sequence object

% convert units. Crucially important for the units to be correct before
% pulseq.
% kHz/m
maxgrad_sb = param.gradLimit * 1e-2; % grad limit is G/cm -> 
Gx = Gx * (1 / 32767) * (maxgrad_sb) * 42.58 * (1e6);
Gy = Gy * (1 / 32767) * (maxgrad_sb) * 42.58 * (1e6);
Gz = Gz * (1 / 32767) * (maxgrad_sb) * 42.58 * (1e6);

% add delay to all gradients due to RF dead time.
Gx = [zeros(floor(RFDeadTime/gradRasterTime), size(Gx,2)); Gx;];
Gy = [zeros(floor(RFDeadTime/gradRasterTime), size(Gy,2)); Gy;];
Gz = [zeros(floor(RFDeadTime/gradRasterTime), size(Gz,2)); Gz;];

% truncate Gx and Gy and Gz amplitudes.
extra_samples = 10;
idxGx = find(Gx(:,1), 1)-extra_samples;
idxGy = find(Gy(:,1), 1)-extra_samples;
idxGz = find(Gz(:,1), 1, 'last')+extra_samples;

delayGx = idxGx * gradRasterTime;
delayGy = idxGy * gradRasterTime;

Gx = Gx(idxGx:end,:);
Gy = Gy(idxGy:end,:);
Gz = Gz(1:idxGz,:);

% truncate RF amplitudes.
RF_trunc_end = find(RF_amp(:,1), 1, 'last')+extra_samples;
% RF_trunc_start = find(RF_amp(:,1), 1, 'first')-10;
% RF_amp((1:RF_trunc_start), :) = [];
RF_amp(RF_trunc_end:end, :) = [];
% RF_delay = RF_trunc_start * dt;
Nint = size(Gx, 2);

%% Spiral Readout
gx = {};
gy = {};
gz = {};
rf = {};
adc = {};

adc_discard = 0;
for i = 1:Nint
    gx{i} = mr.makeArbitraryGrad('x', Gx(:,i)', 'delay', delayGx, 'system', sys);
    gy{i} = mr.makeArbitraryGrad('y', Gy(:,i)', 'delay', delayGy, 'system', sys);
    gz{i} = mr.makeArbitraryGrad('z', Gz(:, i)', sys);
    rf{i} = mr.makeArbitraryRf(RF_amp(:,i)', (param.FA*pi/180), 'system', sys);
    
    % adc_start = find((Gx(:,i) .^2 + Gy(:, i) .^2 > 0), 1, 'first');
    
    adc_start = min(idxGx, idxGy) + extra_samples;
    
    % adc_last = find((Gx(:,i) .^2 + Gy(:, i) .^2 > 0), 1, 'last');
    % adc_sb{i} = mr.makeAdc(floor((adc_last-adc_start) * (dt/dt_adc)), 'Dwell', dt_adc, 'Delay', (adc_start) * dt, 'system', sys);

    adc{i} = mr.makeAdc(floor((param.readoutDuration * 1e-3) * (1/dt_adc)) + adc_discard, 'Dwell', dt_adc, 'Delay', (adc_start * gradRasterTime) - (adc_discard * dt_adc), 'system', sys);
end

% Loop over phase encodes and define sequence blocks

% concatenate number of reps to repeat.
% we are doing 2 loops because this guarantees SSFP flipping
% phase on even RF increments.
repeat_loops = 2;

for j=1:(repeat_loops * Nint)
    i = mod(j, Nint);
    if i == 0
        i = Nint;
    end
    rf{i}.phaseOffset=pi*mod(j,2);
    adc{i}.phaseOffset=pi*mod(j,2);
    % seq.addBlock(mr.makeLabel('SET','LIN', i));
    seq.addBlock(rf{i},gz{i},gx{i},gy{i}, adc{i});
end

fprintf('Sequence ready\n');

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% prepare export
seq.setDefinition('FOV', [param.fov*1e-2 param.fov*1e-2 param.thickness*1e-3]);
seq.setDefinition('Name', 'trufi');

seq.write(filename + '.seq')       % Write to pulseq file

%seq.install('siemens');
%% plots and displays

seq.plot('timeDisp','us');

%% plot entire interpolated wave forms -- good for debugging of gaps, jumps,
% etc, but is relatively slow
%gw=seq.gradient_waveforms();
%figure; plot(gw'); % plot the entire gradient waveform

%{
gw_data=seq.waveforms_and_times();
figure; plot(gw_data{1}(1,:),gw_data{1}(2,:),gw_data{2}(1,:),gw_data{2}(2,:),gw_data{3}(1,:),gw_data{3}(2,:)); % plot the entire gradient waveform
title('gradient waveforms');
%}

%% trajectory calculation 
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

ktraj_adc = ktraj_adc(:, 1:(size(ktraj_adc,2)/repeat_loops));

% deal with the discarded adc samples.
ktraj_adc = reshape(ktraj_adc, 3, adc{1}.numSamples, param.repetitions);
ktraj_adc = ktraj_adc(:, adc_discard+1:end, :);
ktraj_adc = reshape(ktraj_adc, 3, []);

Nsamp = adc{1}.numSamples - (adc_discard);

% plot k-spaces
%{
figure; plot(t_ktraj,ktraj'); title('k-space components as functions of time'); % plot the entire k-space trajectory
figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
title('2D k-space');
%}

% maybe flipping needs to hppen, but could be not.
% ktraj_adc(2,:)=-ktraj_adc(2,:);

kx = ktraj_adc(1,:);
ky = ktraj_adc(2,:);
Nlines = param.repetitions;

k_max = max(sqrt(kx(:).^2 + ky(:).^2));
k_ = kx(:) / k_max + 1i * ky(:) / k_max;
w = voronoidens(k_);
w = [w; zeros(Nsamp-length(w),1)];
w = w(1:Nsamp);
w = w / (max(w));
w(w > 0.5) = 0.5;
w = w / max(w);

save("/server/home/pkumar/pulseq_metadata/" + filename + ".mat", "kx", "ky", "w", "param");

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  

rep = seq.testReport;
fprintf([rep{:}]);


