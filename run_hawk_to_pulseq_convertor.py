from hawk_to_pulseq import hawk_to_pulseq
import numpy as np
from scipy.io import savemat
import sigpy.mri as mri

application_path = "Sequences/RTHawk/RealTimeSpiralSSFP055_FOV24_n13/"

# convert the sequence. 
# note that the hardware settings must be encoded in the .spv files, not specified by pypulseq directly.
seq, sys, repetitions, apd = hawk_to_pulseq(application_path)

# write the sequence file
seq.write("Sequences/Pulseq_Export/" + application_path.split("/")[-2] + '.seq')

# seq.plot()

# export k-space trajectory
k_traj_adc, k_traj, t_excitation, t_refocusing, t_adc = seq.calculate_kspacePP()

Nsample = int(k_traj_adc.shape[1]/repetitions)
kx = k_traj_adc[0,:]
ky = k_traj_adc[1,:]
k_max = np.max(np.abs(kx + 1j * ky))
k = (kx / k_max) + (1j * ky / k_max)
# w = dcf_voronoi(k)

# calculate density compensation weights using Pipe and Menon's method
# This isn't working so well... need to figure out why.
Nsample = int(k_traj_adc.shape[1]/repetitions)
w = mri.dcf.pipe_menon_dcf(k_traj_adc.T)
w = w[Nsample+1:2*Nsample+1]
w = w / (np.max(w));
w[w > 0.4] = 0.4;
w = w / np.max(w);   
w[int(w.shape[0]*2/3):w.shape[0]] = 1

apd.app_params['repetitions'] = repetitions

savemat("ismrmrd-pypulseq-matlab-server/pulseq_metadata/" + application_path.split("/")[-2] + ".mat", {'kx': kx, 'ky': ky, 'w' : w, 'param': apd.app_params})