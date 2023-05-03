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

seq.plot()

# export k-space trajectory
k_traj_adc, k_traj, t_excitation, t_refocusing, t_adc = seq.calculate_kspacePP()

kx = k_traj_adc[0,:]
ky = k_traj_adc[1,:]

# solve for density compensation weights. Assume it is the same each repetition for now.
# this dcf computation is not really great, but it is a start. 
# TODO: use voronoi density compensation instead.
Nsample = int(k_traj_adc.shape[1]/repetitions)
w = mri.dcf.pipe_menon_dcf(k_traj_adc[:,0:Nsample+1].T)
w = w / (np.max(w));
w[w > 0.6] = 0.6;
w = w / np.max(w);   
w[int(w.shape[0]*2/3):w.shape[0]] = 1

apd.app_params['repetitions'] = repetitions

savemat("/server/home/pkumar/pulseq_metadata/test.mat", {'kx': kx, 'ky': ky, 'w' : w, 'param': apd.app_params})