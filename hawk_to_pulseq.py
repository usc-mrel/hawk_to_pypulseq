import extract_sb as sb
import numpy as np
import pypulseq as pp
import os

EXTRA_PAD = 10

class APD:
    def __init__(self):
        # an APD file contains a list of sections. 
        # the application section is the list of lists of lines that start with [Application]
        # the sequencing sections are the list of lists of lines that do not start with [Application] (likely sequence blocks)
        # len(SPVS) == len(loop_indices) ==  len(loop_block_max_indices) == len(linear_phase_increments) ==         # == len(sequencing_sections)
        # len(quadratic_phase_increments) == len(sequencing_sections).
        # for any lists that do not apply (for example waveform parameters for loop blocks), they are None.

        self.sequencing_sections = []
        self.application_section = []
        self.SPVS = []
        self.loop_indices = []
        self.loop_block_max_indices = []
        self.linear_phase_increments = []
        self.quadratic_phase_increments = []
        self.repetitions = []
        self.app_params = {}

class SPV:
    def __init__(self, waveform_params, waveform_raw_data):
        self.waveform_params = waveform_params
        self.RF_phase = waveform_raw_data[0]
        self.RF_amp = waveform_raw_data[1]
        self.Gx = waveform_raw_data[2]
        self.Gy = waveform_raw_data[3]
        self.Gz = waveform_raw_data[4]

    def get_waveform_params(self):
        return self.waveform_params
    
    def get_RF_phase(self):
        return self.RF_phase
    
    def get_RF_amp(self):
        return self.RF_amp
    
    def get_Gx(self):
        return self.Gx

    def get_Gy(self):
        return self.Gy
    
    def get_Gz(self):
        return self.Gz

def read_apd(apd_file):
    # read the .apd file into a list and return it.
    with open(apd_file, 'r') as f:
        lines = f.readlines()
    
    # find all sections delimited by "[text]"
    sections = []
    section = []
    for line in lines:
        if line[0] == '[':
            sections.append(section)
            section = []
        section.append(line)
    sections.append(section)

    # remove any empty sections
    sections = [section for section in sections if section != []]

    # remove any empty lines
    sections = [[line for line in section if line != '\n'] for section in sections]

    # remove all the newlines
    sections = [[line.strip() for line in section] for section in sections]

    # remove all lines that start with a semicolon (these are comments)
    sections = [[line for line in section if line[0] != ';'] for section in sections]

    # extract the sections that are not [Application]
    sequencing_sections = [section for section in sections if section[0] != '[Application]']

    # extract the section that is [Application]
    application_section = [section for section in sections if section[0] == '[Application]'][0]

    return [sequencing_sections, application_section]

def read_spv(spv_file):
    # read the .spv file as an XML file and return it.
    import xml.etree.ElementTree as ET
    tree = ET.parse(spv_file)
    root = tree.getroot()
    return root

def extract_repetitions(waveform):
    # sometimes we don't use a "repeat" block but instead use a "steps" block.
    # limitation: only one of these can be used at a time, and there can only be one.
    # For now we prefer "repetitions", 
    # but if it is not specified, then we use "steps".
    repetitions, _ = sb.extract_sb_parameters(waveform, ["repetitions"])
    steps, _ = sb.extract_sb_parameters(waveform, ["steps"])
    if repetitions[0] == -1:
        repetitions = steps
    else:
        repetitions = repetitions[0]
    # convert the repetitions to an integer
    repetitions = int(repetitions)
    return repetitions

def load_waveform(waveform_raw_file, waveform):
    # extract the application parameters
    repetitions = extract_repetitions(waveform)
    rfRateMultiplier = sb.extract_sb_parameters(waveform, ["rfRateMultiplier"])[0][0]
    gradLimit = sb.extract_sb_parameters(waveform, ["gradLimit"])[0][0]

    # read the binary file using ieee big endian int 16
    waveform_raw_data = np.fromfile(waveform_raw_file, dtype='>i2')

    total_len = len(waveform_raw_data)
    divisor = 3 + 2*rfRateMultiplier

    nGradSamples = int(total_len / divisor)
    nRFSamples = int(nGradSamples * rfRateMultiplier)

    # extract the gradient data
    RF_phase = waveform_raw_data[0:nRFSamples]
    RF_amp = waveform_raw_data[nRFSamples:2*nRFSamples]
    Gx = waveform_raw_data[2*nRFSamples:2*nRFSamples+nGradSamples]
    Gy = waveform_raw_data[2*nRFSamples+nGradSamples:2*nRFSamples+2*nGradSamples]
    Gz = waveform_raw_data[2*nRFSamples+2*nGradSamples:2*nRFSamples+3*nGradSamples]

    RF_phase = np.reshape(RF_phase, (repetitions, -1))
    RF_amp = np.reshape(RF_amp, (repetitions, -1))
    Gx = np.pad(np.reshape(Gx, (repetitions, -1)), (0, 1))
    Gy = np.pad(np.reshape(Gy, (repetitions, -1)), (0, 1))
    Gz = np.pad(np.reshape(Gz, (repetitions, -1)), (0, 1))

    # Convert units. Crucially important for the units to be correct.
    # kHz/m
    maxgrad_sb = gradLimit * 1e-2
    Gx = Gx * (1 / 32767) * (maxgrad_sb) * 42.58 * (1e6);
    Gy = Gy * (1 / 32767) * (maxgrad_sb) * 42.58 * (1e6);
    Gz = Gz * (1 / 32767) * (maxgrad_sb) * 42.58 * (1e6);

    return [RF_phase, RF_amp, Gx, Gy, Gz]

def read_HAWK(application_path):

    apd_file = application_path + "/application.apd"

    # parse the .apd file to get the sequence structure
    sequencing_sections, application_sections = read_apd(apd_file)

    # for each section, if it is type = waveform, then find the waveform name "file = <name>"
    apd = APD()

    for section in sequencing_sections:
        # check if the section has "type = waveform" in it (spaces are ignored)
        if "type=waveform" in [line.replace(' ', '') for line in section]:
            # find the line that has "file = <name>" (spaces are ignored)
            waveform_name = [line.replace(' ', '') for line in section if "file=" in line.replace(' ', '')][0].split('=')[1]
            waveform_raw_name = waveform_name.split('.')[0] + '.wfm'

            # find the waveform file
            waveform_file = application_path + waveform_name
            waveform_raw_file = application_path + waveform_raw_name
            # read the waveform file
            waveform = read_spv(waveform_file)
            # read the .wfm file
            s = SPV(waveform, load_waveform(waveform_raw_file, waveform))
            apd.SPVS.append(s)

        else:
            apd.SPVS.append(None)
        
        # check if the section has a "loopIndex=[<item>]" in it (spaces are ignored)
        if any(["loopIndex=" in line.replace(' ', '') for line in section]):
            # find the number after "loopIndex = " in the line.
            loop_index = [line.replace(' ', '') for line in section if "loopIndex=" in line.replace(' ','')][0].split('=')[1]
            # add the loop index to the list (int)
            apd.loop_indices.append(int(loop_index))
        else:
            # if not specified, then the loop index is 0
            apd.loop_indices.append(0)
        
        # check if the section has a "type=loop" in it (spaces are ignored)
        if "type=loop" in [line.replace(' ', '') for line in section]:
            # find the maximum number of iterations in the loop block by maximumIndex="x" (spaces are ignored)
            max_index = [line.replace(' ', '') for line in section if "maximumIndex=" in line.replace(' ', '')][0].split('=')[1]
            apd.loop_block_max_indices.append(int(max_index))
        else:
            # if not specified, then the maximum index is 0
            apd.loop_block_max_indices.append(0)
        
        if any(["linearPhaseIncrement=" in line.replace(' ', '') for line in section]):
            # find the number after "linearPhaseIncrement = " in the line.
            linear_phase_increment = [line.replace(' ', '') for line in section if "linearPhaseIncrement=" in line.replace(' ','')][0].split('=')[1]
            # add the linear phase increment to the list (float)
            linear_phase_increment = float(linear_phase_increment)
            apd.linear_phase_increments.append(linear_phase_increment)
        else:
            apd.linear_phase_increments.append(0)

        if any(["quadraticPhaseIncrement=" in line.replace(' ', '') for line in section]):
            # find the number after "quadraticPhaseIncrement = " in the line.
            quadratic_phase_increment = [line.replace(' ', '') for line in section if "quadraticPhaseIncrement=" in line.replace(' ','')][0].split('=')[1]
            # add the quadratic phase increment to the list (float)
            quadratic_phase_increment = float(quadratic_phase_increment)
            apd.quadratic_phase_increments.append(quadratic_phase_increment)
        else:
            apd.quadratic_phase_increments.append(0)

    # for each waveform, read the spinbench file and extract the number of repetitions.
    apd.repetitions =  [0 for _ in range(len(apd.SPVS))]
    for idx, spv in enumerate(apd.SPVS):
        if spv is not None:
            # extract the number of repetitions
            repetition = extract_repetitions(spv.get_waveform_params())
            apd.repetitions[idx] = repetition
        elif apd.loop_block_max_indices[idx] != 0:
            apd.repetitions[idx] = apd.loop_block_max_indices[idx]
        else:
            apd.repetitions[idx] = 0 

    # Grab parameters from the .spv files.
    # NOTE: this assumes that all .spv files have the same parameters. 
    # there is NO way to reconcile which .spv file is the "correct" one.
    # we will loop through all of them and use the first one we find.

    apd.app_params = {
        "gradLimit": None,
        "slewRateLimit": None,
        "gradSamplingRate": None,
        "rfRateMultiplier": None,
        "FA": None,
        "spatialResolution": None,
        "thickness": None,
        "fov": None,
        "rf_phase_initial": None,
        "B0": None
    }

    for spv in apd.SPVS:
        if spv is not None:
            [gradLimitFile, slewRateLimitFile, gradSamplingRateFile, rfRateMultiplierFile, FAFile, spatialResolutionFile, thicknessFile, fovFile, rf_phase_initialFile, B0File], _ = sb.extract_sb_parameters(spv.get_waveform_params(), ["gradLimit", "slewRateLimit", "gradSamplingRate", "rfRateMultiplier", "tip", "spatialResolution", "thickness", "fov", "rfPhase_0_0", "fieldStrength"])
            for param in apd.app_params:
                # check if the parameter is in the file
                if apd.app_params[param] == None and locals()[param + "File"] != None:
                    # if it is, then set the parameter
                    apd.app_params[param] = locals()[param + "File"]
    
    return apd

def check_and_increment(apd, waveform_list, time_index, rep_index):
    # recursive function to check the loop index list of lists and increment the current index
    # based on the "time_index" and "rep_index" variables.

    # copy the previous waveform_list to the current one, but only at indexes above rep_index
    # this is because we are incrementing the current index, so we don't want to overwrite it.
    for i in range(rep_index, len(waveform_list[time_index])):
        waveform_list[time_index][i] = waveform_list[time_index-1][i]

    #waveform_list[time_index] = waveform_list[time_index-1].copy()
    
    # increment the current index
    waveform_list[time_index][rep_index] += 1

    if (waveform_list[time_index][rep_index]) >= apd.repetitions[rep_index]:
        waveform_list[time_index][rep_index] = 0
        if rep_index < len(apd.repetitions)-1:
            rep_index = check_and_increment(apd, waveform_list, time_index, rep_index+1)
            return rep_index
    else:
        return rep_index

def truncate_zero_padded_amplitudes(data, extra_pad=10, start_only=False):
    # truncate the zero-padded amplitudes to the first non-zero value with a buffer of "extra_pad" samples
    # return the truncated data and the index of the first non-zero value (delay)
    # data 2D array of [repetition, time]

    if start_only:
        first_nonzero = 0
    else:
        # find the first non-zero index
        first_nonzero = np.max([np.nonzero(data)[1][0]-extra_pad,0])

    # find the last non-zero index
    last_nonzero = np.min([np.nonzero(data)[1][-1]+extra_pad-1, data.shape[1]-1])
    
    # remove the zeros from the beginning and end, with 10 samples of padding
    data = data[:,first_nonzero:last_nonzero]

    # add one to the first non-zero index to account for the delay (zero indexed)
    return data, first_nonzero

def hawk_to_pulseq(application_path):
    # main function to convert a HAWK sequence to a Pulseq sequence
    # takes input applicatio path and creates and a Pulseq sequence
    # returns the Pulseq sequence object and system object.

    apd = read_HAWK(application_path)

    # ensure the gradient sampling rate is set to the correct value
    assert apd.app_params["gradSamplingRate"] == 100, "Gradient sampling rate must be 100 us!"

    rfRasterTime = 10e-6/apd.app_params["rfRateMultiplier"]
    gradRasterTime = 10e-6
    rfDeadTime = 100e-6
    dt_adc = 1e-6

    # create a new sequence object.
    system = pp.Opts(max_grad=apd.app_params["gradLimit"]*10, grad_unit='mT/m', max_slew=apd.app_params["slewRateLimit"]*10, slew_unit='mT/m/ms',
                    rf_ringdown_time=30e-6, rf_dead_time=rfDeadTime, adc_dead_time=10e-6, rf_raster_time=rfRasterTime, grad_raster_time=gradRasterTime,
                    B0=apd.app_params["B0"])

    seq = pp.Sequence(system=system)

    # max loop index
    loop_index_range = range(0, max(apd.loop_indices)+1)

    # pythonic way to get all the repetitions for sections corresponding to each loop index
    repeats = [[apd.repetitions[i] for j in range(0, len(apd.sequencing_sections)) if apd.loop_indices[j] == i] for i in range(0, len(loop_index_range))]
    from operator import mul
    from functools import reduce
    num_blocks = reduce(mul, [min([apd.repetitions[i]]) for i in range(0,len(apd.repetitions))],1)

    # add delays to gradients due to RF dead time and compute delay times.
    for spv in apd.SPVS:
        if spv is not None:
            # add zeros to beginning of spv.Gx, spv.Gy, and spv.Gz to account for RF dead time
            if spv.Gx is not None:
                spv.Gx = np.hstack((np.zeros((spv.Gx.shape[0], int(rfDeadTime/gradRasterTime))), spv.Gx))
            if spv.Gy is not None:
                spv.Gy = np.hstack((np.zeros((spv.Gy.shape[0], int(rfDeadTime/gradRasterTime))), spv.Gy))
            if spv. Gz is not None:
                spv.Gz = np.hstack((np.zeros((spv.Gz.shape[0], int(rfDeadTime/gradRasterTime))), spv.Gz))
            
            # remove the zeros from beginning and end, but store number removed for later
            spv.Gx, spv.Gx_delayTime = truncate_zero_padded_amplitudes(spv.Gx, extra_pad=EXTRA_PAD)
            spv.Gy, spv.Gy_delayTime = truncate_zero_padded_amplitudes(spv.Gy, extra_pad=EXTRA_PAD)
            spv.Gz, spv.Gz_delayTime = truncate_zero_padded_amplitudes(spv.Gz, start_only=True, extra_pad=EXTRA_PAD)
            spv.RF_amp, spv.RF_delayTime = truncate_zero_padded_amplitudes(spv.RF_amp, start_only=True, extra_pad=EXTRA_PAD)

    # the loop index list is a list of lists that contains the loop indices for each time step
    # across all sequencing sections. This is similar to the loop index in the RTHAWK tutorials
    # under sequencing.
    loop_index_list = [[0 for i in range(len(apd.repetitions))] for j in range(num_blocks)]

    for i in range(0, len(loop_index_list)-1):
        check_and_increment(apd, loop_index_list, i+1, 0)


    phase_offset = np.zeros(len(apd.SPVS))

    for i in range(0, len(loop_index_list)):
        # start a new block, based on the loop index.
        # iterate through spvs
        for idx, spv in enumerate(apd.SPVS):
            if spv is not None:

                # increment phase.
                phase_offset[idx] = phase_offset[idx] + ((apd.linear_phase_increments[idx] / 180) * np.pi)
                # wrap phase offset to -pi to pi
                phase_offset[idx] = np.mod(phase_offset[idx] + np.pi, 2*np.pi) - np.pi

                # readout time
                readout_duration, readout_position = sb.extract_sb_parameters(spv.get_waveform_params(), ["readoutDuration"])
                if spv.Gx is not None:
                    Gx = pp.make_arbitrary_grad.make_arbitrary_grad(channel='x', waveform=spv.Gx[loop_index_list[i][idx]], system=system, delay=apd.SPVS[idx].Gx_delayTime*gradRasterTime)
                if spv.Gy is not None:
                    Gy = pp.make_arbitrary_grad.make_arbitrary_grad(channel='y', waveform=spv.Gy[loop_index_list[i][idx]], system=system, delay=apd.SPVS[idx].Gy_delayTime*gradRasterTime)
                if spv.Gz is not None:
                    Gz = pp.make_arbitrary_grad.make_arbitrary_grad(channel='z', waveform=spv.Gz[loop_index_list[i][idx]], system=system, delay=apd.SPVS[idx].Gz_delayTime*gradRasterTime)
                if spv.RF_amp is not None:
                    tip, tip_pos = sb.extract_sb_parameters(spv.waveform_params, ["tip"])
                    tip = tip[0]
                    RF = pp.make_arbitrary_rf(signal=spv.RF_amp[loop_index_list[i][idx]], flip_angle=tip * np.pi / 180, system=system, delay=apd.SPVS[idx].RF_delayTime*rfRasterTime)
                    RF.phase_offset = phase_offset[idx]
                if readout_duration[0] > 0:
                    ADC = pp.make_adc(num_samples=int(readout_duration[0] * 1e-3 / dt_adc), system=system, delay=(round(readout_position[0]*1e-3 / gradRasterTime) + EXTRA_PAD) * gradRasterTime, dwell=dt_adc)
                    ADC.phase_offset = phase_offset[idx]
                    # ADC = pp.make_adc(num_samples=spv.ADC.numSamples, system=system, delay=spv.ADC.delayTime*dt_adc)
                # add the block to the sequence
                seq.add_block(Gx, Gy, Gz, RF, ADC)

    print("--------- Timing Check --------")
    if seq.check_timing():
        print("-------------------------------")
        print("Timing check passed!")
        print("-------------------------------")
    else:
        print("-------------------------------")
        print("Timing check failed!")
        print(seq.timing_report())
        print("-------------------------------")

    # return seq, system, number of blocks, apd-structure.
    return seq, system, len(loop_index_list), apd