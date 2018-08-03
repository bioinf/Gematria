# Build GMS-track

import sys
import time
import numpy as np
import scipy.io
import makegms
import os
from math import sqrt, pi, exp
from Bio import SeqIO
from collections import defaultdict
from GMS_aux_lib import track_write, parse_inputs

if __name__ == '__main__':    

    start_time = time.time()          
    current_dir = os.path.dirname(os.path.realpath(__file__))
    parsed_settings = parse_inputs(current_dir) # Parse command-line arguments        
    
    if parsed_settings == -1:        
        sys.exit()
    else:
        in_file, project_name, read_lng, insertion_size, output_formats,  mat_flag, threads = parsed_settings

    if not( os.path.isdir("./{}".format('Cache')) ):
        os.system('mkdir {}'.format('Cache'))
   
    if not( os.path.isdir("./{}".format(project_name)) ):
        os.system('mkdir {}'.format(project_name))
    
    makegms.run(in_file, read_lng, os.getcwd(), threads)
    #run_cmd = "{0}/makeRawGMS {1} {2} {3} {4}".format(current_dir, in_file, read_lng, threads, os.getcwd())
    #os.system(run_cmd)
        
    map_bin_track = np.unpackbits( np.fromfile('./Cache/GMS_track.bin', dtype = "uint8") )
    
    print ('Splitting GMS-track into chromosoms and windowing was started at {}.'.format(time.ctime(int(time.time()))) )

    length_array = []
     
    chromosome_len = []
    chromosome_names = []
    chromosome_GMS = []

    if (insertion_size['reads_type'] == "U"):
        min_dist = insertion_size['param1']
        sums = np.ones(insertion_size['param2'] - insertion_size['param1']) / (insertion_size['param2'] - insertion_size['param1'])
    elif (insertion_size['reads_type'] == "N"):
        mu, sigma = insertion_size['param1'], insertion_size['param2']
        min_dist = mu - 3 * sigma
        sums = [exp(-(float(i) - mu)**2 / (2 * sigma**2)) / (sqrt(2*pi)*sigma) for i in range(mu - 3 * sigma, mu + 3 * sigma + 1)]

    sliding_window = 100*np.ones(read_lng)/read_lng    
    index = 0
    with open(in_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_lng = len(record.seq) # Length of the sequence
                        
            print( 'Processing {0}. Sequence length is {1} bps.'.format(record.name,seq_lng) )    
            
            chromosome_len.append(seq_lng)
            chromosome_names.append( record.name )
            num_of_reads = seq_lng-read_lng+1 # Number of reads in the chromosome
            subseq = map_bin_track[index:index+num_of_reads]
            index = index+seq_lng
            
            if insertion_size['reads_type'] == "S":
                final_subseq = subseq
            else:
                num_of_reads = len(subseq)
                unique_pairs = np.convolve(subseq, sums)
                unique_l = np.concatenate((np.zeros(min_dist), unique_pairs))[:num_of_reads]
                unique_r = np.concatenate((unique_pairs, np.zeros(min_dist)))[-num_of_reads:]                       
                uniques = (unique_l + unique_r) / 2            
                final_subseq = subseq + (np.ones(num_of_reads) - subseq) * uniques                
                
            GMS = np.round( np.convolve(final_subseq, sliding_window) )
            chromosome_GMS.append(GMS)
            
            if mat_flag:
                scipy.io.savemat('./{0}/{0}_{1}.mat'.format(project_name, record.name), mdict={'GMS': GMS})
    
    print ('All GMS-tracks were calculated at {}.'.format(time.ctime(int(time.time()))))
        
    flag = track_write(chromosome_len, chromosome_names, chromosome_GMS, in_file, project_name, output_formats, threads, current_dir)
    
    os.system('rm -r ./Cache/') # Clear cache

    stop_time = time.time()
    
    print ('GeMaTrIA has finished at {}.'.format(time.ctime(int(time.time()))) )
    print ('Elapsed time: {}'.format(stop_time-start_time) ) 
    print ('See {} directory.'.format(project_name) )