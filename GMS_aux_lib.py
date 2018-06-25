# Auxiliary functions for GMS builder
import argparse
import os
import sys
import time
import numpy as np
import multiprocessing
from functools import reduce

def track_make(data):

    chrom,reads_stat = data   
    
    print ('Making GMS-track. Chromosome ' + chrom + '.') 
    
    seq_lng = len( reads_stat )
    
    lng = min(5000000, seq_lng)
    start_arr = np.zeros((lng))
    span_arr = np.zeros((lng))
    current_value_arr = np.zeros((lng))
    chrom_arr = ['' for _ in range(lng)]     
    
    start = 0
    span = 0
    current_value = reads_stat[0]
    series_counter = 0
        
    for i in range(seq_lng):
        span += 1
           
        if (i == seq_lng-1):
            series_counter += 1
            start_arr[series_counter-1] = start+1
            span_arr[series_counter-1] = span
            current_value_arr[series_counter-1] = current_value
            chrom_arr[series_counter-1] = chrom
            continue
                
        if not( reads_stat[i] == reads_stat[i + 1] ):
            
            series_counter += 1
            if series_counter == lng+1:
                start_arr = np.append(start_arr, np.zeros((lng)) )
                span_arr = np.append(span_arr, np.zeros((lng)))  
                current_value_arr = np.append(current_value_arr, np.zeros((lng)))
                chrom_arr+=['' for _ in range(lng)]
                lng += lng 
            start_arr[series_counter-1] = start+1
            span_arr[series_counter-1] = span
            current_value_arr[series_counter-1] = current_value
            chrom_arr[series_counter-1] = chrom
            start = i + 1
            span = 0
            current_value = reads_stat[i+1]
            
    start_arr = start_arr[:series_counter] 
    span_arr = span_arr[:series_counter] 
    current_value_arr = current_value_arr[:series_counter]
    chrom_arr= chrom_arr[:series_counter]

    current_chr = zip(start_arr, span_arr, current_value_arr, chrom_arr)
        
    print (str(series_counter) + ' series are found.')
    return(current_chr)
    
# Write data to Wig, BigWig, Bed, BigBed, Tdf file
def track_write(chromosome_len, chromosome_names, chromosome_GMS, in_file, project_name, output_formats, threads, current_dir):
            
    if output_formats['Wig'] or output_formats['TDF']:
        handle = open('./{0}/{0}.wig'.format(project_name), 'w')
    if output_formats['Bed'] or output_formats['bigBed']:
        bed_handle = open('./{0}/{0}.bed'.format(project_name), 'w')
    if output_formats['bigWig']:
        import pyBigWig
        bw = pyBigWig.open('./{0}/{0}.bw'.format(project_name), 'w')

    # Sorting
    chromosome_len, chromosome_names, chromosome_GMS = zip(*sorted(zip(chromosome_len, chromosome_names, chromosome_GMS), reverse = True))
    print ('Track-making was started at {}.'.format(time.ctime(int(time.time()))) )

    if threads > 1:
        pool = multiprocessing.Pool(processes= threads)
        result_chr = pool.map(track_make, list( zip(chromosome_names, chromosome_GMS)) )
        pool.close()
    else:
        result_chr = map(track_make, list( zip(chromosome_names, chromosome_GMS)) )        
    
    if output_formats['bigWig']:
         bw.addHeader( list(zip( chromosome_names, chromosome_len )) ) 
    
    for current_chr in result_chr:
        for entry in current_chr:
            if output_formats['Wig'] or output_formats['TDF']:
                handle.write("fixedStep chrom={0} start={1} step=1 span={2}\n".format(entry[3], int(entry[0]), int(entry[1])))
                handle.write((str(entry[2])) + "\n")
           
            if output_formats['bigWig']:
                 bw.addEntries(entry[3], [int((entry[0]-1))] ,values = [entry[2]], span = int(entry[1]), step = 1)
            
            if output_formats['Bed'] or output_formats['bigBed']:
                color = '{0},{0},255'.format(str(int(255-np.round(255*entry[2]/100))))
                    
                bed_handle.write("{0}\t{1}\t{2}\t.\t{3}\t.\t{1}\t{2}\t{4}\n".format(entry[3], int(entry[0]-1), int(entry[0]-1 + entry[1]), int(entry[2]), color))
    if output_formats['Wig'] or output_formats['TDF']:
        handle.close()
    if output_formats['bigWig']:
        bw.close()
    if output_formats['Bed'] or output_formats['bigBed']:
        bed_handle.close()

    print ('Track-making was finished at {}.'.format(time.ctime(int(time.time()))) )

    if output_formats['bigBed']:   
        print('Converting BED to BigBED.')
        
        bs_cmd = "sort -k1,1 -k2,2n ./{0}/{0}.bed > sorted.bed".format(project_name)
        os.system(bs_cmd)

        bs_cmd1 = "{0}/utils/faSize {1} -detailed > ./Cache/{2}.sizes".format(current_dir, in_file, project_name)
        os.system(bs_cmd1)
       
        bs_cmd2 = "{0}/utils/bedToBigBed -type=bed9 sorted.bed ./Cache/{1}.sizes ./{1}/{1}.bb".format(current_dir,project_name)
        os.system(bs_cmd2)
        
        if output_formats['Bed'] == False:
            bs_cmd3 = "rm ./{0}/{0}.bed".format(project_name)
            os.system(bs_cmd3)

    if output_formats['TDF']:
        print('Converting Wig to TDF')
        bs_cmd = "{0}/utils/faSize {1} -detailed > ./Cache/{2}.chrom.sizes".format(current_dir, in_file, project_name)
        os.system(bs_cmd)
        bs_cmd1 = "{0}/utils/igvtools toTDF ./{1}/{1}.wig ./{1}/{1}.tdf ./Cache/{1}.chrom.sizes".format(current_dir, project_name)
        os.system(bs_cmd1)
        if output_formats['Wig'] == False:
            bs_cmd3 = "rm ./{0}/{0}.wig".format(project_name)
            os.system(bs_cmd3)  
    
    return 0
        
# Parse input arguments
def parse_inputs(current_dir):
    parser = argparse.ArgumentParser(description='Parser for GMS-track builder')
    parser.add_argument('-i', '--input' , help='Input file' , metavar='Str',
                        type=str, required=True)
    parser.add_argument('-o', '--output'  , help='Output file', metavar='Str',
                        type=str, required=True)
    parser.add_argument('-l', '--read_lng', help='Read length (default 100)', metavar='Int',
                        type=int, default=100)
    parser.add_argument('-r', '--pe', help='Reads type parameters: N:mu:sigma (for Normal distribution of insertions size) or U:a:b (for Uniform distribution of insertion size) or S (for single-end reads).',
                        type=str, default = 'S')
    parser.add_argument('-f', '--formats', help='List of required output formats in a way format1,format2,... etc. Supported formats are Wig, bigWig, Bed, bigBed, TDF.',
                        type=str, default = 'bigWig')
    parser.add_argument('-m', '--mat_file', help='Save debug data in Matlab MAT format', action='store_true')  
    parser.add_argument('-t', '--threads', help = 'Number of threads', metavar = 'Int', type = int, default = 1)   

    args = parser.parse_args()
    
    print ('Checking requested formats.')
    requirements = {'Bed':{},'bigBed':{'faSize','bedToBigBed'},'Wig':{},'bigWig':{},'TDF':{'faSize','igvtools','igvtools.jar'}}
    ouptput_formats = {'Bed':False,'bigBed':False,'Wig':False,'bigWig':False,'TDF':False}
    
    requested_formats = args.formats.split(',')
    

    if not( os.path.isdir(current_dir+'/utils') ):
        existing_utils = []
    else:
        existing_utils = set(os.listdir(current_dir+'/utils'))
        
    for frmt in requested_formats:
        if frmt == 'bigWig':
            try:
                import pyBigWig
                ouptput_formats[frmt] = True
                print('pyBigWig library for bigWig was found. Will try to generate.')
            except:
                print('Cannot load pyBigWig library. BigWig output will be omitted. You still can generate Wig file and then convert it with UCSC wigToBigWig utility manually.')
            continue             

        if frmt in requirements:
            if requirements[frmt].issubset(existing_utils):
                ouptput_formats[frmt] = True
                print('Utilities for {} are found. Will try to generate.'.format(frmt))
            else:
                print('Utilities for {} weren\'t found. Will be ignored.'.format(frmt))                
        else:
            print ('Format {} is unsupported. Will be ignored.'.format(frmt))
    
    if not reduce((lambda x,y: x or y),ouptput_formats.values()):
        print ('There is no output formats available.')
        sys.exit()
        
    print ('Checking reads paramters.')
    pe_correctness_flag = False
    insertion_size = {'reads_type':'','param1':0,'param2':0}
    if args.pe == 'S':
        insertion_size['reads_type'] = 'S'
        pe_correctness_flag = True
    else:
        tmp = args.pe.split(':')
        if (len(tmp) == 3) and (tmp[0] in {'U','N'}):
            try:
                x = int(tmp[1])
                x = int(tmp[2])                
            except:
                pass
            else:
                insertion_size['reads_type'] = tmp[0]            
                insertion_size['param1'] = int(tmp[1])
                insertion_size['param2'] = int(tmp[2])
                pe_correctness_flag = True
        
    if not pe_correctness_flag:
        print('Reads type is wrong (see help for explanations).')
        sys.exit()
        
    return args.input, args.output, args.read_lng, insertion_size, ouptput_formats, args.mat_file, args.threads