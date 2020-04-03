#!/usr/bin/python3

import sys

#filename = sys.argv[1] + 'sorted_freq.txt'
basename_file = sys.argv[1]

basename_list = []
with open(basename_file, 'r') as basenames:
    for basename in basenames:
        filename = basename.strip() + '_sorted_freq.txt'
        counter = 1
        total = 0
        smallest_val = None
        largest_val = None
        with open(filename, 'r') as current_analysis_file:
            for line in current_analysis_file:
                if counter == 1:
                    smallest_val = line.strip()
                elif counter == 93:
                    largest_val = line.strip()
                total += int(line)
                counter += 1
            print(basename.strip() + ':')
            print('SMALLEST VALUE:', smallest_val)
            print('AVERAGE:', total/93)
            print('LARGEST VALUE', largest_val)
            print('\n')
