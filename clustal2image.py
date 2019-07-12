import os
from os import listdir
from os.path import isfile, join
from Bio import AlignIO
from functions import get_color, get_most_common
import numpy as np
import matplotlib.pyplot as plt
from textwrap import wrap
import sys
import getopt


def main(argv):
    input_path = ''
    output_path = ''
    file = ''

    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["inputfile=", "outputpath="])
    except getopt.GetoptError:
        print('main.py -i <inputfile> -o <outputpath>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('main.py -i <inputfile> -o <outputpath>')
            sys.exit()
        elif opt in ("-i", "--inputfile"):
            file = arg
        elif opt in ("-o", "--outputpath"):
            output_path = os.path.dirname(os.path.abspath(__file__))

    # I/O Config
    #input_path = os.getcwd() + '\\input\\'
    #output_path = os.getcwd() + '\\output\\'

    #onlyfiles = [f for f in listdir(input_path) if isfile(join(input_path, f))]

    #for file in onlyfiles:
    try:
        os.mkdir(output_path + file)
    except FileExistsError:
        print("Directory ", file, " already exists")

    print(file)
    alignment = AlignIO.read(open(input_path + file), 'clustal')

    # Limits
    repetitions = 0
    coordinates = 0
    for record in alignment:
        repetitions += 1
        coordinates = len(record)

    canvas = np.zeros((repetitions, coordinates, 3))
    raw_data = np.zeros((repetitions, coordinates, 1), dtype=object)

    x_rep = 0
    y_coo = 0
    repetitions_id = []

    # Generate plot array from alignment file
    for record in alignment:
        repetitions_id.append(record.id)
        for nucleotide in record.seq:
            canvas[x_rep, y_coo] = get_color(nucleotide)
            raw_data[x_rep, y_coo] = nucleotide
            y_coo += 1
        x_rep += 1
        y_coo = 0

    # Calculate conservation
    conservation = np.zeros(coordinates, dtype=object)
    for coo in range(0, coordinates):
        conservation[coo] = get_most_common(raw_data, coo, repetitions)

    conservation_amount = [(item[1] / repetitions) * 100 for item in conservation]

    # Draw conservation
    plt.title("\n".join(wrap('Nucleotide conservation of: ' + file, 60)))
    plt.xlabel('Coordinates')
    plt.ylabel('Percentage')
    # plt.plot(np.arange(coordinates), conservation_amount)

    #fig, ax = plt.subplots()
    # ax.plot(x, y, label="Original")
    # Compute moving averages using different window sizes
    #window_lst = [int(coordinates * 0.02), int(coordinates * 0.03), int(coordinates * 0.04), int(coordinates * 0.05)]
    #y_avg = conservation_amount
    plt.ylim(0, 100)
    windows_size = int(coordinates * 0.02)
    avg_mask = np.ones(windows_size) / windows_size
    y_avg = np.convolve(conservation_amount, avg_mask, 'same')
    plt.plot(np.arange(coordinates), y_avg)
    # Save
    plt.savefig(output_path + '_conservation.png')
    plt.show()

    # Draw plot
    plt.imshow(canvas, aspect='auto')
    # Decorate
    plt.title("\n".join(wrap('Multiple alignment of repetitions and their borders of: ' + file, 60)))
    plt.xlabel('Coordinates')
    plt.ylabel('Repetitions')
    plt.yticks(np.arange(repetitions), repetitions_id)
    # Save
    plt.savefig(output_path + '_alignment.png')
    plt.show()

