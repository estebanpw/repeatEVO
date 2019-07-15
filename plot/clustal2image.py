import os
from Bio import AlignIO
from functions import get_color, get_most_common
import numpy as np
import matplotlib.pyplot as plt
from textwrap import wrap
import sys
import getopt
import matplotlib


def main(argv):
    input_path = ''
    # path of execution
    output_path = os.path.dirname(os.path.abspath(__file__))
    file = ''

    # Get all arguments given
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["inputfile=", "outputfolder="])
    except getopt.GetoptError:
        print('clustal2image.py -i <inputfile> -o <outputfolder>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('clustal2image.py -i <inputfile> -o <outputfolder>')
            sys.exit()
        elif opt in ("-i", "--inputfile"):
            if os.path.isabs(arg):
                file = os.path.basename(arg)
                input_path = os.path.dirname(arg)
            else:
                path = os.path.join(os.getcwd(), arg)
                file = os.path.basename(path)
                input_path = os.path.dirname(path)
                # print('path: ' + path + ' file: ' + file + ' inputp: ' + input_path)
        elif opt in ("-o", "--outputfolder"):
            if os.path.isabs(arg):
                output_path = arg
            else:
                output_path = os.path.join(os.getcwd(), arg)
                # print('arg: ' + arg + ' outputp: ' + output_path)

    # Batch config
    # input_path = os.getcwd() + slash + 'input' + slash
    # output_path = os.getcwd() + slash + 'output' + slash
    # onlyfiles = [f for f in listdir(input_path) if isfile(join(input_path, f))]
    # for file in onlyfiles:

    print('Clustal file: ' + file)
    alignment = AlignIO.read(open(os.path.join(input_path, file)), 'clustal')

    # Declare variables
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
    fig_size = (10, 10)

    plt.figure(figsize=fig_size)
    plt.title("\n".join(wrap('Nucleotide conservation of: ' + file, 60)))
    plt.xlabel('Coordinates')
    plt.ylabel('Percentage')

    # Compute moving averages using different window sizes
    plt.ylim(0, 100)
    windows_size = int(coordinates * 0.02)
    avg_mask = np.ones(windows_size) / windows_size
    y_avg = np.convolve(conservation_amount, avg_mask, 'same')

    plt.plot(np.arange(coordinates), y_avg)
    # Save
    plt.savefig(os.path.join(output_path, file + '_conservation.png'))
    print('Output: ' + os.path.join(output_path, file + '_conservation.png'))
    plt.close()

    # Draw plot
    plt.figure(figsize=fig_size)
    plt.imshow(canvas, aspect='auto')

    plt.title("\n".join(wrap('Multiple alignment of repetitions and their borders of: ' + file, 60)))
    plt.xlabel('Coordinates')
    plt.ylabel('Repetitions')
    plt.yticks(np.arange(repetitions), repetitions_id)
    # Save
    plt.savefig(os.path.join(output_path, file + '_alignment.png'))
    print('Output: ' + os.path.join(output_path, file + '_alignment.png'))
    plt.close()


matplotlib.use('Agg')

if __name__ == "__main__":
    print('Clustal2Image v0.8')
    main(sys.argv[1:])
