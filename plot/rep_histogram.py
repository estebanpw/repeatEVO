from Bio.Blast import NCBIWWW
import os
import matplotlib.pyplot as plt
from textwrap import wrap
import sys
import getopt
import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main(argv):
    input_path = ''
    # path of execution
    output_path = os.path.dirname(os.path.abspath(__file__))
    file = ''

    # Get all arguments given
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["inputfile=", "outputfolder="])
    except getopt.GetoptError:
        print('rep_histogram.py -index <inputfile> -o <outputfolder>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('rep_histogram.py -index <inputfile> -o <outputfolder>')
            sys.exit()
        elif opt in ("-i", "--inputfile"):
            if os.path.isabs(arg):
                file = os.path.basename(arg)
                input_path = os.path.dirname(arg)
            else:
                path = os.path.join(os.getcwd(), arg)
                file = os.path.basename(path)
                input_path = os.path.dirname(path)
                print('path: ' + path + ' file: ' + file + ' inputp: ' + input_path)
        elif opt in ("-o", "--outputfolder"):
            if os.path.isabs(arg):
                output_path = arg
            else:
                output_path = os.path.join(os.getcwd(), arg)
                # print('arg: ' + arg + ' outputp: ' + output_path)

    #print(f'Blast file: {file}')
    df = pd.read_csv(os.path.join(input_path, file),
                     sep='\t',
                     names=['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'START', 'END', 'C7', 'C8', 'C9',
                            'C10'])

    max_length = 0
    for index, row in df.iterrows():
        #print(row['START'], row['END'])
        if row['END'] > max_length:
            max_length = row['END']
    print(max_length)
    canvas = np.zeros(max_length)

    for index, row in df.iterrows():
        for hit in range(row['START'], row['END']):
            canvas[hit] += 1
    print(canvas)

    plt.figure(figsize=(10, 10))
    plt.plot(np.arange(max_length), canvas)
    plt.title("\n".join(wrap('Blast significant profile of: ' + file, 60)))
    plt.xlabel('Repetitions length')
    plt.ylabel('Repetitions')
    #plt.show()
    plt.savefig(os.path.join(output_path, file + '_histogram.png'))

    '''
        Each plot in hist function is a subplot so we have to load them as a subplots() 
        
        fig, ax = plt.subplots()
        df.hist(column='START', ax=ax, alpha=0.5)
        plt.show()
        fig.savefig(os.path.join(output_path, file + '_histogram.png'))
        #ax.plot()
        #ax.show()
        # print(df[[6]])
    '''


matplotlib.use('Agg')

if __name__ == "__main__":
    print('RepHistogram v0.1')
    main(sys.argv[1:])
