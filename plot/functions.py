import numpy as np


def get_color(nucleotide):
    if nucleotide == 'A':
        # red
        return [1.0, 0.0, 0.0]

    elif nucleotide == 'C':
        # green
        return [0.0, 1.0, 0.0]

    elif nucleotide == 'G':
        # blue
        return [0.0, 0.0, 1.0]

    elif nucleotide == 'T':
        # yellow
        return [1.0, 1.0, 0.0]
	elif nucleotide == 'N':
		# white
		return [1.0, 1.0, 1.0]

    else:
        # black
        return [0.0, 0.0, 0.0]


def get_most_common(data, coordinate, repetitions):
    nucleotide_count = [0, 0, 0, 0]
    # print(repetitions)
    repetitions = repetitions - 1
    for rep in range(0, repetitions):
        # print(data[coordinate, rep])
        #print(rep)
        if data[rep, coordinate] == 'A':
            nucleotide_count[0] += 1
        elif data[rep, coordinate] == 'C':
            nucleotide_count[1] += 1
        elif data[rep, coordinate] == 'G':
            nucleotide_count[2] += 1
        elif data[rep, coordinate] == 'T':
            nucleotide_count[3] += 1

    if nucleotide_count[0] > nucleotide_count[1] and nucleotide_count[0] > nucleotide_count[2] and \
            nucleotide_count[0] > nucleotide_count[3]:
        return 'A', nucleotide_count[0]
    elif nucleotide_count[1] > nucleotide_count[0] and nucleotide_count[1] > nucleotide_count[2] and \
            nucleotide_count[1] > nucleotide_count[3]:
        return 'C', nucleotide_count[1]
    elif nucleotide_count[2] > nucleotide_count[1] and nucleotide_count[2] > nucleotide_count[0] and \
            nucleotide_count[2] > nucleotide_count[3]:
        return 'G', nucleotide_count[2]
    elif nucleotide_count[3] > nucleotide_count[1] and nucleotide_count[3] > nucleotide_count[2] and \
            nucleotide_count[3] > nucleotide_count[0]:
        return 'T', nucleotide_count[3]
    else:
        return '-', 0
