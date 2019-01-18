import os
import numpy as np
import pickle
import datetime
import shutil


MASSBANK_LOCATION = 'MassBank-data-master'



def parse_record_file(filepath):
    data = [''] * 3
    with open(filepath, 'r') as f:
        line = f.readline()
        while line:
            if line[:18] == 'AC$INSTRUMENT_TYPE':
                data[0] = line.split()[1]
            elif line[:20] == 'AC$MASS_SPECTROMETRY':
                elements = line.split()
                if elements[1] == 'MS_TYPE':
                    data[1] = elements[2]
                elif elements[1] == 'ION_MODE':
                    data[2] = elements[2]
                    break
            line = f.readline()
    return data


def get_info(filepath):
    intensities = list()
    with open(filepath, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if len(line) > 1 and not line[0].isalpha():
                elements = line.split()
                if len(elements) > 1:
                    try:
                        intensities.append(float(elements[1]))
                    except ValueError:
                        pass

    relint = np.median(intensities) / np.max(intensities)
    return relint


def load_records():
    records_filepath = 'records.pkl'
    if not os.path.isfile(records_filepath):
        records = dict()
        for subdir in os.listdir(MASSBANK_LOCATION):
            folder = os.path.join(MASSBANK_LOCATION, subdir)
            if os.path.isdir(folder):
                for filename in os.listdir(folder):
                    filepath = os.path.join(folder, filename)
                    if filepath.lower().endswith('.txt'):
                        data = parse_record_file(filepath)
                        records[filepath] = data
        with open(records_filepath, 'wb') as f:
            pickle.dump(records, f)
    else:
        with open(records_filepath, 'rb') as f:
            records = pickle.load(f)
    return records


records = load_records()

shutil.rmtree('dataset1'); os.mkdir('dataset1')
i = 1
for filepath in records.keys():
    instrument, ms_type, ion_mode = records[filepath]
    if ms_type == 'MS' and ion_mode == 'POSITIVE' and 'ESI-QTOF' in instrument:
        shutil.copyfile(filepath, os.path.join('dataset1', '%i.txt' % i))
        print(i)
        i += 1

shutil.rmtree('dataset2'); os.mkdir('dataset2')
i = 1
for filepath in records.keys():
    instrument, ms_type, ion_mode = records[filepath]
    if ms_type == 'MS2' and ion_mode == 'POSITIVE' and 'ESI-QTOF' in instrument:
        relint = get_info(filepath)
        if abs(relint - 0.1) < 0.02:
            shutil.copyfile(filepath, os.path.join('dataset2', '%i.txt' % i))
            i += 1
