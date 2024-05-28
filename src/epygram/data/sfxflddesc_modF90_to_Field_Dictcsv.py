#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) Météo France (2014-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info
from __future__ import print_function, absolute_import, unicode_literals, division

import sys
import io
from bronx.datagrip.misc import read_dict_in_CSV
import argparse

type2nature = {'X':'float', 'L':'bool', 'C':'str', 'N':'int', 'Y':'float',
               'T':'?T?'}


def main(mod, csv, out=sys.stdout):
    """
    Read sfxflddesc_mod.F90 and Field_Dict_FA.csv, then compare them.
    """
    # 1. parsing sfxflddesc_mod.F90
    with io.open(mod, 'r') as d:
        lines = d.readlines()
    fields = {}
    data_indexes = []
    for i, l in enumerate(lines):
        if l.startswith('DATA') and 'CDESC' in l:
            data_indexes.append(i)
    data_indexes.append(-1)
    for i in range(len(data_indexes[:-1])):
        for l in range(data_indexes[i] + 1, data_indexes[i + 1]):
            line = lines[l].replace('"', '').strip()
            if line.startswith("1234567890"):  # we skip the first line of CDESC
                continue
            if line.endswith(",&"):
                line = line.split(".")
                line = [e for e in line if e != '']
                if len(line) < 3:
                    continue
                else:
                    name = line[0]
                    vartype = line[2]
                    if vartype[0] != 'X':  # H2D fields
                        fields[name] = {'type':'Misc',
                                        'nature':type2nature[vartype[0]],
                                        'dimension':vartype[1]}
            elif lines[l].endswith('/'):
                break
    mod = {'SFX.' + f:{'nature':fields[f]['nature'],
                       'dimension':fields[f]['dimension']}
           for f in fields.keys() if fields[f]['type'] == 'Misc'}
    # 2. parsing csv
    csv = read_dict_in_CSV(csv)[0]
    csv = {f['name']:{'nature':f['nature'],
                      'dimension':str(f['dimension'])}
           for f in csv if f['type'] == 'Misc'}
    # 3. compare
    inconsistencies = {}
    missing_in_csv = []
    missing_in_mod = []
    for f in mod.keys():
        if f in csv.keys():
            if mod[f] != csv[f]:
                inconsistencies[f] = {'mod':mod[f], 'csv':csv[f]}
        else:
            missing_in_csv.append(f)
    for f in csv.keys():
        if f not in mod.keys():
            missing_in_mod.append(f)
    # 4. print out
    out.write('-' * 80 + '\n')
    out.write('inconsistencies: ' + str(len(inconsistencies)) + '\n')
    out.write('-' * 80 + '\n')
    for f, d in sorted(inconsistencies.items()):
        out.write('{}: mod => {}, csv => {}\n'.format(f,
                                                      str(d['mod']),
                                                      str(d['csv'])))
    out.write('-' * 80 + '\n')
    out.write('missing_in_csv: ' + str(len(missing_in_csv)) + '\n')
    for f in sorted(missing_in_csv):
        out.write('name:{},type:Misc,nature:{nature},dimension:{dimension}\n'.
                  format(f, **mod[f]))
    out.write('-' * 80 + '\n')
    out.write('missing_in_mod: ' + str(len(missing_in_mod)) + '\n')
    out.write('-' * 80 + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read sfxflddesc_mod.F90 and Field_Dict_FA.csv, then compare them.")
    parser.add_argument('--csv', required=True, help='the Field_Dict_FA.csv')
    parser.add_argument('--mod', required=True, help='the sfxflddesc_mod.F90')
    args = parser.parse_args()
    main(csv=args.csv, mod=args.mod)
