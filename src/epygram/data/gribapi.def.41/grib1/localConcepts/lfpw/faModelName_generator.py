#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re

from jinja2 import Template

re_comment = re.compile(r'^\s*([#%!]|$)')

if __name__ == "__main__":

    defs = dict(
        n_asscom_slots=20,
        n_members_pearp=100,
        n_members_pearo=100,
    )

    # tuning
    align_equal_signs = True

    # read the template file, removing comments and empty lines asap
    with open('faModelName.source') as fp:
        contents = ''.join([
            s.strip(' \t')
            for s in fp.readlines()
            if not re_comment.match(s)
        ])

    # define the template object, asking for no empty lines generation
    t = Template(
        contents,
        trim_blocks=True,
        lstrip_blocks=True,
    )

    # apply the template with our variables
    # defs is here a dictionary (easier to find in the template)
    # we could use **defs or any variable definition
    contents = t.render(defs=defs).split('\n')

    # sort lines for easy comparison between versions
    contents = sorted(contents)

    if align_equal_signs:
        pairs = [line.split('=', 1) for line in contents]
        maxlen = max([len(lhs) for lhs in list(zip(*pairs))[0]])
        fmt = '{' + ':{}s'.format(maxlen) + '}={}'
        contents = [fmt.format(u, v) for (u, v) in pairs]

    with open('faModelName.def.with_jinja', 'w') as fp:
        fp.write('\n'.join(contents) + '\n')
