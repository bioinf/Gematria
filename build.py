#!/usr/bin/env python3
export = 'gematria.standalone.py'


def importer(filename, head=True):
    print ('\033[31m > {0}\033[0m'.format(filename))
    if head:
        data = '# ' + '-' * 75 + ' #\n'
        data += '# filename: ' + filename + '\n\n\n'
    else:
        data = ''

    with open(filename, 'r') as f:
        for line in f:
            if line[0:len("from include")] == "from include":
                name = line.split(' ')[1].replace('.', '/') + '.py'
                data += importer(name)
            else:
                data += line
    return data

code = open(export, 'w')
code.write(importer('gematria.py', False))

print ('\033[32m > Done: {0}\033[0m'.format(export))
