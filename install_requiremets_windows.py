import importlib, os, sys
"""
Install requirements for Agdas.
"""

upgrade = True

# path of agdas project
script_path = os.path.dirname(os.path.realpath(__file__))

# upgrade pip
print('=======================================')
print('Upgrade pip')
cmd = '{} -m pip install pip --upgrade'.format(sys.executable)
print(cmd)

# install pipreqs
try:
    importlib.import_module('pipreqs')
    # create requirements.txt by pipreqs
    os.system('{} -m pipreqs.pipreqs {} --encoding utf-8-sig --force'.format(sys.executable, script_path))
except ModuleNotFoundError:
    print('=======================================')
    mess = '{} is not installed'.format('pipreqs')
    print(mess)
    cmd = '{} -m pip install {}'.format(sys.executable, 'pipreqs')
    print(cmd)
    os.system(cmd)
    print('=======================================')

# open requirements
with open('requirements.txt', 'r') as f:

    for l in f.read().splitlines():
        library = l.split('=')[0]

        # try load the library
        try:
            im = importlib.import_module(library)
            mess = '{} is installed'.format(library)
            print(mess)

        # install the library if it doesn't
        except ModuleNotFoundError:
            print('=======================================')
            mess = '{} is not installed'.format(library)
            print(mess)
            cmd = '{} -m pip install {}'.format(sys.executable, library)
            print(cmd)
            os.system(cmd)
            print('=======================================')

# upgrade libraries
if upgrade:
    print('=======================================')
    print('Upgrading of libraries')
    with open('requirements.txt', 'r') as f:

        for l in f.read().splitlines():
            library = l.split('=')[0]
            print('=======================================')
            cmd = '{} -m pip install {} --upgrade'.format(sys.executable, library)
            print(cmd)