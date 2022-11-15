import importlib, os, sys
"""
Install requirements for Agdas.
"""

# path of agdas project
script_path = os.path.dirname(os.path.realpath(__file__))

ins = False
# open requirements
with open('requirements.txt', 'r') as f:

    for l in f.read().splitlines():
        library = l.split('=')[0]

        try:
            importlib.import_module(library)

        except ModuleNotFoundError:
            ins = True

if ins:
    os.system('{} -m pip install pipreqs'.format(sys.executable))
    os.system('{} -m pipreqs.pipreqs {} --encoding utf-8-sig --force'.format(sys.executable, script_path))
    os.system('{} -m pip install -r requirements.txt'.format(sys.executable))