from pathlib import Path

print(Path('.').absolute())
import pip

def import_or_install(package):
    try:
        __import__(package)
    except ImportError:
        pip.main(['install', package])      
import sys, os

#for p in sys.path:
#    print(p)

import pylatexenc