from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES
from os import sys, path
import os,shutil,re
from glob import glob
for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']

#  test
from imp import find_module
try: find_module('numpy')
except: sys.exit('### Error: python module numpy not found')
    
try: find_module('astropy')
except: sys.exit('### Error: python module astropy not found')

try: find_module('pyraf')
except: sys.exit('### Error: python module pyraf not found')

try: find_module('matplotlib')
except: sys.exit('### Error: python module matplotlib not found')


verstr = "unknown"
try:
    parentdir=os.getcwd()+'/'
    verstrline = open(parentdir+'/src/lichshane/_version.py', "rt").read()
except EnvironmentError:
    pass #  Okay, there is no version file.
else:
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    mo = re.search(VSRE, verstrline, re.M)
    if mo:
        verstr = mo.group(1)
    else:
        raise RuntimeError("unable to find version in "+parentdir+"+src/lickshane/_version.py")


setup(
    name='floyds',
    version=verstr,#'0.1.3',
    author='S. Valenti',
    author_email='stfn.valenti@gmail.com',
    scripts=['bin/lickspec'],
    url='',
    license='LICENSE.txt',
    description='lickshane is a package to reduce lick shane spectra',
    long_description=open('README.txt').read(),
    requires=['numpy','astropy','pyraf','matplotlib'],
    packages=['lickshane'],
    package_dir={'':'src'},
    package_data = {'lickshane' : ["standard/MAB/*","standard/ident/*","standard/cat/*","standard/extinction/*",\
                                   "standard/fits/*","standard/sex/*","standard/stdlist/*","standard/flux/*",\
                                   "archive/blu/*/*fits","archive/red/*/*fits",\
                                   "archive/red/arc/*fits","archive/red/*/*/id*"]}
)
