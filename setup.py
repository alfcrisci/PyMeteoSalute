#!/usr/bin/env python

from distutils.core import setup

classifiers = ['Development Status :: 5 - Production/Stable',
    'Environment :: Console',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License (GPL)',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
    'Topic :: Scientific/Engineering :: Biometeorology',
    'Topic :: Software Development :: Libraries :: Python Modules']


setup(name='pymeteosalute',
    version='0.0.1',
    description='Collection of Python libraries to estimate index for indoor/environmental thermal comfort ',
    author='Alfonso Crisci,Marco Morabito, Graziano Giuliani,Valerio Capecchi',
    author_email='alfcrisci@gmail.com',
    license = 'GNU General Public License (GPL)',
    url='http://www.ibimet.cnr.it',
    py_modules=['pymeteosalute'],
    requires = ['time', 'math', 'datetime', 'solar', 'pytz'],
    )

