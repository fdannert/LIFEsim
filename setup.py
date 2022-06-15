import os
from setuptools import setup

setup(
    name='LIFEsim',
    version='0.2.17',
    description='Simulator software for the Large Interferometer For Exoplanets (LIFE)',
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    author='Felix Dannert, Maurice Ottiger & Sascha Quanz',
    author_email='fdannert@ethz.ch',
    url='https://github.com/fdannert/LIFEsim',
    project_urls={'Documentation': 'https://lifesim.readthedocs.io'},
    packages=['lifesim',
              'lifesim.core',
              'lifesim.gui',
              'lifesim.instrument',
              'lifesim.optimize',
              'lifesim.util'],
    include_package_data=True,
    install_requires=['astropy==4.2.1',
                      'matplotlib==3.4.2',
                      'numpy==1.20.3',
                      'pandas==1.2.4',
                      'PyQt5==5.15.4',
                      'tqdm==4.61.0',
                      'tables==3.6.1'
                      ],
    license='GPLv3',
    zip_safe=False,
    keywords='lifesim',
    python_requires='~=3.8',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.8',
    ]
)