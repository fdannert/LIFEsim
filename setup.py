import os
from setuptools import setup

def read(rel_path: str) -> str:
    here = os.path.abspath(os.path.dirname(__file__))
    # intentionally *not* adding an encoding option to open, See:
    #   https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    with open(os.path.join(here, rel_path)) as fp:
        return fp.read()

def get_version(rel_path: str) -> str:
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")

setup(
    name='LIFEsim',
    version=get_version("lifesim/__init__.py"),
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
    install_requires=['astropy>=5.2.1',
                      'matplotlib>=3.7.0',
                      'numpy>=1.24.2',
                      'pandas>=1.5.3',
                      'PyQt5==5.15.4',
                      'tqdm>=4.64.1',
                      'tables>=3.8.0',
                      'GitPython>=3.1.32'],
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
