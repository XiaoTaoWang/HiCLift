# Author: Tarik Salameh

"""
Setup script for pairLiftOver.

This is a free software under GPLv3. Therefore, you can modify, redistribute
or even mix it with other GPL-compatible codes. See the file LICENSE
included with the distribution for more details.

"""
import os, sys, pairLiftOver, glob
import setuptools

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if (sys.version_info.major!=3) or (sys.version_info.minor<6):
    print('PYTHON 3.6+ IS REQUIRED. YOU ARE CURRENTLY USING PYTHON {}'.format(sys.version.split()[0]))
    sys.exit(2)

# Guarantee Unix Format
for src in glob.glob('scripts/*'):
    text = open(src, 'r').read().replace('\r\n', '\n')
    open(src, 'w').write(text)

setuptools.setup(
    name = 'pairLiftOver',
    version = pairLiftOver.__version__,
    author = pairLiftOver.__author__,
    author_email = 'wangxiaotao686@gmail.com',
    url = 'https://github.com/XiaoTaoWang/pairLiftOver',
    description = 'Convert genomic coordinates of contact pairs from one assembly to another.',
    keywords = 'liftover pairs Hi-C 4DN',
    packages = setuptools.find_packages(),
    package_data = {
        '': ['data/*']
    },
    scripts = glob.glob('scripts/*'),
    long_description = read('README.rst'),
    long_description_content_type='text/x-rst',
    classifiers = [
        'Programming Language :: Python :: 3.6',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ]
    )

