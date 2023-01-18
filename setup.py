# Author: Tarik Salameh

"""
Setup script for HiCLift.

This is a free software under GPLv3. Therefore, you can modify, redistribute
or even mix it with other GPL-compatible codes. See the file LICENSE
included with the distribution for more details.

"""
import os, sys, HiCLift, glob
import setuptools

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if (sys.version_info.major!=3) or (sys.version_info.minor<7):
    print('PYTHON 3.7+ IS REQUIRED. YOU ARE CURRENTLY USING PYTHON {}'.format(sys.version.split()[0]))
    sys.exit(2)

# Guarantee Unix Format
for src in glob.glob('scripts/*'):
    text = open(src, 'r').read().replace('\r\n', '\n')
    open(src, 'w').write(text)

setuptools.setup(
    name = 'HiCLift',
    version = HiCLift.__version__,
    author = HiCLift.__author__,
    author_email = 'wangxiaotao686@gmail.com',
    url = 'https://github.com/XiaoTaoWang/HiCLift',
    description = 'A fast and efficient tool for converting chromatin interaction data between genome assemblies',
    keywords = 'Hi-C LiftOver',
    packages = setuptools.find_packages(),
    package_data = {
        '': ['data/*']
    },
    scripts = glob.glob('scripts/*'),
    long_description = read('README.rst'),
    long_description_content_type='text/x-rst',
    classifiers = [
        'Programming Language :: Python :: 3 :: Only',
        'Operating System :: POSIX',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ]
    )

