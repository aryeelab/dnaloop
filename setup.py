"""
A preprocessing and QC pipeline for ChIA-PET data.
"""
from setuptools import find_packages, setup

dependencies = ['click', 'pyyaml', 'cutadapt', 'MACS2']

setup(
    name='dnaloop',
    version='0.5.5',
    url='https://github.com/aryeelab/dnaloop',
    license='BSD',
    author='Caleb Lareau and Martin Aryee',
    author_email='aryee.martin@mgh.harvard.edu',
    description='A preprocessing and QC pipeline for ChIA-PET data.',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    entry_points={
        'console_scripts': [
            'preprocess_chiapet = dnaloop.cli:main',
        ],
    },
    scripts=['dnaloop/preprocess_chiapet_fastq.sh', 'dnaloop/preprocess_chiapet_set.sh'],
    classifiers=[
        # As from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        # 'Development Status :: 1 - Planning',
        # 'Development Status :: 2 - Pre-Alpha',
        'Development Status :: 3 - Alpha',
        # 'Development Status :: 4 - Beta',
        # 'Development Status :: 5 - Production/Stable',
        # 'Development Status :: 6 - Mature',
        # 'Development Status :: 7 - Inactive',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)
