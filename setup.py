from setuptools import setup
setup(
    name='julia_rings',
    version='0.1.0',    
    description='Ring counting in amorphous structure implemented in Julia',
    url='https://github.com/MorrowChem/julia_rings',
    author='Joe Morrow',
    author_email='joe.morrow@queens.ox.ac.uk',
    license='BSD 2-clause',
    packages=['julia_rings'],
    install_requires=[
                      'ase',
                      'numpy',
                      'julia',                
                      ],

    classifiers=[
        'Development Status :: 2 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux :: OSX',        
        'Programming Language :: Python :: 3'
    ],
)