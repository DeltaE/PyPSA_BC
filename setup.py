from setuptools import setup, find_packages

setup(
    name='pypsa_bc', 
    version='1.0',                
    author='Elias Islam, Pierre M',
    author_email='elias_islam@sfu.ca',  # Corresponding author email
    description='A combined package tool to serve the requirement for the scripts of PypSA-BC model',
    long_description=open('README.md').read(),  
    long_description_content_type='text/markdown', 
    url='https://github.com/DeltaE/PyPSA_BC', 
    packages=find_packages(),    # Automatically find all packages
    install_requires=[           # Specify dependencies here
        # Example: 'numpy>=1.20.0', 'pandas>=1.3.0'
    ],
    # entry_points={               # For creating CLI commands
    #     'console_scripts': [
    #         'bccm=bc_combined_modelling.cli:main',  # Maps command `bccm` to `bc_combined_modelling/cli.py:main`
    #     ],
    # },
    classifiers=[                # Optional: Metadata for PyPI
        'Programming Language :: Python :: 3.12',
        'License :: MIT License',
        'Operating System :: Linux (Ubuntu)',
    ],
    python_requires='>=3.10',     # Minimum Python version
    project_urls={              
        'Research Lab': 'https://deltaeplus.sfu.ca',  
        'Documentation': 'https://github.com/DeltaE/PyPSA_BC/wiki',
    },
)
