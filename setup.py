from setuptools import setup

setup(
    name= 'CAVA',
    version = '1.2.3',
    description = 'CAVA (Clinical Annotation of VAriants)',
    url = 'https://github.com/RahmanTeam/CAVA',
    author = 'RahmanTeam',
    author_email = 'rahmanlab@icr.ac.uk',
    license = 'MIT',
    packages=['cava_', 'ensembldb'],
    scripts=[
        'bin/CAVA.py',
        'bin/cava',
        'bin/EnsemblDB.py',
        'bin/ensembl_db',
        'bin/dbSNPDB.py',
        'bin/dbsnp_db'
    ],
    zip_safe=False,
    include_package_data=True
)
