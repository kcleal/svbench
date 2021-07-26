from setuptools import setup, find_packages

setup(
    name="svbench",
    version='0.5.1',
    python_requires='>=3.7',
    install_requires=[
            'numpy',
            'pandas',
            'networkx>=2.4',
            'ncls',
            'pyvcf'
        ],
    packages=find_packages(where="."),
    include_package_data=True,
    # entry_points='''
    #     [console_scripts]
    #     svbench=svbench.io_tools:apply_model
    # ''',
)
