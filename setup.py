from setuptools import setup, find_packages

setup(
    name="svbench",
    version='0.5.2',
    license='MIT',
    python_requires='>=3.7',
    install_requires=[
            'numpy',
            'pandas',
            'networkx>=2.4',
            'ncls',
            'pyvcf',
            'matplotlib'
        ],
    packages=find_packages(where="."),
    include_package_data=True,
    # entry_points='''
    #     [console_scripts]
    #     svbench=svbench.io_tools:apply_model
    # ''',
)
