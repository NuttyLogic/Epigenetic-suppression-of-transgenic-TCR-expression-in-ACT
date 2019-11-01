from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(name='DetectDiscord',
      version='0.0.1',
      description='Bisulfite Vector Integration Detection',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/NuttyLogic/RebisMethylation',
      project_urls={'Documentation': ''},
      author='Colin P. Farrell',
      author_email='colinpfarrell@gmail.com',
      license='GPLv3',
      packages=[],
      classifiers=['Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   'Programming Language :: Python :: 3.8'],
      platforms=["Linux", "Mac OS-X", "Unix"],
      requires=['pysam', 'numpy', 'tqdm'],
      install_requires=['pysam==0.15.2', 'numpy>=1.16.3', 'tqdm>=4.31.1'],
      entry_points={'console_scripts': ['BSBolt = BSBolt.__main__:launch_bsb']},
      python_requires='>=3.6',
      test_suite='tests',
      include_package_data=True,
      cmdclass={'develop': DevelopCmd,
                'build_py': BuildCmd,
                'bdist_wheel': bdist_wheel},
      zip_safe=False
      )
