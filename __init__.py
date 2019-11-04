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
      python_requires='>=3.6',
      test_suite='tests',
      include_package_data=True,
