from setuptools import setup


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='BSIntegration',
      version='0.0.1',
      description='Bisulfite Sequencing Vector Integration Discovery Tool',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='',
      project_urls={'Documentation': ''},
      author='Colin P. Farrell',
      author_email='colinpfarrell@gmail.com',
      license='MIT',
      packages=['IntegrationSiteSearch'],
      entry_points={'console_scripts': ['BSBolt = BSBolt.__main__:launch_bsb']},
      python_requires='>=3.6',
      test_suite='tests',
      include_package_data=True,
      )
