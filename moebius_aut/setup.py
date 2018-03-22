'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Jan 27, 2017
@author: Niels Lubbes

https://python-packaging.readthedocs.io/en/latest/minimal.html
https://pypi.python.org/pypi?%3Aaction=list_classifiers
'''


from setuptools import setup


def readme():
    with open( 'README.md' ) as f:
        return f.read()


setup( name = 'moebius_aut',
       version = '0',
       description = 'Moebius automorphisms of surfaces',
       long_description = readme(),
       classifiers = [
           'Development Status :: 3 - Alpha',
           'License :: OSI Approved :: MIT License',
           'Programming Language :: Python :: 3.6',
           'Topic :: Scientific/Engineering :: Mathematics',
           ],
      keywords = 'surfaces automorphisms circles',
      url = 'http://github.com/niels-lubbes/moebius_aut',
      author = 'Niels Lubbes',
      license = 'MIT',
      package_dir = {'': 'src'},
      packages = ['moebius_aut'],
      # install_requires = ['markdown'],
      # dependency_links = ['http://github.com/niels-lubbes/linear_series/tarball/master#egg=package-1.0'],
      test_suite = 'nose.collector',
      tests_require = ['nose'],
      entry_points = {
          'console_scripts': ['run-moebius_aut=moebius_aut.__main__:main'],
      },
      include_package_data = True,
      zip_safe = False )


