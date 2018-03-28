'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Jan 27, 2017
@author: Niels Lubbes

https://python-packaging.readthedocs.io/en/latest/minimal.html
https://pypi.python.org/pypi?%3Aaction=list_classifiers
'''


from setuptools import setup


setup( name = 'moebius_aut',
       version = '0',
       description = 'Moebius automorphisms of surfaces',
       classifiers = [
           'Development Status :: 3 - Alpha',
           'License :: OSI Approved :: MIT License',
           'Programming Language :: Python :: 2',
           'Programming Language :: Python :: 3',
           'Topic :: Scientific/Engineering :: Mathematics',
           ],
      keywords = 'surfaces automorphisms circles',
      url = 'http://github.com/niels-lubbes/moebius_aut',
      author = 'Niels Lubbes',
      license = 'MIT',
      package_dir = {'': 'src'},
      packages = ['moebius_aut'],
      test_suite = 'nose.collector',
      tests_require = ['nose'],
      entry_points = {
          'console_scripts': ['run-moebius_aut=moebius_aut.__main__:main'],
      },
      include_package_data = True,
      zip_safe = False )


