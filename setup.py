try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
import sys, os

config = {
    'descripton': 'Molecular Electronic Structure and other Stuff DB.',
    'url': 'https://github.com/vamin/messdb',
    'author': 'Victor Amin',
    'author_email': 'victor.amin@gmail.com',
    'version': '0.1',
    'install_requires': ['pybel'],
    'packages': ['mess'],
    'scripts': [],
    'name': 'messdb',
    'entry_points': {'console_scripts': ['mess=mess.__main__:main']}
}

setup(**config)