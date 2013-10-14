try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
import sys, os

config = {
    'description': 'Molecular Electronic Structure and other Stuff DB.',
    'url': 'https://github.com/vamin/messdb',
    'author': 'Victor Amin',
    'author_email': 'victor.amin@gmail.com',
    'version': '0.3',
    'install_requires': [],
    'packages': ['mess'],
    'scripts': [],
    'name': 'messdb',
    'entry_points': {'console_scripts': ['mess=mess.__main__:main']}
}

setup(**config)
