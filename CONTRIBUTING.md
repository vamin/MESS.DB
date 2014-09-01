# Coding Standards #

All code commited should, at minimum, adhere to the PEP8 style guide. Empty
lines in indented code blocks should be indented to facilitate copy/pasting
code into e.g. python terminals. The [pep8][] python package should be used to
check code conformance by running `pep8 --ignore=W293`. 

It is also recommended to run [pylint][] on new code, though full compliance is
neither expected nor necessary. Use `pylint --disable=C0303` to silence errors
about indents in empty lines.

[pep8]: https://pypi.python.org/pypi/pep8
[pylint]: http://www.pylint.org/

# Pull Requests #

For now, refer to
[these guidelines](http://contribute.jquery.org/commits-and-pull-requests/). If
that doesn't make sense to you, feel free to [contact me](http://vamin.net) and
I'll be happy to work with you.
