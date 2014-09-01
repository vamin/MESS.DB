# Running Locally #

Each MESS.DB directory is completely self-contained, including molecule
sources, backups, logs, data, and the scripts used to generate them. No
installation is required if you run your scripts within the directory.

To run locally, select a directory you want to work in on a filesystem that has
enough room to store all the data you will generate. Then, clone the MESS.DB
[repository][] to that directory.


From within a MESS.DB directory, you can run mess with `python mess` or
`bin/mess`. For convenience, you can alias the mess executable:

```bash
alias mess=${PWD}/bin/mess
```

and simply run `mess`.

This project structure has the advantage of being completely portable. Simply
copying the entire MESS.DB directory to a new location preserves all data.

[repository]: http://github.com/vamin/MESS.DB

# Installing as a Module #

In the future, MESS.DB may provide a setup.py which will make `mess` available
in the system path. Before that can happen, MESS.DB will need to be modified to
handle separation of MESS.DB scripts and database files.

One of the goals of MESS.DB is to provide a system that maintains total
accountability for every datum generated so that all calculations can be
repeated exactly. For now, running locally is the recommended option, as it
keeps the apparatus used to generate data together with that data.
