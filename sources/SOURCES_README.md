# Manual Source Setup #

Each source must have its own directory containing:

1. molecule files or multi-molecule file in a format/file extension that Open
   Babel understands (run `obabel -L formats` for an exhaustive list)
2. an SQL file (called dirname.sql) describing the source for import into the
   db, like so:

        INSERT OR IGNORE INTO source
        (name, dirname, url, url_template, last_update)
        VALUES
        (
            'LONG/HUMAN READABLE NAME',
            'DIRECTORY/BASE NAME',
            null,
            null,
            null
        );
        UPDATE source SET
            url='PRIMARY SORUCE URL', --can be null
            url_template='URL TEMPLATE FOR PARTICULAR MOLECULE', --can be null, use regex in square brackets to describe where to put catalog no (e.g., [(.+), \1] for direct replace)
            last_update='DATE OF LAST UPDATE, YYYY-MM-DD'
        WHERE dirname='DIRECTORY/BASE NAME';

    For example, in sialbb/sialbb.sql:

        INSERT OR IGNORE INTO source
        (name, dirname, url, url_template, last_update)
        VALUES
        (
            'Sigma Aldrich Building Blocks',
            'sialbb',
            null,
            null,
            null
        );
        UPDATE source SET
            url='http://www.sigmaaldrich.com/',
            url_template='http://www.sigmaaldrich.com/catalog/product/[(.+?)\|(.+), \2/\1]',
            last_update='2013-08-13'
        WHERE dirname='sialbb';

# Automated Source Setup #

In the future, source setup will be automated with setup_source.sh.
