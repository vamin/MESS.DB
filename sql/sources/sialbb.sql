--insert new source file only if it does not exist
INSERT OR IGNORE INTO source
(name, filename, url, url_template, last_update)
VALUES
(
    'Sigma Aldrich Building Blocks',
    'sialbb.smi',
    Null,
    Null,
    Null
);

--update instead of replace to preserve id,
--which is typically used as foreign key
UPDATE source SET
    url='http://www.sigmaaldrich.com/',
    url_template='http://www.sigmaaldrich.com/catalog/product/[(.+?)\|(.+), \2/\1]',
    last_update='2013-08-13'
WHERE filename='sialbb.smi';