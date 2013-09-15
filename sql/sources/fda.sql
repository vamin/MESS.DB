--insert new source file only if it does not exist
INSERT OR IGNORE INTO source
(name, filename, url, url_template, last_update)
VALUES
(
    'FDA-approved drugs (via DSSTOX)',
    'fda.smi',
    null,
    null,
    null
);

--update instead of replace to preserve id,
--which is typically used as foreign key
UPDATE source SET
    url='http://www.epa.gov/nheerl/dsstox/',
    url_template=null,
    last_update='2012-07-25'
WHERE filename='fda.smi';