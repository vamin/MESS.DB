INSERT OR IGNORE INTO source
(name, dirname, url, url_template, last_update)
VALUES
(
    'FDA-approved drugs (via DSSTOX)',
    'fda',
    null,
    null,
    null
);
UPDATE source SET
    url='http://www.epa.gov/nheerl/dsstox/',
    url_template=null,
    last_update='2012-07-25'
WHERE dirname='fda';