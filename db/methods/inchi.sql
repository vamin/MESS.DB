INSERT OR IGNORE INTO program
(name, version, url)
VALUES
('Open Babel', '', 'http://openbabel.org/wiki/Main_Page');

INSERT OR IGNORE INTO method
(program_id, name, note)
SELECT
program.program_id, 'inchi', 'initial import'
FROM program
WHERE
program.version=''
AND
program.name='Open Babel';

INSERT OR IGNORE INTO method_edge
(parent_method_id, child_method_id)
SELECT
method_id, method_id
FROM method
WHERE
name='inchi';
