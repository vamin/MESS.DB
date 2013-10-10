--insert program
INSERT OR IGNORE INTO program
(name, version, url)
VALUES
('MOPAC', '2012', 'http://openmopac.net/');

--make sure theory level is in there
INSERT OR IGNORE INTO level
(name)
VALUES
('semiempirical');

--insert method row
INSERT OR IGNORE INTO method
(level_id, program_id, geop, name)
SELECT
level.level_id, program.program_id, 1, 'pm7_mopac2012' 
FROM level, program
WHERE level.name='semiempirical'
AND program.name='MOPAC'
AND program.version='2012';

--insert parameters
--PM7 PRECISE LARGE ALLVEC BONDS LOCALIZE AUX GNORM=0.0 T=1W
INSERT OR IGNORE INTO parameter
(name)
VALUES
('PM7');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, ''
FROM method, parameter
WHERE method.name='pm7_mopac2012'
AND parameter.name='PM7';

INSERT OR IGNORE INTO parameter
(name)
VALUES
('PRECISE');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, ''
FROM method, parameter
WHERE method.name='pm7_mopac2012'
AND parameter.name='PRECISE';

INSERT OR IGNORE INTO parameter
(name)
VALUES
('LARGE');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, ''
FROM method, parameter
WHERE method.name='pm7_mopac2012'
AND parameter.name='LARGE';

INSERT OR IGNORE INTO parameter
(name)
VALUES
('ALLVEC');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, ''
FROM method, parameter
WHERE method.name='pm7_mopac2012'
AND parameter.name='ALLVEC';

INSERT OR IGNORE INTO parameter
(name)
VALUES
('BONDS');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, ''
FROM method, parameter
WHERE method.name='pm7_mopac2012'
AND parameter.name='BONDS';

INSERT OR IGNORE INTO parameter
(name)
VALUES
('LOCALIZE');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, ''
FROM method, parameter
WHERE method.name='pm7_mopac2012'
AND parameter.name='LOCALIZE';

INSERT OR IGNORE INTO parameter
(name)
VALUES
('AUX');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, ''
FROM method, parameter
WHERE method.name='pm7_mopac2012'
AND parameter.name='AUX';

INSERT OR IGNORE INTO parameter
(name)
VALUES
('GNORM');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, '0.0'
FROM method, parameter
WHERE method.name='pm7_mopac2012'
AND parameter.name='GNORM';

INSERT OR IGNORE INTO parameter
(name)
VALUES
('T');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, '1W'
FROM method, parameter
WHERE method.name='pm7_mopac2012'
AND parameter.name='T';