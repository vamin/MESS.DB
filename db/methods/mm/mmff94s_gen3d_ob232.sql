--insert program
INSERT OR IGNORE INTO program
(name, version, url)
VALUES
('Open Babel', '2.3.2', 'http://openbabel.org/wiki/Main_Page');

--make sure theory level is in there
INSERT OR IGNORE INTO level
(name)
VALUES
('mm');

--insert method row
INSERT OR IGNORE INTO method
(level_id, program_id, geop, name)
SELECT
level.level_id, program.program_id, 1, 'mmff94s_gen3d_ob232' 
FROM level, program
WHERE level.name='mm'
AND program.name='Open Babel'
AND program.version='2.3.2';

--insert parameters
INSERT OR IGNORE INTO parameter
(name)
VALUES
('forcefield');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, 'mmmff94s'
FROM method, parameter
WHERE method.name='mmff94s_gen3d_ob232'
AND parameter.name='forcefield';

INSERT OR IGNORE INTO parameter
(name)
VALUES
('steps');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, '1024'
FROM method, parameter
WHERE method.name='mmff94s_gen3d_ob232'
AND parameter.name='steps';

INSERT OR IGNORE INTO parameter
(name)
VALUES
('localopt');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, ''
FROM method, parameter
WHERE method.name='mmff94s_gen3d_ob232'
AND parameter.name='localopt';

INSERT OR IGNORE INTO parameter
(name)
VALUES
('make3D');

INSERT OR IGNORE INTO method_parameter
(method_id, parameter_id, setting)
SELECT
method.method_id, parameter.parameter_id, ''
FROM method, parameter
WHERE method.name='mmff94s_gen3d_ob232'
AND parameter.name='make3D';
