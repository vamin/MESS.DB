INSERT OR IGNORE INTO program
(name, version, url)
VALUES
('xlogs', '', '');

INSERT OR IGNORE INTO level
(name)
VALUES
('empirical');

INSERT OR IGNORE INTO empirical_method
(name)
VALUES
('xlogs');

INSERT OR IGNORE INTO method
(level_id, empirical_id, program_id, name, note)
SELECT
level.level_id, empirical_method.empirical_id, program.program_id, 'xlogs', ''
FROM level, mm_method, job_type, program
WHERE level.name='empirical'
AND empirical_method.name='xlogs'
AND program.name='xlogs';