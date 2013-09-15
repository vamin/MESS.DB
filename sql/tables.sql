--if your sqlite3 installation supports foreign keys,
--uncommenting the foreign key declarations is recommended

--DROP TABLE molecule
CREATE TABLE IF NOT EXISTS
    molecule( --contains molecule names and empirical parameters
        inchikey TEXT PRIMARY KEY,
        inchi TEXT,
        smiles TEXT,
        cactvs_hashisy TEXT,
        formula TEXT,
        iupac TEXT
    );

--DROP TABLE molecule_synonym
CREATE TABLE IF NOT EXISTS
    molecule_synonym(
        inchikey TEXT,
        --FOREIGN KEY(inchikey) REFERENCES molecule(inchikey),
        name TEXT,
        UNIQUE(inchikey, name)
    );
CREATE INDEX IF NOT EXISTS ix_molecule_synonym_name ON molecule_synonym (name);

--DROP TABLE source
CREATE TABLE IF NOT EXISTS
    source(
        source_id INTEGER PRIMARY KEY,
        name TEXT,
        filename TEXT,
        url TEXT,
        url_template TEXT,
        last_update TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_source_name ON source (name);
CREATE UNIQUE INDEX IF NOT EXISTS ux_source_filename ON source (filename);

--DROP TABLE molecule_source
CREATE TABLE IF NOT EXISTS
    molecule_source(
        inchikey TEXT,
        --FOREIGN KEY(inchikey) REFERENCES molecule(inchikey),
        source_id INTEGER,
        --FOREIGN KEY(source_id) REFERENCES source(source_id),
        identifier TEXT
    );
CREATE INDEX IF NOT EXISTS ix_molecule_souce_identifier ON molecule_source (identifier);

--DROP TABLE program
CREATE TABLE IF NOT EXISTS
    program(
        program_id INTEGER PRIMARY KEY,
        name TEXT,
        version REAL
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_program_name ON program (name);

--DROP TABLE program
CREATE TABLE IF NOT EXISTS
    job_type(
        job_id INTEGER PRIMARY KEY,
        name TEXT,
        geop INT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_job_type_name ON job_type (name);

--DROP TABLE theory_level
CREATE TABLE IF NOT EXISTS
    theory_level(
        theory_id INTEGER PRIMARY KEY,
        name TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_theory_level_name ON theory_level (name);

--DROP TABLE mm_method
CREATE TABLE IF NOT EXISTS
    mm_method(
        mm_id INTEGER PRIMARY KEY,
        name TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_mm_method_name ON mm_method (name);

--DROP TABLE md_method
CREATE TABLE IF NOT EXISTS
    md_method(
        md_id INTEGER PRIMARY KEY,
        name TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_md_method_name ON md_method (name);

--DROP TABLE empirical_method
CREATE TABLE IF NOT EXISTS
    empirical_method(
        empirical_id INTEGER PRIMARY KEY,
        name TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_empirical_method_name ON empirical_method (name);

--DROP TABLE semiempirical_method
CREATE TABLE IF NOT EXISTS
    semiempirical_method(
        semiempirical_id INTEGER PRIMARY KEY,
        name TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_semiempirical_method_name ON semiempirical_method (name);

--DROP TABLE mp_method
CREATE TABLE IF NOT EXISTS
    mp_method(
        mp_id INTEGER PRIMARY KEY,
        name TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_mp_method_name ON mp_method (name);

--DROP TABLE cc_method
CREATE TABLE IF NOT EXISTS
    cc_method(
        cc_id INTEGER PRIMARY KEY,
        name TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_cc_method_name ON cc_method (name);

--DROP TABLE basis_set
CREATE TABLE IF NOT EXISTS
    basis_set(
        basis_id INTEGER PRIMARY KEY,
        name TEXT,
        valence INTEGER,
        zeta INTEGER,
        p INTEGER,
        d INTEGER,
        f INTEGER,
        diff INTEGER
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_basis_set_name ON basis_set (name);

--DROP TABLE ecp
CREATE TABLE IF NOT EXISTS
    ecp(
        ecp_id INTEGER PRIMARY KEY,
        name TEXT,
        core_size INTEGER
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_ecp_name ON ecp (name);

--DROP TABLE dispersion
CREATE TABLE IF NOT EXISTS
    dispersion(
        dispersion_id INTEGER PRIMARY KEY,
        name TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_dispersion_name ON dispersion (name);

--DROP TABLE functional
CREATE TABLE IF NOT EXISTS
    functional(
        functional_id INTEGER PRIMARY KEY,
        name TEXT,
        rung REAL
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_functional_name ON functional (name);

--DROP TABLE method
CREATE TABLE IF NOT EXISTS
    method(
        method_id INTEGER PRIMARY KEY,
        theory_id INTEGER,
        --FOREIGN KEY(theory_id) REFERENCES theory_level(theory_id),
        mm_id INTEGER,
        --FOREIGN KEY(mm_id) REFERENCES mm_method(mm_id),
        md_id INTEGER,
        --FOREIGN KEY(md_id) REFERENCES md_method(md_id),
        empirical_id INTEGER,
        --FOREIGN KEY(empirical_id) REFERENCES empirical_method(empirical_id),
        semiempirical_id INTEGER,
        --FOREIGN KEY(semiempirical_id) REFERENCES semiempirical_method(semiempirical_id),
        mp_id INTEGER,
        --FOREIGN KEY(mp_id) REFERENCES mp_method(mp_id),
        cc_id INTEGER,
        --FOREIGN KEY(cc_id) REFERENCES cc_method(cc_id),
        basis_id INTEGER,
        --FOREIGN KEY(basis_id) REFERENCES basis_set(basis_id),
        ecp_id INTEGER,
        --FOREIGN KEY(ecp_id) REFERENCES ecp(ecp_id),
        functional_id INTEGER,
        --FOREIGN KEY(functional_id) REFERENCES functional(functional_id),
        dispersion_id INTEGER,
        --FOREIGN KEY(dispersion_id) REFERENCES dispersion(dispersion_id),
        job_id INTEGER,
        --FOREIGN KEY(run_id) REFERENCES job_type(job_id)
        program_id INTEGER,
        --FOREIGN KEY(program_id) REFERENCES program(program_id),
        parent_method_id INTEGER,
        --FOREIGN KEY(program_id) REFERENCES program(program_id),
        name TEXT,
        note TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_method_name ON method (name);

--DROP TABLE property
CREATE TABLE IF NOT EXISTS
    property(
        property_id INTEGER PRIMARY KEY,
        name TEXT,
        description TEXT,
        format TEXT
    );

--DROP TABLE molecule_property
CREATE TABLE IF NOT EXISTS
    molecule_property(
        inchikey TEXT PRIMARY KEY,
        --FOREIGN KEY(inchikey) REFERENCES molecule(inchikey),
        property_id INTEGER,
        --FOREIGN KEY(property_id) REFERENCES property(property_id),
        method_id INTEGER,
        --FOREIGN KEY(method_id) REFERENCES method(method_id)
        result BLOB
    );

--DROP TABLE blacklist
--CREATE TABLE IF NOT EXISTS
--    blacklist(
--        inchikey TEXT PRIMARY KEY,
--        --FOREIGN KEY(inchikey) REFERENCES molecule(inchikey),
--        method_id INTEGER,
--        --FOREIGN KEY(method_id) REFERENCES method(method_id)
--        reason TEXT
--    );

