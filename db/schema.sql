/*
If your sqlite3 installation supports foreign keys,
uncommenting the foreign key declarations is recommended.
*/

/*
Describing molecules, their names, and where to get them.
*/

--DROP TABLE molecule;
CREATE TABLE IF NOT EXISTS
    molecule( --contains molecule name and fundamental identifiers
        inchikey TEXT PRIMARY KEY,
        inchi TEXT,
        smiles TEXT,
        formula TEXT,
        iupac TEXT
    );

--DROP TABLE molecule_synonym;
CREATE TABLE IF NOT EXISTS
    molecule_synonym( --any synonyms, common names, identifiers, etc.
        inchikey TEXT NOT NULL,
        --FOREIGN KEY(inchikey) REFERENCES molecule(inchikey),
        name TEXT NOT NULL,
        UNIQUE(inchikey, name)
    );
CREATE INDEX IF NOT EXISTS ix_molecule_synonym_name ON molecule_synonym (name);

--DROP TABLE state;
CREATE TABLE IF NOT EXISTS
    state( --descriptors for various electronic/conformational states that can be applied to a molecule
        state_id INTEGER PRIMARY KEY,
        name TEXT NOT NULL, --e.g. anion, cation, triplet, protonated, etc.
        UNIQUE(name)
    );

--DROP TABLE source;
CREATE TABLE IF NOT EXISTS
    source( --catalogs where molecules exist and might be purchased
        source_id INTEGER PRIMARY KEY,
        name TEXT NOT NULL,
        filename TEXT NOT NULL,
        url TEXT,
        url_template TEXT,
        last_update TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_source_name ON source (name);
CREATE UNIQUE INDEX IF NOT EXISTS ux_source_filename ON source (filename);

--DROP TABLE molecule_source;
CREATE TABLE IF NOT EXISTS
    molecule_source( --links molecules to sources
        inchikey TEXT,
        --FOREIGN KEY(inchikey) REFERENCES molecule(inchikey),
        source_id INTEGER,
        --FOREIGN KEY(source_id) REFERENCES source(source_id),
        identifier TEXT
    );
CREATE INDEX IF NOT EXISTS ix_molecule_souce_identifier ON molecule_source (identifier);

/*
Describing methods to calculate properties of molecules and their relationships.
*/

--DROP TABLE level;
CREATE TABLE IF NOT EXISTS
    level( --level of theory, used to classify methods
        level_id INTEGER PRIMARY KEY,
        name TEXT NOT NULL --e.g. empirical, semiempirical, dft, hf, etc
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_level_name ON level (name);

--DROP TABLE program;
CREATE TABLE IF NOT EXISTS
    program( --programs used to calculate properties
        program_id INTEGER PRIMARY KEY,
        name TEXT NOT NULL,
        version REAL NOT NULL,
        url TEXT,
        UNIQUE(name, version)
    );
--CREATE UNIQUE INDEX IF NOT EXISTS ux_program_name ON program (name);

--DROP TABLE method;
CREATE TABLE IF NOT EXISTS
    method( --methods for calculating properties
        method_id INTEGER PRIMARY KEY,
        level_id INTEGER,
        --FOREIGN KEY(level_id) REFERENCES level(level_id),
        program_id INTEGER,
        --FOREIGN KEY(program_id) REFERENCES program(program_id),
        geop INTEGER, --flag, indicates whether method generates new geometry
        name TEXT NOT NULL,
        note TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_method_name ON method (name);

--DROP TABLE parameter
CREATE TABLE IF NOT EXISTS
    parameter( --description of parameters that are specified in methods
        parameter_id INTEGER PRIMARY KEY,
        name TEXT NOT NULL,
        description TEXT
    );
CREATE UNIQUE INDEX IF NOT EXISTS ux_parameter_name ON parameter (name);

--DROP TABLE method_parameter;
CREATE TABLE IF NOT EXISTS
    method_parameter( --links methods to their parameters and records the parameter setting
        method_id INTEGER NOT NULL,
        parameter_id INTEGER NOT NULL,
        setting TEXT,
        UNIQUE(method_id, parameter_id)
    );

--DROP TABLE method_edge;
CREATE TABLE IF NOT EXISTS
    method_edge( --edges in DAG describing method relationships
        method_edge_id INTEGER PRIMARY KEY,
        parent_method_id INTEGER NOT NULL,
        --FOREIGN KEY(parent_method_id) REFERENCES method(method_id),
        child_method_id INTEGER NOT NULL,
        --FOREIGN KEY(child_method_id) REFERENCES method(method_id),
        UNIQUE(parent_method_id, child_method_id)
    );

--DROP TABLE method_path;
CREATE TABLE IF NOT EXISTS
    method_path( --valid paths, e.g. ordered collections of edges, and their lengths
        method_path_id INTEGER PRIMARY KEY,
        length INTEGER NOT NULL
    );

--DROP TABLE method_path_edge;
CREATE TABLE IF NOT EXISTS
    method_path_edge( --links paths to their edges and their order (i.e. 'distance')
        method_path_id INTEGER NOT NULL,
        --FOREIGN KEY(method_path_id) REFERENCES method_path(method_path_id),
        method_edge_id INTEGER NOT NULL,
        --FOREIGN KEY(method_edge_id) REFERENCES method_edge(method_edge_id),
        distance INTEGER NOT NULL,
        UNIQUE(method_path_id, method_edge_id)
    );

--DROP TABLE method_path_parent;
CREATE TABLE IF NOT EXISTS
    method_path_parent( --this table tracks path parents and the ending method id for convenience of some queries
        method_path_id INTEGER NOT NULL,
        --FOREIGN KEY(method_path_id) REFERENCES method_path(method_path_id),
        parent_method_path_id INTEGER NOT NULL,
        --FOREIGN KEY(parent_method_path_id) REFERENCES method_path(method_path_id),
        method_id INTEGER NOT NULL, --method id of the lastest method in then path
        --FOREIGN KEY(method_id) REFERENCES method(method_id),
        UNIQUE(method_path_id, parent_method_path_id, method_id)
    );

/*
Describing properties of molecules.
*/

--DROP TABLE property;
CREATE TABLE IF NOT EXISTS
    property( --describes properties and their formats
        property_id INTEGER PRIMARY KEY,
        name TEXT NOT NULL,
        description TEXT NOT NULL,
        format TEXT NOT NULL,
        UNIQUE(name, description, format)
    );

--DROP TABLE molecule_state_method_property;
CREATE TABLE IF NOT EXISTS
    molecule_state_method_property( --links molecule, state, method_path, and property; assigns value (i.e. 'result')
        inchikey TEXT NOT NULL,
        --FOREIGN KEY(inchikey) REFERENCES molecule(inchikey),
        state_id INTEGER NOT NULL,
        --FOREIGN KEY(state_id) REFERENCES state(state_id),
        method_path_id INTEGER NOT NULL,
        --FOREIGN KEY(method_path_id) REFERENCES method(method_path_id),
        property_id INTEGER NOT NULL,
        --FOREIGN KEY(property_id) REFERENCES property(property_id),
        result BLOB,
        UNIQUE(inchikey, state_id, method_path_id, property_id)
    );
