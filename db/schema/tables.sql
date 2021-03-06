/*
Describing molecules, their names, and where to get them
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
  molecule_synonym( --any synonyms, common names, etc.
    inchikey TEXT NOT NULL,
    --FOREIGN KEY(inchikey) REFERENCES molecule(inchikey),
    name TEXT NOT NULL,
    UNIQUE(inchikey, name)
  );
CREATE INDEX IF NOT EXISTS ix_molecule_synonym_name ON molecule_synonym (name);

--DROP TABLE molecule_fingerprint;
CREATE TABLE IF NOT EXISTS
  molecule_fingerprint(
    inchikey TEXT NOT NULL,
    --FOREIGN KEY(inchikey) REFERENCES molecule(inchikey),
    name TEXT NOT NULL,
    settings TEXT,
    fingerprint TEXT NOT NULL,
    method_path_id INTEGER,
    --FOREIGN KEY(method_path_id) REFERENCES method_path(method_path_id),
    UNIQUE(inchikey, name, settings, method_path_id)
  );

--DROP TABLE source;
CREATE TABLE IF NOT EXISTS
  source( --catalogs where molecules exist and might be purchased
    source_id INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    dirname TEXT NOT NULL,
    url TEXT,
    url_template TEXT,
    citation TEXT,
    last_update TEXT
  );
CREATE UNIQUE INDEX IF NOT EXISTS ux_source_name ON source (name);
CREATE UNIQUE INDEX IF NOT EXISTS ux_source_dirname ON source (dirname);

--DROP TABLE molecule_source;
CREATE TABLE IF NOT EXISTS
  molecule_source( --links molecules to sources
    inchikey TEXT,
    --FOREIGN KEY(inchikey) REFERENCES molecule(inchikey),
    source_id INTEGER,
    --FOREIGN KEY(source_id) REFERENCES source(source_id),
    identifier TEXT,
    UNIQUE(inchikey, source_id, identifier)
  );
CREATE INDEX IF NOT EXISTS ix_molecule_souce_identifier
  ON molecule_source (identifier);

/*
Describing methods to calculate properties of molecules
*/

--DROP TABLE program;
CREATE TABLE IF NOT EXISTS
  program( --programs used to calculate properties
    program_id INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    version TEXT NOT NULL,
    url TEXT,
    citation TEXT,
    UNIQUE(name, version)
  );

--DROP TABLE parameter
CREATE TABLE IF NOT EXISTS
  parameter( --description of parameters that are specified in programs
    parameter_id INTEGER PRIMARY KEY,
    name TEXT NOT NULL
  );
CREATE UNIQUE INDEX IF NOT EXISTS ux_parameter_name ON parameter (name);

--DROP TABLE method;
CREATE TABLE IF NOT EXISTS
  method( --methods for calculating properties
    method_id INTEGER PRIMARY KEY,
    program_id INTEGER,
    --FOREIGN KEY(program_id) REFERENCES program(program_id),
    name TEXT,
    shortdesc TEXT,
    citation TEXT,
    geop INTEGER, --flag, indicates whether method generates new geometry
    hash TEXT NOT NULL --sha1 hash of method name and parameter settings
  );
CREATE UNIQUE INDEX IF NOT EXISTS ux_method_hash ON method (hash);

--DROP TABLE method_parameter;
CREATE TABLE IF NOT EXISTS
  method_parameter( --links methods to parameters and records parameter setting
    method_id INTEGER NOT NULL,
    --FOREIGN KEY(method_id) REFERENCES method(method_id),
    parameter_id INTEGER NOT NULL,
    --FOREIGN KEY(parameter_id) REFERENCES parameter(parameter_id),
    setting TEXT,
    UNIQUE(method_id, parameter_id, setting)
  );

/*
Describing method relationships.
*/

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
  method_path( --valid paths (ordered collections of edges) and their lengths
    method_path_id INTEGER PRIMARY KEY,
    length INTEGER NOT NULL
  );

--DROP TABLE method_path_edge;
CREATE TABLE IF NOT EXISTS
  method_path_edge( --links paths to their edges and their order ('distance')
    method_path_id INTEGER NOT NULL,
    --FOREIGN KEY(method_path_id) REFERENCES method_path(method_path_id),
    method_edge_id INTEGER NOT NULL,
    --FOREIGN KEY(method_edge_id) REFERENCES method_edge(method_edge_id),
    distance INTEGER NOT NULL,
    UNIQUE(method_path_id, method_edge_id)
  );

/*
Describing properties of molecules
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

--DROP TABLE molecule_method_property;
CREATE TABLE IF NOT EXISTS
  molecule_method_property( --link value (result) to  molecule, method_path,
                            --and property
    inchikey TEXT NOT NULL,
    --FOREIGN KEY(inchikey) REFERENCES molecule(inchikey),
    method_path_id INTEGER NOT NULL,
    --FOREIGN KEY(method_path_id) REFERENCES method_path(method_path_id),
    property_id INTEGER NOT NULL,
    --FOREIGN KEY(property_id) REFERENCES property(property_id),
    units TEXT NOT NULL,
    result BLOB,
    UNIQUE(inchikey, method_path_id, property_id, units)
  );
