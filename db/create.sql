CREATE TABLE std_inchi_v1 (
    id BIGSERIAL PRIMARY KEY,
    inchis TEXT NOT NULL,
    inchik TEXT NOT NULL,
    time_created TIMESTAMP NOT NULL DEFAULT NOW(),

    CONSTRAINT uq_inchis UNIQUE (inchis)
);

CREATE TABLE smiles_non_iso (
    id BIGSERIAL PRIMARY KEY,
    std_inchi_v1_id BIGINT NOT NULL,
    smiles TEXT NOT NULL,
    inchis TEXT NOT NULL,
    inchik TEXT NOT NULL,
    hac SMALLINT NOT NULL,
    rac SMALLINT NOT NULL,
    time_created TIMESTAMP NOT NULL DEFAULT NOW(),

    CONSTRAINT fk_non_iso2std_inchi_v1 FOREIGN KEY (std_inchi_v1_id) REFERENCES std_inchi_v1(id) ON DELETE CASCADE
);

CREATE TABLE smiles_iso (
    id BIGSERIAL PRIMARY KEY,
    smiles_non_iso_id BIGINT NOT NULL,
    smiles TEXT NOT NULL,
    time_created TIMESTAMP NOT NULL DEFAULT NOW(),

    CONSTRAINT fk_iso2non_iso FOREIGN KEY (smiles_non_iso_id) REFERENCES smiles_non_iso(id) ON DELETE CASCADE
);

CREATE TABLE source (
    id SERIAL PRIMARY KEY,
    name TEXT NOT NULL,
    version TEXT NOT NULL,
    time_created TIMESTAMP NOT NULL DEFAULT NOW(),

    CONSTRAINT uq_source UNIQUE (name, version)
);

CREATE TABLE mol_source (
    id SERIAL PRIMARY KEY,
    source_id INTEGER NOT NULL,
    smiles_non_iso_id BIGINT,
    smiles_iso_id BIGINT,
    source_code TEXT,
    source_smiles TEXT,
    time_created TIMESTAMP NOT NULL DEFAULT NOW(),

    CONSTRAINT fk_mol_source2source FOREIGN KEY (source_id) REFERENCES source(id) ON DELETE CASCADE,
    CONSTRAINT fk_mol_source2smiles_iso FOREIGN KEY (smiles_iso_id) REFERENCES smiles_iso(id) ON DELETE CASCADE,
    CONSTRAINT fk_mol_source2smiles_non_iso FOREIGN KEY (smiles_non_iso_id) REFERENCES smiles_non_iso(id) ON DELETE CASCADE
);
