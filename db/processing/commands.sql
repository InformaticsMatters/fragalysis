INSERT INTO employees (
   employee_id,
   full_name,
   manager_id
)
VALUES
   (1, 'Michael North', NULL),
   (2, 'Megan Berry', 1),
   (3, 'Sarah Berry', 1),
   (4, 'Zoe Black', 1),
   (5, 'Tim James', 1),
   (6, 'Bella Tucker', 2),
   (7, 'Ryan Metcalfe', 2),
   (8, 'Max Mills', 2),
   (9, 'Benjamin Glover', 2),
   (10, 'Carolyn Henderson', 3),
   (11, 'Nicola Kelly', 3),
   (12, 'Alexandra Climo', 3),
   (13, 'Dominic King', 3),
   (14, 'Leonard Gray', 4),
   (15, 'Eric Rampling', 4),
   (16, 'Piers Paige', 7),
   (17, 'Ryan Henderson', 7),
   (18, 'Frank Tucker', 8),
   (19, 'Nathan Ferguson', 8),
   (20, 'Kevin Rampling', 8);



WITH RECURSIVE subordinates AS (
   SELECT
      employee_id,
      manager_id,
      full_name
   FROM
      employees
   WHERE
      employee_id = 2
   UNION
      SELECT
         e.employee_id,
         e.manager_id,
         e.full_name
      FROM
         employees e
      INNER JOIN subordinates s ON s.employee_id = e.manager_id
) SELECT
   *
FROM
   subordinates;





SELECT n1.id as p, n2.id as c, n2.smiles
FROM nonisomol n1
JOIN edge e ON n1.id = e.parent_id
JOIN nonisomol n2 ON e.child_id = n2.id
WHERE n1.id = 1;

DROP VIEW v_edge;

CREATE OR REPLACE VIEW v_edge AS 
SELECT e.id, e.parent_id, e.child_id,
 p.smiles AS parent_smiles, c.smiles AS child_smiles, 
 p.hac - c.hac AS hac,
 p.rac - c.rac AS rac,
 e.label
FROM edge e
INNER JOIN nonisomol p ON e.parent_id = p.id
INNER JOIN nonisomol c ON e.child_id = c.id;

# parent_id =1 : CN(N)C(=O)c1cnc(Cl)cn1

WITH RECURSIVE fragments AS (
 SELECT parent_id, child_id, parent_smiles, child_smiles, hac, rac, label
 FROM v_edge 
 WHERE parent_id IN (1,2,3,4,5,6,7,8,9,10)
 UNION
   SELECT c.parent_id, c.child_id, c.parent_smiles, c.child_smiles, c.hac, c.rac, c.label
   FROM v_edge c
   INNER JOIN fragments p ON c.parent_id = p.child_id
) SELECT * FROM fragments;




SELECT n.id, n.smiles FROM nonisomol n
INNER JOIN mol_source ms1 ON ms1.nonisomol_id = n.id
INNER JOIN isomol iso ON iso.noniso_id = n.id
INNER JOIN mol_source ms2 ON ms2.isomol_id = iso.id
WHERE ms1.source_id = 2 OR ms2.source_id = 2
LIMIT 10;


SELECT n.id, n.smiles FROM nonisomol n
INNER JOIN mol_source ms ON ms.nonisomol_id = n.id
WHERE ms.source_id = 2
UNION
SELECT n.id, n.smiles FROM nonisomol n
INNER JOIN isomol iso ON iso.noniso_id = n.id
INNER JOIN mol_source ms ON ms.isomol_id = iso.id
WHERE ms.source_id = 2
LIMIT 10;


WITH RECURSIVE fragments AS (
 SELECT parent_id, child_id, parent_smiles, child_smiles, hac, rac, label
 FROM v_edge 
 WHERE parent_id IN (
  SELECT n.id FROM nonisomol n
   INNER JOIN mol_source ms ON ms.nonisomol_id = n.id
   WHERE ms.source_id = 2
  UNION
  SELECT n.id FROM nonisomol n
   INNER JOIN isomol iso ON iso.noniso_id = n.id
   INNER JOIN mol_source ms ON ms.isomol_id = iso.id
   WHERE ms.source_id = 2
  LIMIT 1)
 UNION
   SELECT c.parent_id, c.child_id, c.parent_smiles, c.child_smiles, c.hac, c.rac, c.label
   FROM v_edge c
   INNER JOIN fragments p ON c.parent_id = p.child_id
) SELECT * FROM fragments;

select n.id, n.smiles from mol_source ms join nonisomol n on n.id = ms.nonisomol_id where ms.source_id = 6
union
select n.id, n.smiles from mol_source ms join isomol i on i.id = ms.isomol_id join nonisomol n on i.nonisomol_id = n.id where ms.source_id = 6;

# export to CSV
\COPY (SELECT parent_smiles, child_smiles, label, 'FRAG' FROM v_edge) TO '/tmp/edges.csv' DELIMITER ',' CSV;
 

DROP TABLE i_edge;
DROP TABLE i_node;

CREATE TABLE i_node (
  smiles TEXT,
  present BOOLEAN,
  hac SMALLINT,
  rac SMALLINT,
  scaff TEXT,
  x TEXT,
  labels TEXT
);

CREATE TABLE i_edge (
  p_smiles TEXT,
  c_smiles TEXT,
  present BOOLEAN,
  data TEXT,
  labels TEXT
);
,
  CONSTRAINT fk_i_edge_p FOREIGN KEY (p_smiles) REFERENCES i_node(smiles) ON DELETE CASCADE,
  CONSTRAINT fk_i_edge_c FOREIGN KEY (c_smiles) REFERENCES i_node(smiles) ON DELETE CASCADE
);

CREATE INDEX ix_i_edge_p_smiles on i_edge using hash (p_smiles);
CREATE INDEX ix_i_edge_c_smiles on i_edge using hash (c_smiles);


ALTER TABLE edge ADD CONSTRAINT uq_edge UNIQUE (parent_id, child_id, label);

ALTER TABLE edge DROP CONSTRAINT uq_edge;

\COPY i_node(smiles, hac, rac, scaff, x, labels) FROM '/home/timbo/data/xchem/standardised/nodes.csv' DELIMITER ',' CSV;
\COPY i_edge(p_smiles, c_smiles, data, labels) FROM '/home/timbo/data/xchem/standardised/edges.csv' DELIMITER ',' CSV;

\COPY i_node(smiles, hac, rac, scaff, x, labels) FROM '/home/timbo/tmp/chemspace/10000x0/nodes-unique.csv' DELIMITER ',' CSV;

\COPY i_edge(p_smiles, c_smiles, data, labels) FROM '/home/timbo/tmp/chemspace/10000x0/edges-unique.csv' DELIMITER ',' CSV;


UPDATE i_node i SET present = TRUE WHERE EXISTS 
  (SELECT 1 FROM nonisomol n WHERE n.smiles = i.smiles);
SELECT present, count(*) FROM i_node GROUP BY present;

INSERT INTO nonisomol (smiles, hac, rac)
  SELECT smiles, hac, rac FROM i_node 
  WHERE present IS NULL;


UPDATE i_edge i SET present = TRUE WHERE EXISTS
  (SELECT 1 FROM edge e 
  JOIN nonisomol np ON np.id = e.parent_id 
  JOIN nonisomol nc ON nc.id = e.child_id 
  WHERE i.p_smiles = np.smiles AND i.c_smiles = nc.smiles);
SELECT present, count(*) FROM i_edge GROUP BY present;

INSERT INTO edge (parent_id, child_id, label)
SELECT np.id, nc.id, i.data FROM i_edge i
  JOIN nonisomol np ON np.smiles = i.p_smiles 
  JOIN nonisomol nc ON nc.smiles = i.c_smiles
  WHERE i.present IS NULL;

INSERT INTO edge (parent_id, child_id, label)
SELECT np.id, nc.id, i.data FROM i_edge i
  JOIN nonisomol np ON np.smiles = i.p_smiles 
  JOIN nonisomol nc ON nc.smiles = i.c_smiles
ON CONFLICT ON CONSTRAINT uq_edge DO NOTHING;

\COPY (SELECT np.smiles, nc.smiles, e.label, 'FRAG' FROM edge e JOIN nonisomol np ON e.parent_id = np.id JOIN nonisomol nc ON e.child_id = nc.id) TO '/tmp/edges.csv' DELIMITER ',' CSV;

-- loading nonisomol, isomol and mol_source for
-- ~/data/chemspace/2019-12/00/standardised/standardised-compounds-no-header.tab

DROP TABLE i_mols;
CREATE TABLE i_mols (
  osmiles TEXT,
  isosmiles TEXT,
  nonisosmiles TEXT,
  hac SMALLINT,
  cmpd_id TEXT,
  price INTEGER,
  isomol_id INTEGER,
  nonisomol_id INTEGER
);

\COPY i_mols(osmiles, isosmiles, nonisosmiles, hac, cmpd_id, price) FROM '/home/timbo/data/chemspace/2019-12/00/standardised/standardised-compounds-no-header.tab';


INSERT INTO nonisomol (smiles, hac)
  SELECT nonisosmiles, hac from i_mols
  ON CONFLICT ON CONSTRAINT nonisomol_smiles_key DO NOTHING;

UPDATE i_mols i SET nonisomol_id = n.id
  FROM nonisomol n
    WHERE n.smiles = i.nonisosmiles;

INSERT INTO isomol (smiles, nonisomol_id)
  SELECT isosmiles, nonisomol_id from i_mols
  WHERE isosmiles != nonisosmiles
  ON CONFLICT ON CONSTRAINT isomol_smiles_key DO NOTHING;

UPDATE i_mols i SET isomol_id = n.id
  FROM isomol n
    WHERE n.smiles = i.isosmiles
    AND i.isosmiles != i.nonisosmiles;

INSERT INTO mol_input (name, started_date, source_id)
 VALUES ('chemspace-bb_00.tab', now(), 3);

INSERT INTO mol_source (smiles, code, source_id, input_id, nonisomol_id, isomol_id)
  (SELECT osmiles, cmpd_id, 3, 14, nonisomol_id, isomol_id FROM i_mols);



------------------------

SELECT count(*) FROM nonisomol n WHERE EXISTS
  (SELECT 1 FROM mol_source ms WHERE ms.source_id = 6
    AND ms.nonisomol_id = n.id);

SELECT count(*) FROM nonisomol n WHERE EXISTS
  (SELECT 1 FROM mol_source ms 
    JOIN isomol i ON i.id = ms.isomol_id
    WHERE ms.source_id = 6
      AND ms.isomol_id = i.id
      AND i.nonisomol_id = n.id);

SELECT count(*) FROM nonisomol n WHERE EXISTS
  (SELECT 1 FROM mol_source ms 
    JOIN isomol i ON i.id = ms.isomol_id
    WHERE ms.source_id = 3
      AND (ms.nonisomol_id = n.id
        OR (ms.isomol_id = i.id AND i.nonisomol_id = n.id)
      )
  );

SELECT count(*) FROM nonisomol n WHERE EXISTS
  (SELECT 1 FROM mol_source ms
      WHERE ms.source_id = 6
        AND ms.nonisomol_id = n.id
    UNION ALL
    SELECT 1 FROM mol_source ms 
      JOIN isomol i ON i.id = ms.isomol_id
      WHERE ms.source_id = 6
        AND ms.isomol_id = i.id
        AND i.nonisomol_id = n.id);


SELECT COUNT(*) FROM nonisomol n WHERE EXISTS
  (SELECT 1 FROM mol_source ms WHERE ms.source_id = 6
    AND ms.nonisomol_id = n.id)
UNION
SELECT COUNT(*) FROM nonisomol n WHERE EXISTS
  (SELECT 1 FROM mol_source ms 
    JOIN isomol i ON i.id = ms.isomol_id
    WHERE ms.source_id = 6
      AND ms.isomol_id = i.id
      AND i.nonisomol_id = n.id);

SELECT COUNT(*) FROM nonisomol n WHERE EXISTS
  (SELECT 1 FROM mol_source ms WHERE ms.source_id = 6
    AND ms.nonisomol_id = n.id
  );


SELECT count(*) FROM nonisomol n WHERE EXISTS
  (SELECT 1 FROM mol_source ms 
    JOIN isomol i ON i.id = ms.isomol_id
    WHERE ms.source_id = 6
      AND ((ms.isomol_id = i.id AND i.nonisomol_id = n.id)
      OR ms.nonisomol_id = n.id));


SELECT count(*) FROM nonisomol n WHERE EXISTS
  (SELECT 1 FROM mol_source ms WHERE ms.source_id = 6
    AND ms.nonisomol_id = n.id) OR EXISTS
  (SELECT 1 FROM mol_source ms 
    JOIN isomol i ON i.id = ms.isomol_id
    WHERE ms.source_id = 6
      AND ms.isomol_id = i.id
      AND i.nonisomol_id = n.id);


\COPY (SELECT n.smiles FROM nonisomol n WHERE NOT EXISTS (SELECT 1 FROM edge e WHERE e.parent_id = n.id)) TO '/data/nonisomol.smi';


