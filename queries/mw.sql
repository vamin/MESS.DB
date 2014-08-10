SELECT m.inchikey, m.smiles, mmp.result as mw
FROM molecule m
JOIN molecule_method_property mmp ON mmp.inchikey=m.inchikey
JOIN property p ON p.property_id=mmp.property_id
WHERE p.name='MW' AND mmp.result BETWEEN 0 AND 10000;
