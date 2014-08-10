SELECT m.inchikey, mmp.result as ip
FROM molecule m
JOIN molecule_method_property mmp ON mmp.inchikey=m.inchikey
JOIN property p ON p.property_id=mmp.property_id
WHERE p.name='IONIZATION POTENTIAL';
