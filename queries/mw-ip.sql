SELECT m.inchikey, 
(case when p.name = 'MW' then mmp.result end) mw, 
(case when p.name = 'IONIZATION POTENTIAL' then mmp.result end) ip 
FROM molecule m
JOIN molecule_method_property mmp ON mmp.inchikey=m.inchikey
JOIN property p ON p.property_id=mmp.property_id
WHERE
p.name = 'MW'
OR p.name = 'IONIZATION POTENTIAL'
GROUP BY m.inchikey
