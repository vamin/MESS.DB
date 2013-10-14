SELECT mp.inchikey, theory_level.name, method.name, mp.result
FROM
	molecule_property mp
	JOIN property ON property.property_id = mp.property_id
	JOIN method ON method.method_id = mp.method_id
    JOIN theory_level ON theory_level.theory_id = method.theory_id
WHERE
	property.name = "logp" AND
	mp.result < 0;
