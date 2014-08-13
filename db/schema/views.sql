--DROP VIEW molecule_method_property_denorm;
CREATE VIEW IF NOT EXISTS molecule_method_property_denorm AS
  SELECT mpp.inchikey,
         mpp.method_path_id,
         p.name,
         p.description,
         p.format,
         mpp.units,
         mpp.result
  FROM molecule_method_property AS mpp
  JOIN property AS p ON mpp.property_id = p.property_id;
