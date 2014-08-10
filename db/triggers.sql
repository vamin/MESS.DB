CREATE TRIGGER property
BEFORE INSERT ON molecule_method_property
FOR EACH ROW
BEGIN
    INSERT OR IGNORE INTO property (name, description, format)
    VALUES (NEW.name, NEW.description, NEW.format);
END;